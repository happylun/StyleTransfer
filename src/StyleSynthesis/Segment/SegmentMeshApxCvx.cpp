//============================================================================
//
// This file is part of the Style Transfer project.
//
// Copyright (c) 2016
// -Zhaoliang Lun, Evangelos Kalogerakis (authors of the code) / UMass-Amherst
//
// This is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this software.  If not, see <http://www.gnu.org/licenses/>.
//
//============================================================================

#include "SegmentMeshApxCvx.h"

#include <queue>

#include "Utility/PlyExporter.h"

#include "Mesh/MeshUtil.h"
#include "Segment/SegmentUtil.h"

#include "Sample/SampleSimplePoissonDisk.h"

#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

#define OUTPUT_PROGRESS

SegmentMeshApxCvx::SegmentMeshApxCvx() {

}

SegmentMeshApxCvx::~SegmentMeshApxCvx() {
}

bool SegmentMeshApxCvx::initSegmentation(TTriangleMesh &mesh) {

	vector<int> tmpIndices;
	if (!MeshUtil::weldMeshFaces(mesh, mMesh, tmpIndices)) return false;

	if (!extractPatch()) return false;
	if (!initPatchGraph()) return false;
	if (!computePatchVisibility()) return false;
	if (!processConvexNeighbors()) return false;

	return true;
}

bool SegmentMeshApxCvx::extractPatch() {

#ifdef OUTPUT_PROGRESS
	cout << "Extracting patch..." << endl;
#endif

	int numFaces = (int)mMesh.indices.size();
	int numVertices = mMesh.amount;

	// build face graph

	if (!MeshUtil::buildFaceGraph(mMesh, mFaceGraph)) return false;

	// compute face normal

	vector<vec3d> faceNormals(numFaces);
	for (int faceID = 0; faceID < numFaces; faceID++) {
		vec3i faceIdx = mMesh.indices[faceID];
		vec3d facePos[3];
		for (int k = 0; k < 3; k++) facePos[k] = mMesh.positions[faceIdx[k]];
		faceNormals[faceID] = cml::normalize(cml::cross(facePos[1] - facePos[0], facePos[2] - facePos[0]));
	}

	// extract patches

	mPatches.clear();

	vector<bool> visitedFaceFlags(numFaces, false);
	for (int seedFaceID = 0; seedFaceID < numFaces; seedFaceID++) {
		if (visitedFaceFlags[seedFaceID]) continue;
		visitedFaceFlags[seedFaceID] = true;

		// extract planar patch (BFS)

		vector<int> patch;
		if (true) {
			vector<int> queue(1, seedFaceID);
			int head = 0;
			while (head < (int)queue.size()) {
				int currentFace = queue[head];
				for (int neighborFace : mFaceGraph[currentFace]) {
					if (visitedFaceFlags[neighborFace]) continue;

					double angle = acos(fabs(cml::dot(faceNormals[currentFace], faceNormals[neighborFace])));
					if (angle > cml::rad(20.0)) continue; // UNDONE: param normal angle threshold

					queue.push_back(neighborFace);
					visitedFaceFlags[neighborFace] = true;
				}

				head++;
			}
			patch.swap(queue);
		}

		//mPatches.push_back(patch);
		//continue;

		// further split patch into convex polygons

		vector<vector<int>> patches;
		if (!splitPlanarPatch(patch, patches)) return false;

		mPatches.insert(mPatches.end(), patches.begin(), patches.end());
	}
	mNumPatches = (int)mPatches.size();

	// initialize each patch as individual segment

	mSegments.resize(mNumPatches);
	for (int patchID = 0; patchID < mNumPatches; patchID++) {
		mSegments[patchID].assign(1, patchID);
	}
	mNumSegments = mNumPatches;

	// re-orient each patch

	if (true) {
		vec3 bbMin, bbMax;
		if (!MeshUtil::computeAABB(mMesh, bbMin, bbMax)) return false;
		float eps = (bbMax - bbMin).length() * 1e-4f;

		TKDTree meshTree;
		TKDTreeData meshTreeData;
		if (!MeshUtil::buildKdTree(mMesh, meshTree, meshTreeData)) return false;

		for (int patchID = 0; patchID < mNumPatches; patchID++) {
			bool flipFlag;
			if (!MeshUtil::reorientPatch(meshTree, mMesh, mPatches[patchID], flipFlag, eps)) return false;
			if (flipFlag) {
				for (int faceID : mPatches[patchID]) {
					vec3i &faceIdx = mMesh.indices[faceID];
					swap(faceIdx[1], faceIdx[2]);
				}
			}
		}
	}
	/*
	if (true) {
		vector<vec3i> faceColor(numFaces, vec3i(127, 127, 127));
		for (int patchID = 0; patchID < mNumPatches; patchID++) {
			for (int faceID : mPatches[patchID]) {
				faceColor[faceID] = SegmentUtil::colorMapping(patchID);
			}
		}
		PlyExporter pe;
		if (!pe.addMesh(&mMesh.positions, &mMesh.normals, &mMesh.indices, &faceColor)) return false;
		if (!pe.output("patch.ply")) return false;
		system("pause");
	}
	*/

	return true;
}

bool SegmentMeshApxCvx::splitPlanarPatch(vector<int> &inPatch, vector<vector<int>> &outPatches) {

	int numPatchFaces = (int)inPatch.size();
	int numVertices = mMesh.amount;

	// build face mapping

	map<int, int> patchFaceMap; // face ID => patch face ID
	for (int patchFaceID = 0; patchFaceID < numPatchFaces; patchFaceID++) {
		int faceID = inPatch[patchFaceID];
		patchFaceMap[faceID] = patchFaceID;
	}

	// compute face info

	vector<vec3d> faceAngles(numPatchFaces);
	vector<double> faceAreas(numPatchFaces);
	for (int patchFaceID = 0; patchFaceID < numPatchFaces; patchFaceID++) {
		int faceID = inPatch[patchFaceID];
		vec3i faceIdx = mMesh.indices[faceID];
		vec3d facePos[3];
		for (int k = 0; k < 3; k++) facePos[k] = mMesh.positions[faceIdx[k]];
		for (int k = 0; k < 3; k++) {
			faceAngles[patchFaceID][k] = cml::unsigned_angle(facePos[(k + 1) % 3] - facePos[k], facePos[(k + 2) % 3] - facePos[k]);
		}
		faceAreas[patchFaceID] = cml::cross(facePos[1] - facePos[0], facePos[2] - facePos[0]).length();
	}

	// find boundary vertices

	vector<bool> boundaryVertexFlags(numVertices, false);
	if (true) {
		map<vec2i, int> edgeFlags; // (vertex ID, vertex ID) => # of faces touching the edge
		for (int patchFaceID = 0; patchFaceID < numPatchFaces; patchFaceID++) {
			vec3i faceIdx = mMesh.indices[inPatch[patchFaceID]];
			for (int k = 0; k < 3; k++) {
				vec2i edgeKey = makeKey(faceIdx[k], faceIdx[(k + 1) % 3]);
				auto it = edgeFlags.find(edgeKey);
				if (it == edgeFlags.end()) edgeFlags[edgeKey] = 1;
				else it->second++;
			}
		}
		for (auto &edge : edgeFlags) {
			if (edge.second == 1) {
				boundaryVertexFlags[edge.first[0]] = true;
				boundaryVertexFlags[edge.first[1]] = true;
			}
		}
	}

	// find convex polygon

	outPatches.clear();

	vector<bool> visitedPatchFaceFlags(numPatchFaces, false);
	for (int patchFaceID = 0; patchFaceID < numPatchFaces; patchFaceID++) {
		if (visitedPatchFaceFlags[patchFaceID]) continue;
		visitedPatchFaceFlags[patchFaceID] = true;
		int faceID = inPatch[patchFaceID];

		// initialize vertex sum angle
		vector<double> vertexSumAngle(numVertices, 0);
		for (int k = 0; k < 3; k++) {
			int vi = mMesh.indices[faceID][k];
			if (boundaryVertexFlags[vi]) vertexSumAngle[vi] = faceAngles[patchFaceID][k];
		}

		// extract convex face set (BFS)
		vector<int> queue(1, patchFaceID);
		int head = 0;
		while (head < (int)queue.size()) {
			int currentPatchFace = queue[head];
			int currentFace = inPatch[currentPatchFace];
			vec3i currentFaceIdx = mMesh.indices[currentFace];
			for (int neighborFace : mFaceGraph[currentFace]) {
				auto it = patchFaceMap.find(neighborFace);
				if (it == patchFaceMap.end()) continue;
				int neighborPatchFace = it->second;
				if (visitedPatchFaceFlags[neighborPatchFace]) continue;
				vec3i neighborFaceIdx = mMesh.indices[neighborFace];

				// check sum angle on boundary vertices
				if (true) {
					bool passVertexSumAngleCheck = true;
					vec3i vi(-1, -1, -1);
					vec3d va(0.0, 0.0, 0.0);
					for (int k = 0; k < 3; k++) {
						if (!boundaryVertexFlags[neighborFaceIdx[k]]) continue;
						vi[k] = neighborFaceIdx[k];
						va[k] = vertexSumAngle[vi[k]] + faceAngles[neighborPatchFace][k];
						if (va[k] > cml::rad(200.0)) { // UNDONE: param vertex sum angle threshold
							passVertexSumAngleCheck = false;
						}
					}
					if (!passVertexSumAngleCheck) continue;
					for (int k = 0; k < 3; k++) if (vi[k] >= 0) vertexSumAngle[vi[k]] = va[k];
				}

				queue.push_back(neighborPatchFace);
				visitedPatchFaceFlags[neighborPatchFace] = true;
			}

			head++;
		}

		// export patch
		vector<int> convexPatch;
		double patchArea = 0;
		for (int id : queue) {
			convexPatch.push_back(inPatch[id]);
			patchArea += faceAreas[id];
		}
		if(patchArea > 1e-7) outPatches.push_back(convexPatch);
	}

	return true;
}

bool SegmentMeshApxCvx::initPatchGraph() {

#ifdef OUTPUT_PROGRESS
	cout << "Initializing patch graph..." << endl;
#endif

	int numFaces = (int)mMesh.indices.size();
	vector<int> faceMap(numFaces, -1);
	for (int patchID = 0; patchID < mNumPatches; patchID++) {
		for (int faceID : mPatches[patchID]) faceMap[faceID] = patchID;
	}

	mPatchGraph.assign(mNumPatches, set<int>());
	mPatchAdjacency.resize(mNumPatches, mNumPatches);
	mPatchAdjacency.setZero();
	for (int faceID = 0; faceID < numFaces; faceID++) {
		int patchID = faceMap[faceID];
		if (patchID < 0) continue;
		for (int nbFaceID : mFaceGraph[faceID]) {
			int nbPatchID = faceMap[nbFaceID];
			if (nbPatchID < 0) continue;
			if (patchID != nbPatchID) {
				mPatchGraph[patchID].insert(nbPatchID);
				mPatchGraph[nbPatchID].insert(patchID);
				mPatchAdjacency(patchID, nbPatchID) = 1;
				mPatchAdjacency(nbPatchID, patchID) = 1;
			}
		}
	}

	return true;
}

bool SegmentMeshApxCvx::computePatchVisibility() {

#ifdef OUTPUT_PROGRESS
	cout << "Computing patch visibility..." << endl;
#endif

	mPatchVisibility.resize(mNumPatches, mNumPatches);
	mPatchVisibility.setIdentity();

	// gnerate sample points on patches

	int numSamples = 20; // UNDONE: param number of samples on patch

	vector<TSampleSet> patchSamplePoints(mNumPatches);
	mPatchSize.resize(mNumPatches);
	for (int patchID = 0; patchID < mNumPatches; patchID++) {
		TTriangleMesh patch;
		patch.positions = mMesh.positions;
		patch.normals = mMesh.normals;
		patch.indices.clear();
		for (int faceID : mPatches[patchID]) {
			patch.indices.push_back(mMesh.indices[faceID]);
		}
		patch.amount = (int)patch.positions.size();
		SampleSimplePoissonDisk sspd(&patch);
		if (!sspd.runSampling(numSamples)) return false;
		if (!sspd.exportSample(patchSamplePoints[patchID])) return false;

		if (!MeshUtil::computeFaceArea(patch, mPatchSize[patchID])) return false;
	}

	// compute visibility for each connected components separately

	vector<bool> patchFlag(mNumPatches, false);
	for (int patchID = 0; patchID < mNumPatches; patchID++) {
		if (patchFlag[patchID]) continue;
		patchFlag[patchID] = true;

		// find connected components of the graph

		vector<int> queue(1, patchID);
		int head = 0;
		while (head < (int)queue.size()) {
			int currentID = queue[head];
			for (int neighborID : mPatchGraph[currentID]) {
				if (!patchFlag[neighborID]) {
					patchFlag[neighborID] = true;
					queue.push_back(neighborID);
				}
			}
			head++;
		}

		// build KD tree for this connected component

		TTriangleMesh component;
		component.positions = mMesh.positions;
		component.normals = mMesh.normals;
		component.indices.clear();
		for (int patchID : queue) {
			for (int faceID : mPatches[patchID]) {
				component.indices.push_back(mMesh.indices[faceID]);
			}
		}
		component.amount = (int)component.positions.size();

		TKDTree componentTree;
		TKDTreeData componentTreeData;
		if (!MeshUtil::buildKdTree(component, componentTree, componentTreeData)) return false;

		// compute visibility between all pairs of patches within connected components

		vector<vec2i> patchPairs;
		for (int i = 0; i < (int)queue.size() - 1; i++) {
			for (int j = i + 1; j < (int)queue.size(); j++) {
				patchPairs.push_back(vec2i(queue[i], queue[j]));
			}
		}

#pragma omp parallel for shared(patchPairs, patchSamplePoints, componentTree, componentTreeData)
		for (int pairID = 0; pairID < (int)patchPairs.size(); pairID++) {

			int patchID1 = patchPairs[pairID][0];
			int patchID2 = patchPairs[pairID][1];
			TSampleSet &sampleSet1 = patchSamplePoints[patchID1];
			TSampleSet &sampleSet2 = patchSamplePoints[patchID2];

			double visibleCount = 0;
			double totalCount = 0;
			for (int id1 = 0; id1 < sampleSet1.amount; id1++) {
				vec3 p1 = sampleSet1.positions[id1];
				vec3 n1 = sampleSet1.normals[id1];
				for (int id2 = 0; id2 < sampleSet2.amount; id2++) {
					vec3 p2 = sampleSet2.positions[id2];
					vec3 n2 = sampleSet2.normals[id2];

					totalCount += 1;

					// check ray intersection

					vec3 rayDir = p2 - p1;
					float rayLength = rayDir.length();
					rayDir /= rayLength;
					float cos1 = cml::dot(n1, rayDir);
					float cos2 = cml::dot(n2, rayDir);
					if (fabs(cos1) < 1e-1f && fabs(cos2) < 1e-1f) {
						visibleCount += 0.9; // UNDONE: param co-planar visbility vote
						continue;
					}
					if (cos1 > 0 || cos2 < 0) continue; // ray is outside		

					float rayEps = rayLength * 1e-3f;
					vec3 rayOrigin = p1 + rayDir * rayEps;
					Thea::Ray3 ray(G3D::Vector3(rayOrigin.data()), G3D::Vector3(rayDir.data()));
					float dist = componentTree.rayIntersectionTime(ray);
					if (dist < -0.5f) {
						// not hit? why? anyway, mark it as visible...
						visibleCount += 0.5;
					} else if (dist + rayEps * 2 > rayLength) {
						visibleCount += 1.0;
					} else {
						// something blocked the ray. mark it as invisible
						continue;
					}
				}
			}

			double visibility = visibleCount / totalCount;
			mPatchVisibility(patchID1, patchID2) = visibility;
			mPatchVisibility(patchID2, patchID1) = visibility;
		}
	}

	return true;
}

bool SegmentMeshApxCvx::processConvexNeighbors() {

#ifdef OUTPUT_PROGRESS
	cout << "Processing convex neighbors..." << endl;
#endif

	// extract edge info

	vector<vec2i> allEdges(0); // (vertex ID, vertex ID) : # of all edges
	map<vec2i, int> allEdgeMap; // (vertex ID, vertex ID) => edge ID
	map<int, set<vec2i>> allEdgeFaces; // edge ID => set of (patchID, face ID)
	vector<set<int>> allPatchEdges(mNumPatches, set<int>()); // set of edge ID : # of patches

	for (int patchID = 0; patchID < mNumPatches; patchID++) {
		auto &patch = mPatches[patchID];
		for (int faceID : patch) {
			vec3i faceIdx = mMesh.indices[faceID];
			for (int k = 0; k < 3; k++) {
				vec2i edgeKey = makeKey(faceIdx[k], faceIdx[(k + 1) % 3]);
				auto itEM = allEdgeMap.find(edgeKey);
				int edgeID;
				if (itEM == allEdgeMap.end()) {
					edgeID = (int)allEdges.size();
					allEdgeMap[edgeKey] = edgeID;
					allEdges.push_back(edgeKey);
				} else {
					edgeID = itEM->second;
				}
				allEdgeFaces[edgeID].insert(vec2i(patchID, faceID));
				allPatchEdges[patchID].insert(edgeID);
			}
		}
	}

	// get all neighboring patches

	set<vec2i> patchPairSet;
	for (int patchID = 0; patchID < mNumPatches; patchID++) {
		for (int neighborID : mPatchGraph[patchID]) {
			patchPairSet.insert(makeKey(patchID, neighborID));
		}
	}
	vector<vec2i> patchPairList(patchPairSet.begin(), patchPairSet.end());

	// detect neighboring convex patches

	TKDTree meshTree;
	TKDTreeData meshTreeData;
	if (!MeshUtil::buildKdTree(mMesh, meshTree, meshTreeData)) return false;

	mPatchConvexAdjacency.resize(mNumPatches, mNumPatches);
	mPatchConvexAdjacency.setZero();

	int numPatchPairs = (int)patchPairList.size();
#pragma omp parallel for
	for (int pairID = 0; pairID < numPatchPairs; pairID++) {

		vec2i patchPair = patchPairList[pairID];
		int srcPatchID = patchPair[0];
		int tgtPatchID = patchPair[1];
		auto &srcEdges = allPatchEdges[srcPatchID];
		auto &tgtEdges = allPatchEdges[tgtPatchID];

		// find common edges
		vector<int> commonEdges;
		set_intersection(
			srcEdges.begin(), srcEdges.end(),
			tgtEdges.begin(), tgtEdges.end(),
			back_inserter(commonEdges));

		// find neighboring faces
		vector<vec2i> facePairs(0);
		facePairs.reserve(commonEdges.size());
		for (int edgeID : commonEdges) {
			vector<int> srcFaces(0);
			vector<int> tgtFaces(0);
			set<vec2i> &edgeFaces = allEdgeFaces[edgeID];
			for (vec2i face : edgeFaces) {
				if (face[0] == srcPatchID) srcFaces.push_back(face[1]);
				if (face[0] == tgtPatchID) tgtFaces.push_back(face[1]);
			}
			for (int srcFaceID : srcFaces) {
				for (int tgtFaceID : tgtFaces) {
					facePairs.push_back(vec2i(srcFaceID, tgtFaceID));
				}
			}
		}
		
		// check visibility
		double invisibleWeight = 0;
		double visibleWeight = 0;
		for (vec2i &facePair : facePairs) {
			vec3i srcFaceIdx = mMesh.indices[facePair[0]];
			vec3i tgtFaceIdx = mMesh.indices[facePair[1]];
			vec3d srcFacePos[3], tgtFacePos[3];
			for (int k = 0; k < 3; k++) {
				srcFacePos[k] = mMesh.positions[srcFaceIdx[k]];
				tgtFacePos[k] = mMesh.positions[tgtFaceIdx[k]];
			}
			vec3d srcCenter = (srcFacePos[0] + srcFacePos[1] + srcFacePos[2]) / 3;
			vec3d tgtCenter = (tgtFacePos[0] + tgtFacePos[1] + tgtFacePos[2]) / 3;
			vec3d srcNormal = cml::cross(srcFacePos[1] - srcFacePos[0], srcFacePos[2] - srcFacePos[0]);			
			vec3d tgtNormal = cml::cross(tgtFacePos[1] - tgtFacePos[0], tgtFacePos[2] - tgtFacePos[0]);
			double srcArea = srcNormal.length();
			double tgtArea = tgtNormal.length();
			srcNormal.normalize();
			tgtNormal.normalize();
			double weight = srcArea + tgtArea;
			double angle;
			if (!MeshUtil::computeDihedralAngle(srcCenter, srcNormal, tgtCenter, tgtNormal, angle)) error("compute dihedral angle");
			bool visible = false;
			if (angle > cml::rad(200.0)) { // UNDONE: neighbor convex angle threshold
				// check ray visibility
				vec3 rayDir = tgtCenter - srcCenter;
				float rayLen = rayDir.length();
				if (rayLen) {
					float rayEps = rayLen * 1e-3f;
					rayDir /= rayLen;
					vec3 rayOrigin = srcCenter + rayDir * rayEps;
					Thea::Ray3 ray(G3D::Vector3(rayOrigin.data()), G3D::Vector3(rayDir.data()));
					float dist = meshTree.rayIntersectionTime(ray);
					if (dist < -0.5f || dist + rayEps * 2 > rayLen) visible = true;
				}
			}
			if (visible) visibleWeight += weight;
			else invisibleWeight += weight;
		}

		// determine convex neighbor
		double visibility = visibleWeight / (visibleWeight + invisibleWeight);
		if (visibility > 0.9) {
			mPatchConvexAdjacency(srcPatchID, tgtPatchID) = 1;
			mPatchConvexAdjacency(tgtPatchID, srcPatchID) = 1;
		}
	}

	return true;
}

bool SegmentMeshApxCvx::computeSegmentVisibility(vector<int> &seg1, vector<int> &seg2, double &visibility) {

	// use average visibility

	double totalCount = 0;
	double visibleCount = 0;
	bool isConvexNeighbor = true; // whether all neighboring patches are convex neighbors
	for (int patchID1 : seg1) {
		for (int patchID2 : seg2) {
			double patchVisibility = mPatchVisibility(patchID1, patchID2);
			double rayAmount = (double)(mPatchSize[patchID1] * mPatchSize[patchID2]);
			totalCount += rayAmount;
			visibleCount += patchVisibility * rayAmount;
			if (mPatchAdjacency(patchID1, patchID2) == 1) {
				if (mPatchConvexAdjacency(patchID1, patchID2) == 0) {
					isConvexNeighbor = false;
				}
			}
		}
	}
	visibility = visibleCount / totalCount;
	if (isConvexNeighbor) {
		visibility = max(0.9, visibility);
	}

	return true;
}

bool SegmentMeshApxCvx::runSegmentation(double visibility) {

#ifdef OUTPUT_PROGRESS
	cout << "Running segmentation..." << endl;
#endif

	// get patch label

	vector<int> patchLabel(mNumPatches, -1);
	for (int segID = 0; segID < mNumSegments; segID++) {
		for (int patchID : mSegments[segID]) {
			patchLabel[patchID] = segID;
		}
	}

	// get initial merge pairs

	map<vec2i, double> mergeScore;
	for (int patchID = 0; patchID < mNumPatches; patchID++) {
		int patchSegment = patchLabel[patchID];
		for (int neighborID : mPatchGraph[patchID]) {
			int neighborSegment = patchLabel[neighborID];
			if (patchSegment == neighborSegment) continue;
			vec2i key = makeKey(patchSegment, neighborSegment);
			if (mergeScore.find(key) != mergeScore.end()) continue; // already computed
			double score;
			if (!computeSegmentVisibility(mSegments[patchSegment], mSegments[neighborSegment], score)) return false;
			mergeScore[key] = score;
		}
	}

	// build priority queue

	typedef pair<vec2i, double> TQueuePair;
	auto LowerScore = [](TQueuePair &p1, TQueuePair &p2) { return p1.second < p2.second; };
	typedef priority_queue<TQueuePair, vector<TQueuePair>, decltype(LowerScore)> TQueue;

	TQueue queue(LowerScore);
	for (auto &it : mergeScore) {
		TQueuePair pair;
		pair.first = it.first;
		pair.second = it.second;
		queue.push(pair);
	}

	// merge in order

	vector<bool> segmentFlags(mNumSegments, true);

	while (!queue.empty()) {
		TQueuePair topPair = queue.top();
		if (topPair.second < visibility) break;
		queue.pop();
		int seg1 = topPair.first[0];
		int seg2 = topPair.first[1];
		if (!segmentFlags[seg1] || !segmentFlags[seg2]) continue; // already merged into other segment
		double previousScore = topPair.second;
		auto currentEntry = mergeScore.find(topPair.first);
		if (currentEntry == mergeScore.end()) continue; // no longer available (?!)
		double currentScore = currentEntry->second;
		if (currentScore != previousScore) continue; // deprecated entry

		// merge it
		if (true) {
			mSegments[seg1].insert(mSegments[seg1].end(), mSegments[seg2].begin(), mSegments[seg2].end());
			mSegments[seg2].clear();
			segmentFlags[seg2] = false;
			vector<TQueuePair> newEntries;
			for (auto it = mergeScore.begin(); it != mergeScore.end();) {
				if (it->first == topPair.first) {
					// erase
					mergeScore.erase(it++);
				} else if (it->first[0] == seg1 || it->first[1] == seg1) {
					// update score
					double newScore;
					if (!computeSegmentVisibility(mSegments[it->first[0]], mSegments[it->first[1]], newScore)) return false;
					TQueuePair newPair(it->first, newScore);
					queue.push(newPair);
					it->second = newScore;
					++it;
				} else if (it->first[0] == seg2 || it->first[1] == seg2) {
					// erase old entry & put in new entry
					vec2i newKey;
					if (it->first[0] == seg2) newKey = makeKey(seg1, it->first[1]);
					else newKey = makeKey(seg1, it->first[0]);
					double newScore;
					if (!computeSegmentVisibility(mSegments[newKey[0]], mSegments[newKey[1]], newScore)) return false;
					TQueuePair newPair(newKey, newScore);
					newEntries.push_back(newPair);
					queue.push(newPair);
					mergeScore.erase(it++);
				} else {
					// skip
					++it;
				}
			}
			for (auto &it : newEntries) {
				mergeScore[it.first] = it.second;
			}
		}
	}

	// remove empty segments (those merged into other segments)

	vector<vector<int>> cleanedSegemnts(0);
	for (int segID = 0; segID < mNumSegments; segID++) {
		if (segmentFlags[segID]) cleanedSegemnts.push_back(mSegments[segID]);
	}
	mSegments.swap(cleanedSegemnts);
	mNumSegments = (int)mSegments.size();

	return true;
}

bool SegmentMeshApxCvx::exportMesh(TTriangleMesh &outMesh) {

	outMesh = mMesh;

	return true;
}

bool SegmentMeshApxCvx::exportSegmentation(vector<vector<int>> &outSegments) {

	outSegments.resize(mNumSegments, vector<int>(0));
	for (int segID = 0; segID < mNumSegments; segID++) {
		auto &segment = outSegments[segID];
		for (int patchID : mSegments[segID]) {
			for (int faceID : mPatches[patchID]) {
				segment.push_back(faceID);
			}
		}
	}

	return true;
}

bool SegmentMeshApxCvx::visualizeSegmentation(string fileName) {

	int numFaces = (int)mMesh.indices.size();
	vector<vec3i> faceColors(numFaces, vec3i(127, 127, 127));

	for (int segID = 0; segID < mNumSegments; segID++) {
		vec3i segColor = SegmentUtil::colorMapping(segID);
		for (int patchID : mSegments[segID]) {
			for (int faceID : mPatches[patchID]) {
				faceColors[faceID] = segColor;
			}
		}
	}

	PlyExporter pe;
	if (!pe.addMesh(&mMesh.positions, &mMesh.normals, &mMesh.indices, &faceColors)) return false;
	if (!pe.output(fileName)) return false;

	return true;
}

