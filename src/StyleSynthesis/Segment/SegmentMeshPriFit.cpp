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

#include "SegmentMeshPriFit.h"

#include <iostream>
#include <unordered_set>

#include "Abstraction/AbstractionSubvolumeTypes.h"

#include "Mesh/MeshUtil.h"
#include "Sample/SampleUtil.h"
#include "Segment/SegmentUtil.h"

#include "Utility/PlyExporter.h"

#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;
using namespace AbstractionSubvolumeNamespace;

#define OUTPUT_PROGRESS

SegmentMeshPriFit::SegmentMeshPriFit() {

}

SegmentMeshPriFit::~SegmentMeshPriFit() {
}

bool SegmentMeshPriFit::loadData(TTriangleMesh &mesh, vector<vector<int>> &patches) {

	mMesh = mesh;
	mPatches = patches;

	mNumVertices = mMesh.amount;
	mNumFaces = (int)mesh.indices.size();
	mNumPatches = (int)patches.size();

	return true;
}

bool SegmentMeshPriFit::runSegmentation() {

	if (!initSegmentation()) return false;
	if (!iterativeMerging()) return false;
	if (!cleanupGroups()) return false;

	return true;
}

bool SegmentMeshPriFit::initSegmentation() {

#ifdef OUTPUT_PROGRESS
	cout << "Initializing segmentation..." << endl;
#endif

	// initialize each patch as individual group

	mPatchGroups.assign(mNumPatches, set<int>());
	mPatchMap.resize(mNumPatches);
	for (int groupID = 0; groupID < mNumPatches; groupID++) {
		mPatchGroups[groupID].insert(groupID);
		mPatchMap[groupID] = groupID;
	}

	// build vertex patch set

	mVertexPatches.assign(mNumVertices, set<int>());
	for (int patchID = 0; patchID < mNumPatches; patchID++) {
		for (int faceID : mPatches[patchID]) {
			vec3i faceIdx = mMesh.indices[faceID];
			for (int k = 0; k < 3; k++) {
				mVertexPatches[faceIdx[k]].insert(patchID);
			}
		}
	}

	// extract vertex set for each patch

	mPatchPoints.resize(mNumPatches);
	if (true) {
		for (int patchID = 0; patchID < mNumPatches; patchID++) {
			unordered_set<int> vertexSet;
			for (int faceID : mPatches[patchID]) {
				vec3i idx = mMesh.indices[faceID];
				for (int k = 0; k < 3; k++) {
					vertexSet.insert(idx[k]);
				}
			}

			TPointSet &points = mPatchPoints[patchID];
			points.positions.clear();
			points.normals.clear();
			for (int vertID : vertexSet) {
				points.positions.push_back(mMesh.positions[vertID]);
				points.normals.push_back(mMesh.normals[vertID]);
			}
			points.amount = (int)points.positions.size();
		}
	}

	return true;
}

bool SegmentMeshPriFit::iterativeMerging() {

#ifdef OUTPUT_PROGRESS
	cout << "Running iterative merging..." << endl;
#endif

	int numIterations = 1000000; // UNDONE: param primitive fitting iterations

	int numPrimTypes = TPrimitive::NUM_TYPES;

	int sampleCount = 0;
	int primCount = 0;
	int fitCount = 0;
	int refitCount = 0;
	int moreCount = 0;
	int mergeCount = 0;

	PlyExporter primCollector;
	
	for (int iterID = 0; iterID < numIterations; iterID++) {

		vec4i sampleVertices;
		set<int> sampleGroups;
		if (!samplePrimitiveVertices(sampleVertices, sampleGroups)) return error("sample vertices");

		for (int typeID = 0; typeID < numPrimTypes; typeID++) {

			// fit primitive
			TPrimitive *prim;
			switch (typeID) {
			case TPrimitive::TYPE_REC_PRISM:
				prim = new TRecPrism(); break;
			case TPrimitive::TYPE_TRI_PRISM:
				prim = new TTriPrism(); break;
			case TPrimitive::TYPE_CYLINDER:
				prim = new TCylinder(); break;
			case TPrimitive::TYPE_TRUNC_CONE:
				prim = new TTruncCone(); break;
			}
			
			if (!prim->fit(&mMesh, sampleVertices)) continue;
			sampleCount++;

			// check primitive
			if (!checkPrimitive(prim)) continue;
			primCount++;

			// check fitting
			bool validFlag = true;
			double allFittingDist = 0;
			TPointSet allFittingPoints;
			if (!checkFitting(prim, sampleGroups, allFittingPoints, allFittingDist, validFlag)) return error("checking fitting");
			if (!validFlag) continue;
			fitCount++;

			// refit primitive
			if (!refitPrimitive(prim, sampleGroups, allFittingPoints, allFittingDist)) return error("refitting");
			refitCount++;

			// check whether more patches can be merged
			if (!checkMerging(prim, sampleGroups)) return error("checking merging");
			moreCount++;

			// merge all groups
			if (!mergeGroups(sampleGroups)) return error("merging group");
			mergeCount++;

			TTriangleMesh primMesh;
			if (!prim->tesselate(primMesh)) return error("tesselation");
			if (!primCollector.addMesh(&primMesh.indices, &primMesh.positions, &primMesh.normals)) return false;

			delete prim;
		}

#ifdef OUTPUT_PROGRESS
		if (iterID % (numIterations/100) == 0) {
			cout << "\rIteration " << iterID << ": ";
			cout << sampleCount << ", " << primCount << ", ";
			cout << fitCount << ", " << refitCount << ", ";
			cout << moreCount << ", " << mergeCount;
			cout << "         ";
		}
#endif
	}
#ifdef OUTPUT_PROGRESS
	cout << endl;
#endif

	mPrimitiveMesh.positions = primCollector.mVertices;
	mPrimitiveMesh.normals = primCollector.mNormals;
	mPrimitiveMesh.indices = primCollector.mFaceIndices;
	mPrimitiveMesh.amount = (int)mPrimitiveMesh.positions.size();
	//if (!visualizePrimitive("prim.ply")) return false;

	return true;
}

bool SegmentMeshPriFit::samplePrimitiveVertices(vec4i &sample, set<int> &groups) {

	bool valid = false;
	while (!valid) {
		set<int> coveredPatches;
		for (int k = 0; k < 4; k++) {
			sample[k] = (int)(SampleUtil::getPreciseRandomNumber() * mNumVertices);
			auto &patches = mVertexPatches[sample[k]];
			coveredPatches.insert(patches.begin(), patches.end());
		}
		groups.clear();
		for (int patchID : coveredPatches) {
			groups.insert(mPatchMap[patchID]);
		}
		if ((int)groups.size() >= 2) {
			valid = true;
		}
	}

	return true;
}

bool SegmentMeshPriFit::checkPrimitive(TPrimitive *primitive) {

	float epsMin = 0.02f; // UNDONE: param primitive eps
	float epsMax = 2.0f;
	if (primitive->mType == TPrimitive::TYPE_TRI_PRISM) {
		auto ptr = (TTriPrism*)primitive;
		float h = (ptr->mCorner[1] - ptr->mCorner[0]).length();
		if (h < epsMin || h > epsMax) return false;
		vec3 c[] = { ptr->mCorner[0], ptr->mCorner[2], ptr->mCorner[3] };
		vec3 v[] = { c[1] - c[0], c[2] - c[1], c[0] - c[2] };
		for (int k = 0; k < 3; k++) {
			float l = v[k].length();
			if (l < epsMin || l > epsMax) return false;
			v[k] /= l;
		}
		float cosA = cos(cml::rad(10.0f));
		for (int k = 0; k < 3; k++) {
			if (fabs(cml::dot(v[k], v[(k + 1) % 3])) > cosA) return false;
		}
	}
	if (primitive->mType == TPrimitive::TYPE_CYLINDER) {
		auto ptr = (TCylinder*)primitive;
		if (ptr->mRadius < epsMin || ptr->mRadius > epsMax) return false;
		float d = (ptr->mCenter[1] - ptr->mCenter[0]).length();
		if (d < epsMin || d > epsMax) return false;
	}
	if (primitive->mType == TPrimitive::TYPE_TRUNC_CONE) {
		auto ptr = (TTruncCone*)primitive;
		if (ptr->mAngle > cml::rad(60.0f)) return false;
		float r = (ptr->mCenter[2] - ptr->mCenter[0]).length() * tan(ptr->mAngle);
		float d = (ptr->mCenter[2] - ptr->mCenter[1]).length();
		if (r < epsMin || d < epsMin || r > epsMax || d > epsMax) return false;
	}

	return true;
}

bool SegmentMeshPriFit::checkFitting(
	TPrimitive *inPrimitive, set<int> &inGroups,
	TPointSet &outFittingPoints, double &outFittingDistance,
	bool &outValidFlag)
{

	double fittingEps = 0.02f; // UNDONE: param fitting distance threshold

	outValidFlag = true;
	outFittingDistance = 0;
	bool gatherPoints = false;
	if (outFittingPoints.amount == 0) gatherPoints = true;
	if (gatherPoints) {
		outFittingPoints.positions.clear();
		outFittingPoints.normals.clear();
	}
	
	for (int groupID : inGroups) {
		for (int patchID : mPatchGroups[groupID]) {
			double dist;
			TPointSet &patchPoints = mPatchPoints[patchID];
			if (!computeFittingDistance(inPrimitive, patchPoints, dist)) return false;
			if (dist > fittingEps) {
				outValidFlag = false;
				return true;
			}
			outFittingDistance = max(outFittingDistance, dist);
			if (gatherPoints) {
				outFittingPoints.positions.insert(outFittingPoints.positions.end(), patchPoints.positions.begin(), patchPoints.positions.end());
				outFittingPoints.normals.insert(outFittingPoints.normals.end(), patchPoints.normals.begin(), patchPoints.normals.end());
			}
		}
	}
	if (gatherPoints) {
		outFittingPoints.amount = (int)outFittingPoints.positions.size();
	}

	return true;
}

bool SegmentMeshPriFit::computeFittingDistance(TPrimitive *primitive, TPointSet &points, double &distance) {

	double totalDist = 0;
	for (vec3 &p : points.positions) {
		totalDist += primitive->distance(p);
	}
	distance = totalDist / points.amount;

	return true;
}

bool SegmentMeshPriFit::refitPrimitive(
	TPrimitive *&primitive, set<int> &groups,
	TPointSet &fittingPoints, double &fittingDistance)
{

	for (int refitID = 0; refitID < 6; refitID++) {
		TPrimitive *newPrim = 0;
		switch (primitive->mType) {
		case TPrimitive::TYPE_REC_PRISM:
			newPrim = new TRecPrism(*(TRecPrism*)primitive);
			break;
		case TPrimitive::TYPE_TRI_PRISM:
			newPrim = new TTriPrism(*(TTriPrism*)primitive);
			break;
		case TPrimitive::TYPE_CYLINDER:
			newPrim = new TCylinder(*(TCylinder*)primitive);
			break;
		case TPrimitive::TYPE_TRUNC_CONE:
			newPrim = new TTruncCone(*(TTruncCone*)primitive);
			break;
		default:
			return error("Unknown primitive");
		}
		if (newPrim->refit(fittingPoints, refitID)) {
			double newFittingDistance;
			bool newValidFlag;
			if (!checkFitting(newPrim, groups, fittingPoints, newFittingDistance, newValidFlag)) return false;
			if (newValidFlag && newFittingDistance < fittingDistance * 1.1) { // relax threshold a bit to promote expanding primitives
				delete primitive;
				primitive = newPrim;
				fittingDistance = newFittingDistance;
			} else {
				delete newPrim;
			}
		} else {
			delete newPrim;
		}
	}

	return true;
}

bool SegmentMeshPriFit::checkMerging(TPrimitive *primitive, set<int> &groups) {

	double mergingEps = 0.02f; // UNDONE: param merging distance threshold

	for (int groupID = 0; groupID < mNumPatches; groupID++) {
		if (groups.find(groupID) != groups.end()) continue;
		bool mergeFlag = true;
		for (int patchID : mPatchGroups[groupID]) {
			double dist;
			TPointSet &patchPoints = mPatchPoints[patchID];
			if (!computeFittingDistance(primitive, patchPoints, dist)) return false;
			if (dist > mergingEps) {
				mergeFlag = false;
				break;
			}
		}
		if (mergeFlag) {
			groups.insert(groupID);
		}
	}

	return true;
}

bool SegmentMeshPriFit::mergeGroups(set<int> &groups) {

	if (groups.empty()) return false;

	// get new group ID (smallest group ID)
	int newGroupID = INT_MAX;
	for (int groupID : groups) {
		newGroupID = min(newGroupID, groupID);
	}

	// update patch map
	for (int patchID = 0; patchID < mNumPatches; patchID++) {
		int oldGroupID = mPatchMap[patchID];
		if (groups.find(oldGroupID) != groups.end()) {
			mPatchMap[patchID] = newGroupID;
		}
	}

	// update groups
	for (int groupID : groups) {
		if (groupID == newGroupID) continue;
		set<int> &oldGroup = mPatchGroups[groupID];
		mPatchGroups[newGroupID].insert(oldGroup.begin(), oldGroup.end());
		oldGroup.clear();
	}

	return true;
}

bool SegmentMeshPriFit::cleanupGroups() {

#ifdef OUTPUT_PROGRESS
	cout << "Cleaning up segmentation groups" << endl;
#endif

	// remove empty groups

	vector<int> newGroupMap(mNumPatches, -1);
	vector<set<int>> newGroup(0);

	int newGroupID = 0;
	for (int oldGroupID = 0; oldGroupID < mNumPatches; oldGroupID++) {
		if (mPatchGroups[oldGroupID].empty()) continue;
		newGroupMap[oldGroupID] = newGroupID;
		newGroup.push_back(mPatchGroups[oldGroupID]);
		newGroupID++;
	}
	mPatchGroups.swap(newGroup);

	for (int patchID = 0; patchID < mNumPatches; patchID++) {
		int oldGroupID = mPatchMap[patchID];
		mPatchMap[patchID] = newGroupMap[oldGroupID];
	}

#ifdef OUTPUT_PROGRESS
	cout << "Extracted " << mPatchGroups.size() << " groups" << endl;
#endif

	return true;
}

bool SegmentMeshPriFit::exportSegmentation(vector<vector<int>> &segments) {

	int numGroups = (int)mPatchGroups.size();
	segments.resize(numGroups);
	for (int groupID = 0; groupID < numGroups; groupID++) {
		segments[groupID].clear();
		for (int patchID : mPatchGroups[groupID]) {
			for (int faceID : mPatches[patchID]) {
				segments[groupID].push_back(faceID);
			}
		}
	}

	return true;
}

bool SegmentMeshPriFit::visualizeSegmentation(string fileName) {

	vector<vec3i> faceColors(mNumFaces, vec3i(127, 127, 127));

	int numGroups = (int)mPatchGroups.size();
	for (int groupID = 0; groupID < numGroups; groupID++) {
		vec3i color = SegmentUtil::colorMapping(groupID);
		for (int patchID : mPatchGroups[groupID]) {
			for (int faceID : mPatches[patchID]) {
				faceColors[faceID] = color;
			}
		}
	}

	PlyExporter pe;
	if (!pe.addMesh(&mMesh.positions, &mMesh.normals, &mMesh.indices, &faceColors)) return false;
	if (!pe.output(fileName)) return false;

	return true;
}

bool SegmentMeshPriFit::visualizePrimitive(string fileName) {

	if (!MeshUtil::saveMesh(fileName, mPrimitiveMesh)) return false;

	return true;
}