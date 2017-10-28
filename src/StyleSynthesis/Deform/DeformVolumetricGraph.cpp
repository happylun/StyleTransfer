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

#include "DeformVolumetricGraph.h"

#include <set>
#include <omp.h>

#include "Library/libiglHelperBegin.h"
#include <igl/tetgen/tetrahedralize.h>
#include "Library/libiglHelperEnd.h"

#include "Mesh/MeshUtil.h"
#include "Sample/SampleUtil.h"

#include "Sample/SamplePoissonDisk.h"

#include "Utility/PlyExporter.h"

#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

#define OUTPUT_VISUALIZATION

DeformVolumetricGraph::DeformVolumetricGraph() {

}

DeformVolumetricGraph::~DeformVolumetricGraph() {

}

bool DeformVolumetricGraph::loadMesh(TTriangleMesh &mesh) {

	mpMesh = &mesh;

	// compute boundign box

	if (!MeshUtil::computeAABB(mesh, mBBMin, mBBMax)) return false;

	// compute average edge length

	double totalLength = 0;
	int totalEdges = 0;
	for (vec3i &face : mesh.indices) {
		for (int k = 0; k < 3; k++) {
			float edgeLength = (mesh.positions[face[k]] - mesh.positions[face[(k + 1) % 3]]).length();
			totalLength += (double)edgeLength;
			totalEdges++;
		}
	}
	if (totalEdges == 0) return error("empty mesh");
	mAverageEdgeLength = totalLength / totalEdges;

	// build Kd tree

	if (!MeshUtil::buildKdTree(mesh, mMeshTree, mMeshTreeData)) return false;

	// add vertices to point cloud

	mPointCloud = mesh.positions;

	return true;
}

bool DeformVolumetricGraph::loadNames(string visualizeFolder) {

	mVisualizeFolder = visualizeFolder;

	return true;
}

bool DeformVolumetricGraph::process() {

	if (!latticeSampling()) return false;
	if (!generateShell()) return false;
	if (!tetrahedralization()) return false;
	if (!extractAlphaShape()) return false;
	if (!simplifyGraph()) return false;
	if (!smoothGraph()) return false;
	if (!washoutInternal()) return false;

	cout << "Done." << endl;

	return true;
}

bool DeformVolumetricGraph::latticeSampling() {

	cout << "Doing lattice sampling" << endl;

	float latticeSize = (float)mAverageEdgeLength; // UNDONE: param lattice size

	// build lattice

	vec2i latticeRange[3]; // (min, max) : {x, y, z}
	int latticeNumber[3]; // amount : {x, y, z}
	for (int k = 0; k < 3; k++) {
		latticeRange[k][0] = (int)(floor(mBBMin[k] / latticeSize));
		latticeRange[k][1] = (int)(ceil(mBBMax[k] / latticeSize));
		latticeNumber[k] = latticeRange[k][1] - latticeRange[k][0] + 1;
	}
	cout << "Lattice: " << latticeNumber[0] << " X " << latticeNumber[1] << " X " << latticeNumber[2] << endl;

	// initialize lattice flag

	vector<vector<vector<bool>>> latticeFlags; // valid : z : y : x
	latticeFlags.resize(latticeNumber[0]);
	for (int x = 0; x < latticeNumber[0]; x++) {
		latticeFlags[x].resize(latticeNumber[1]);
		for (int y = 0; y < latticeNumber[1]; y++) {
			latticeFlags[x][y].assign(latticeNumber[2], false);
		}
	}

	// mark inside samples

	vector<vec3i> sampleList(0);
	for (int x = 0; x < latticeNumber[0]; x++) {
		for (int y = 0; y < latticeNumber[1]; y++) {
			for (int z = 0; z < latticeNumber[2]; z++) {
				sampleList.push_back(vec3i(x, y, z));
			}
		}
	}
	int numSamples = (int)sampleList.size();
#pragma omp parallel for shared(latticeFlags, sampleList)
	for (int sampleID = 0; sampleID < numSamples; sampleID++) {
		vec3i bin = sampleList[sampleID];
		vec3 sample;
		for (int k = 0; k < 3; k++) {
			sample[k] =  (latticeRange[k][0] + bin[k]) * latticeSize;
		}
		if (checkInside(sample)) {
			latticeFlags[bin[0]][bin[1]][bin[2]] = true;
		}
	}

	// insert samples into point cloud

	vector<vec3> sampleSet(0);
	for (int sampleID = 0; sampleID < numSamples; sampleID++) {
		vec3i bin = sampleList[sampleID];
		if (!latticeFlags[bin[0]][bin[1]][bin[2]]) continue;
		vec3 sample;
		for (int k = 0; k < 3; k++) {
			sample[k] = (latticeRange[k][0] + bin[k]) * latticeSize;
		}
		sampleSet.push_back(sample);
	}
	mPointCloud.insert(mPointCloud.end(), sampleSet.begin(), sampleSet.end());

#ifdef OUTPUT_VISUALIZATION
	if (true) {
		// visualize
		PlyExporter pe;
		if (!pe.addPoint(&sampleSet)) return false;
		if (!pe.output(mVisualizeFolder + "lattice.ply")) return false;
	}
#endif

	return true;
}

bool DeformVolumetricGraph::generateShell() {

	cout << "Generating shell" << endl;

	float shellOffset = (float)mAverageEdgeLength; // UNDONE: param shell offset

	if (!MeshUtil::recomputeNormals(*mpMesh)) return false;

	vector<bool> outerFlag(mpMesh->amount, true);
#pragma omp parallel for
	for (int vertID = 0; vertID < mpMesh->amount; vertID++) {
		vec3 point = mpMesh->positions[vertID];
		if (checkInside(point)) outerFlag[vertID] = false;
	}

	vector<vec3> sampleSet(0);
	for (int vertID = 0; vertID < mpMesh->amount; vertID++) {
		if (!outerFlag[vertID]) continue;
		vec3 sample = mpMesh->positions[vertID] + mpMesh->normals[vertID] * shellOffset;
		sampleSet.push_back(sample);
	}
	mPointCloud.insert(mPointCloud.end(), sampleSet.begin(), sampleSet.end());

#ifdef OUTPUT_VISUALIZATION
	if (true) {
		// visualize
		PlyExporter pe;
		if (!pe.addPoint(&sampleSet)) return false;
		if (!pe.output(mVisualizeFolder + "shell.ply")) return false;
	}
#endif

	return true;
}

bool DeformVolumetricGraph::tetrahedralization() {

	cout << "Tetrahedralizing point cloud" << endl;

	int numPoints = (int)mPointCloud.size();
	Eigen::MatrixXd matV(numPoints, 3);
	for (int pointID = 0; pointID < numPoints; pointID++) {
		for (int dim = 0; dim < 3; dim++) {
			matV(pointID, dim) = mPointCloud[pointID][dim];
		}
	}
	Eigen::MatrixXi matNull(0, 0);

	Eigen::MatrixXd matTV;
	Eigen::MatrixXi matTT;
	Eigen::MatrixXi matTF;
	igl::tetgen::tetrahedralize(matV, matNull, "pcQ", matTV, matTT, matTF);

	if (!MeshUtil::mat2tet(mTetraMesh, matTV, matTT)) return false;

#ifdef OUTPUT_VISUALIZATION
	TTriangleMesh visualMesh;
	if (!MeshUtil::tet2mesh(mTetraMesh, visualMesh)) return false;
	if (!MeshUtil::saveMesh(mVisualizeFolder + "tet.ply", visualMesh)) return false;

	TTriangleMesh refMesh = *mpMesh;
	refMesh.positions = mTetraMesh.positions;
	refMesh.positions.resize(mpMesh->amount);
	if (!MeshUtil::saveMesh(mVisualizeFolder + "tet-ref.ply", refMesh)) return false;
#endif

	return true;
}

bool DeformVolumetricGraph::extractAlphaShape() {

	cout << "Extracting alpha shape" << endl;

	int numTets = (int)mTetraMesh.indices.size();

	double maxCircumRadius = mAverageEdgeLength * 2; // UNDONE: param maximum circumradius

	// compute tetrahedron's circumradius & length of longest side
	// ref: http://mathworld.wolfram.com/Tetrahedron.html

	vector<double> tetRadius(numTets);
	vector<double> tetLength(numTets);
#pragma omp parallel for
	for (int tetID = 0; tetID < numTets; tetID++) {
		vec4i idx = mTetraMesh.indices[tetID];

		// get vertices and sides

		vec3 v[4];
		for (int k = 0; k < 4; k++) v[k] = mTetraMesh.positions[idx[k]];
		vec3d a = v[1] - v[0]; double al = a.length();
		vec3d b = v[2] - v[0]; double bl = b.length();
		vec3d c = v[3] - v[0]; double cl = c.length();
		vec3d ap = v[3] - v[2]; double apl = ap.length();
		vec3d bp = v[3] - v[1]; double bpl = bp.length();
		vec3d cp = v[2] - v[1]; double cpl = cp.length();

		// compute circumradius

		double sa = al*apl;
		double sb = bl*bpl;
		double sc = cl*cpl;
		if (sa < sb) swap(sa, sb);
		if (sa < sc) swap(sa, sc);
		if (sb < sc) swap(sb, sc);
		// Kahan's Heron's formula
		double delta = sqrt((sa + (sb + sc)) * (sc - (sa - sb)) * (sc + (sa - sb)) * (sa + (sb - sc))) / 4;
		double volume = fabs(cml::dot(a, cml::cross(b, c))) / 6;
		double radius = delta / (volume * 6);

		// compute length of longest side

		double length = max(al, max(bl, max(cl, max(apl, max(bpl, cpl)))));

		tetRadius[tetID] = radius;
		tetLength[tetID] = length;
	}

	// compute fixed vertex associated tetrahedra counts

	int numFixedPoints = mpMesh->amount;

	vector<int> tetCounts(numFixedPoints, 0);
	for (int tetID = 0; tetID < numTets; tetID++) {
		vec4i idx = mTetraMesh.indices[tetID];
		for (int k = 0; k < 4; k++) {
			if (idx[k] < numFixedPoints) tetCounts[idx[k]]++;
		}
	}
	for (int vertID = 0; vertID < numFixedPoints; vertID++) {
		if (tetCounts[vertID] == 0) return error("found vertex not associated with any tets");
	}

	// sort tetrahedra by decreasing maximum side length

	vector<int> tetOrderList(numTets);
	for (int k = 0; k < numTets; k++) tetOrderList[k] = k;
	sort(tetOrderList.begin(), tetOrderList.end(),
		[&tetLength](int i1, int i2){return tetLength[i1] > tetLength[i2]; });

	// remove tetrahedra one by one

	vector<bool> tetFlags(numTets, true);
	int removedCount = 0;
	int skippedCount = 0;
	for (int orderID = 0; orderID < numTets; orderID++) {
		int tetID = tetOrderList[orderID];
		if (tetRadius[tetID] < maxCircumRadius) continue;
		vec4i tetIdx = mTetraMesh.indices[tetID];
		bool removeFlag = true;
		for (int k = 0; k < 4; k++) {
			if (tetIdx[k] < numFixedPoints && tetCounts[tetIdx[k]] <= 1) {
				removeFlag = false;
				skippedCount++;
				break;
			}
		}
		if (removeFlag) {
			for (int k = 0; k < 4; k++) {
				if (tetIdx[k] < numFixedPoints) tetCounts[tetIdx[k]]--;
			}
			tetFlags[tetID] = false;
			removedCount++;
		}
	}

	cout << "Skipped " << skippedCount << " tetrahedra" << endl;
	cout << "Removed " << removedCount << " tetrahedra" << endl;

	vector<vec4i> newIndices(0);
	for (int tetID = 0; tetID < numTets; tetID++) {
		if (tetFlags[tetID]) newIndices.push_back(mTetraMesh.indices[tetID]);
	}
	mTetraMesh.indices.swap(newIndices);

	if (!MeshUtil::cleanUp(mTetraMesh)) return false;

#ifdef OUTPUT_VISUALIZATION
	TTriangleMesh visualMesh;
	if (!MeshUtil::tet2mesh(mTetraMesh, visualMesh)) return false;
	if (!MeshUtil::saveMesh(mVisualizeFolder + "alpha.ply", visualMesh)) return false;

	TTriangleMesh refMesh = *mpMesh;
	refMesh.positions = mTetraMesh.positions;
	refMesh.positions.resize(mpMesh->amount);
	if (!MeshUtil::saveMesh(mVisualizeFolder + "alpha-ref.ply", refMesh)) return false;
#endif

	return true;
}

bool DeformVolumetricGraph::simplifyGraph() {

	cout << "Simplifying graph" << endl;

	double collapseMinLength = mAverageEdgeLength * 0.5;

	int numFixedPoints = mpMesh->amount;

	// get all edges

	map<vec2i, float> edgeSet; // (id1, id2) => edge length
	for (int tetID = 0; tetID < (int)mTetraMesh.indices.size(); tetID++) {
		vec4i idx = mTetraMesh.indices[tetID];
		for (int k = 0; k < 4; k++) {
			vec2i key = makeKey(idx[k], idx[(k + 1) % 4]);
			auto insertion = edgeSet.insert(make_pair(key, 0.0f));
			auto &iter = insertion.first;
			if (insertion.second) { // newly added
				iter->second = (mTetraMesh.positions[key[0]] - mTetraMesh.positions[key[1]]).length();
			}
		}
	}

	// count tets associated with fixed vertices

	vector<int> tetCounts(numFixedPoints, 0); // # of tets : # of fixed vertices
	for (vec4i idx : mTetraMesh.indices) {
		for (int k = 0; k < 4; k++) {
			if (idx[k] < numFixedPoints) tetCounts[idx[k]]++;
		}
	}
	for (int vertID = 0; vertID < numFixedPoints; vertID++) {
		if (tetCounts[vertID] == 0) return error("found vertex not associated with any tets");
	}

	// sort edges by increasing length

	vector<vec2i> edgeList;
	vector<float> edgeLengths;
	for (auto &edge : edgeSet) {
		edgeList.push_back(edge.first);
		edgeLengths.push_back(edge.second);
	}
	int numEdges = (int)edgeList.size();
	vector<int> edgeOrderList(numEdges);
	for (int k = 0; k < numEdges; k++) edgeOrderList.push_back(k);
	sort(edgeOrderList.begin(), edgeOrderList.end(),
		[&edgeLengths](int i1, int i2){return edgeLengths[i1] < edgeLengths[i2]; });

	// collapse one by one

	vector<bool> tetFlag(mTetraMesh.indices.size(), true);

	vector<int> vertexMap(mTetraMesh.amount);
	for (int k = 0; k < mTetraMesh.amount; k++) vertexMap[k] = k;

	int collapsedCount = 0;
	for (int orderID = 0; orderID < numEdges; orderID++) {
		int edgeID = edgeOrderList[orderID];
		vec2i key = edgeList[edgeID];

		// get real edge key

		key[0] = getVertexMap(vertexMap, key[0]);
		key[1] = getVertexMap(vertexMap, key[1]);
		if (key[0] == key[1]) continue; // already collapsed
		key = makeKey(key);
		if (key[1] < numFixedPoints) continue; // fixed edge

		// check edge length

		vec3 p0 = mTetraMesh.positions[key[0]];
		vec3 p1 = mTetraMesh.positions[key[1]];
		float edgeLength = (p1 - p0).length();
		if (edgeLength > collapseMinLength) break; // no need to collapse any more (not really but anyway...)

		// check fixed vertex collapsed tet count
		// TODO: need a faster implementation...

#ifdef _OPENMP
		int numThreads = omp_get_max_threads();
#else
		int numThreads = 1;
#endif

		vector<bool> ompFlags(numThreads, true); // can collapsed
		vector<map<int, int>> ompTetCounts(numThreads, map<int, int>()); // vert ID => collapsed tets count
		vector<vector<int>> ompTetList(numThreads, vector<int>(0)); // tet ID : # of tets need to be updated

#pragma omp parallel for num_threads(numThreads)
		for (int tetID = 0; tetID < (int)mTetraMesh.indices.size(); tetID++) {

#ifdef _OPENMP
			int threadID = omp_get_thread_num();
#else
			int threadID = 0;
#endif

			if (!ompFlags[threadID]) continue;
			if (!tetFlag[tetID]) continue;
			vec4i tetIdx = mTetraMesh.indices[tetID];
			bool hasFrom = false;
			bool hasTo = false;
			for (int k = 0; k < 4; k++) {
				if (tetIdx[k] == key[0]) hasTo = true;
				if (tetIdx[k] == key[1]) hasFrom = true;
			}
			if (hasFrom) {
				ompTetList[threadID].push_back(tetID);
				if (hasTo) {
					for (int k = 0; k < 4; k++) {
						if (tetIdx[k] >= numFixedPoints) continue;
						auto iter = ompTetCounts[threadID].insert(make_pair(tetIdx[k], 1));
						if (!iter.second) iter.first->second++;
						if (iter.first->second >= tetCounts[tetIdx[k]]) {
							ompFlags[threadID] = false;
							break;
						}
					}
				}
			}
		}

		// merge omp results
		bool canCollapse = true;
		vector<int> updateTetList;
		map<int, int> updateTetCounts;
		for (int threadID = 0; threadID < numThreads; threadID++) {
			if (!ompFlags[threadID]) {
				canCollapse = false;
				break;
			}
			updateTetList.insert(updateTetList.end(), ompTetList[threadID].begin(), ompTetList[threadID].end());
			for (auto it : ompTetCounts[threadID]) {
				auto insertion = updateTetCounts.insert(it);
				if (!insertion.second) {
					insertion.first->second += it.second;
				}
				if (insertion.first->second >= tetCounts[insertion.first->first]) canCollapse = false;
				break;
			}
			if (!canCollapse) break;
		}
		if (!canCollapse) continue;

		// move vertex

		if (key[0] >= numFixedPoints) {
			mTetraMesh.positions[key[0]] = (p0 + p1) / 2;
		}

		// update tet indices

		for (int tetID : updateTetList) {
			vec4i &tetIdx = mTetraMesh.indices[tetID];
			for (int k = 0; k < 4; k++) {
				if (tetIdx[k] == key[0]) tetFlag[tetID] = false;
				if (tetIdx[k] == key[1]) tetIdx[k] = key[0];
			}
		}

		// update tet counts

		for (auto &it : updateTetCounts) {
			tetCounts[it.first] -= it.second;
		}

		// update vertex map
		vertexMap[key[1]] = key[0];

		collapsedCount++;
	}

	cout << "Collapsed " << collapsedCount << " edges" << endl;

	// remove invalid tets

	vector<vec4i> newIndices;
	for (int tetID = 0; tetID < (int)mTetraMesh.indices.size(); tetID++) {
		if (tetFlag[tetID]) newIndices.push_back(mTetraMesh.indices[tetID]);
	}
	mTetraMesh.indices.swap(newIndices);

	if (!MeshUtil::cleanUp(mTetraMesh)) return false;

#ifdef OUTPUT_VISUALIZATION
	TTriangleMesh visualMesh;
	if (!MeshUtil::tet2mesh(mTetraMesh, visualMesh)) return false;
	if (!MeshUtil::saveMesh(mVisualizeFolder + "simplified.ply", visualMesh)) return false;

	TTriangleMesh refMesh = *mpMesh;
	refMesh.positions = mTetraMesh.positions;
	refMesh.positions.resize(mpMesh->amount);
	if (!MeshUtil::saveMesh(mVisualizeFolder + "simplified-ref.ply", refMesh)) return false;
#endif

	return true;
}

bool DeformVolumetricGraph::smoothGraph() {

	cout << "Smoothing graph" << endl;

	int numSmoothingIterations = 3; // UNDONE: param number of smoothing iterations

	// build neighbor set

	vector<set<int>> vertexNeighbors(mTetraMesh.amount, set<int>());
	for (vec4i idx : mTetraMesh.indices) {
		for (int k = 0; k < 4; k++) {
			for (int j = 1; j <= 3; j++) {
				vertexNeighbors[idx[k]].insert(idx[(k + j) % 4]);
			}
		}
	}

	// mark isolated tetrahedra

	int numFixedPoints = mpMesh->amount;

	if (true) {
		vector<bool> vertexFlags(mTetraMesh.amount, false);
		vector<int> queue(numFixedPoints);
		for (int vertID = 0; vertID < numFixedPoints; vertID++) {
			queue[vertID] = vertID;
			vertexFlags[vertID] = true;
		}

		// mark valid vertices (connected to fixed vertices)
		int head = 0;
		while (head < (int)queue.size()) {
			int vertID = queue[head];
			for (int nbID : vertexNeighbors[vertID]) {
				if (!vertexFlags[nbID]) {
					vertexFlags[nbID] = true;
					queue.push_back(nbID);
				}
			}
			head++;
		}

		// clear neighbors for invalid vertices
		for (int vertID = numFixedPoints; vertID < mTetraMesh.amount; vertID++) {
			if (!vertexFlags[vertID]) vertexNeighbors[vertID].clear();
		}

		// extract valid tetrahedra
		vector<vec4i> validTetra;
		for (vec4i &idx : mTetraMesh.indices) {
			if (vertexFlags[idx[0]]) validTetra.push_back(idx);
		}
		mTetraMesh.indices.swap(validTetra);
	}

	// Laplacian smoothing

	for (int iter = 0; iter < numSmoothingIterations; iter++) {

		vector<vec3> newPosition(mTetraMesh.amount);
		for (int vertID = 0; vertID < mTetraMesh.amount; vertID++) {
			if (vertID < numFixedPoints) {
				newPosition[vertID] = mTetraMesh.positions[vertID];
			} else {
				vec3 pos(0.0f, 0.0f, 0.0f);
				for (int neighborID : vertexNeighbors[vertID]) {
					pos += mTetraMesh.positions[neighborID];
				}
				int numNeighbors = (int)vertexNeighbors[vertID].size();
				if (!vertexNeighbors[vertID].empty()) {
					pos /= (float)vertexNeighbors[vertID].size();
				} else {
					pos = mTetraMesh.positions[vertID];
				}
				newPosition[vertID] = pos;
			}
		}

		mTetraMesh.positions.swap(newPosition);
	}

	if (!MeshUtil::cleanUp(mTetraMesh)) return false;

#ifdef OUTPUT_VISUALIZATION
	TTriangleMesh visualMesh;
	if (!MeshUtil::tet2mesh(mTetraMesh, visualMesh)) return false;
	if (!MeshUtil::saveMesh(mVisualizeFolder + "smoothed.ply", visualMesh)) return false;

	TTriangleMesh refMesh = *mpMesh;
	refMesh.positions = mTetraMesh.positions;
	refMesh.positions.resize(mpMesh->amount);
	if (!MeshUtil::saveMesh(mVisualizeFolder + "smoothed-ref.ply", refMesh)) return false;
#endif

	return true;
}

bool DeformVolumetricGraph::washoutInternal() {

	cout << "Washing out internal tets" << endl;

	vector<vec3> surfacePoints = mTetraMesh.positions;
	surfacePoints.resize(mpMesh->amount);

	SKDTree tree;
	SKDTreeData treeData;
	if (!SampleUtil::buildKdTree(surfacePoints, tree, treeData)) return false;

	double validDistance = StyleSynthesisConfig::mDeform_HandleImpactRadius;

	vector<bool> validFlags(mTetraMesh.amount, false);
#pragma omp parallel for
	for (int vertexID = 0; vertexID < mTetraMesh.amount; vertexID++) {
		if (vertexID < mpMesh->amount) {
			validFlags[vertexID] = true;
			continue;
		}

		vec3 point = mTetraMesh.positions[vertexID];
		SKDT::NamedPoint queryPoint(point[0], point[1], point[2]);
		Thea::BoundedSortedArray<SKDTree::Neighbor> queryResult(1);
		tree.kClosestElements<Thea::MetricL2>(queryPoint, queryResult);
		if (!queryResult.isEmpty()) {
			int neighborID = (int)tree.getElements()[queryResult[0].getIndex()].id;
			float dist = (mTetraMesh.positions[neighborID] - point).length();
			if (dist < validDistance) {
				validFlags[vertexID] = true;
			}
		}
	}
	if (true) {
		vector<bool> newValidFlags = validFlags;
		for (vec4i idx : mTetraMesh.indices) {
			bool valid = false;
			for (int k = 0; k < 4; k++) if (validFlags[idx[k]]) valid = true;
			if (valid) for (int k = 0; k < 4; k++) newValidFlags[idx[k]] = true;
		}
		validFlags.swap(newValidFlags);
	}

	TTetrahedralMesh newMesh;
	int numValidVertices = 0;
	vector<int> validVertexMap(mTetraMesh.amount, -1);
	for (int vertexID = 0; vertexID < mTetraMesh.amount; vertexID++) {
		if (validFlags[vertexID]) {
			newMesh.positions.push_back(mTetraMesh.positions[vertexID]);
			validVertexMap[vertexID] = numValidVertices;
			numValidVertices++;
		}
	}
	newMesh.amount = numValidVertices;
	for (vec4i idx : mTetraMesh.indices) {
		vec4i newIdx;
		bool valid = true;
		for (int k = 0; k < 4; k++) {
			newIdx[k] = validVertexMap[idx[k]];
			if (newIdx[k] < 0) {
				valid = false;
				break;
			}
		}
		if (valid) newMesh.indices.push_back(newIdx);
	}

	mTetraMesh = newMesh;

#ifdef OUTPUT_VISUALIZATION
	TTriangleMesh visualMesh;
	if (!MeshUtil::tet2mesh(mTetraMesh, visualMesh)) return false;
	if (!MeshUtil::saveMesh(mVisualizeFolder + "washed.ply", visualMesh)) return false;

	TTriangleMesh refMesh = *mpMesh;
	refMesh.positions = mTetraMesh.positions;
	refMesh.positions.resize(mpMesh->amount);
	if (!MeshUtil::saveMesh(mVisualizeFolder + "washed-ref.ply", refMesh)) return false;
#endif

	return true;
}

bool DeformVolumetricGraph::checkInside(vec3 &point) {

	const int numRays = 16; // UNDONE: param number of rays for visibility checking
	const double minVisibility = 0.1; // UNDONE: param min visibility (0: any visible)
	int minVisibleRays = (int)(numRays*minVisibility);
	float eps = (float)(mAverageEdgeLength * 0.01);

	int unhitCount = 0;
#pragma omp parallel for shared(unhitCount)
	for (int rayID = 0; rayID<numRays; rayID++) {

		if (unhitCount > minVisibleRays) continue;

		double r1 = cml::random_unit();
		double r2 = cml::random_unit();
		vec3 rayDir = vec3d(
			2.0 * cos(cml::constantsd::two_pi()*r1) * cml::sqrt_safe(r2*(1 - r2)),
			2.0 * sin(cml::constantsd::two_pi()*r1) * cml::sqrt_safe(r2*(1 - r2)),
			1.0 - 2.0*r2); // random direction on unit sphere

		vec3 rayOrigin = point + rayDir * eps; // add eps for offset
		Thea::Ray3 ray(G3D::Vector3(rayOrigin.data()), G3D::Vector3(rayDir.data()));
		float dist = (float)mMeshTree.rayIntersectionTime(ray);

		if (dist<-0.5) { // should be -1 if not hit
#pragma omp atomic
			unhitCount++;
		}
	}

	return unhitCount <= minVisibleRays;
}

int DeformVolumetricGraph::getVertexMap(vector<int> &vmap, int id) {

	vector<int> path(1, id);
	while (true) {
		int currentID = path.back();
		int nextID = vmap[currentID];
		if (nextID != currentID) path.push_back(nextID);
		else break;
	}
	int newID = path.back();
	for (int pathID : path) vmap[pathID] = newID;
	return newID;
}

bool DeformVolumetricGraph::output(TTetrahedralMesh &mesh) {

	mesh = mTetraMesh;

	return true;
}