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

#include "SegmentUtil.h"

#include <fstream>
#include <set>

#include "Sample/SampleUtil.h"
#include "Utility/PlyExporter.h"

using namespace StyleSynthesis;

bool SegmentUtil::saveSegmentationData(string fileName, vector<vector<int>> &segmentData) {

	ofstream outFile(fileName);

	outFile << segmentData.size() << endl;
	for (auto &segment : segmentData) {
		outFile << segment.size() << " ";
		for (auto &sampleID : segment) {
			outFile << sampleID << " ";
		}
		outFile << endl;
	}

	outFile.close();

	return true;
}

bool SegmentUtil::loadSegmentationData(string fileName, vector<vector<int>> &segmentData) {

	ifstream inFile(fileName);
	if (!inFile.is_open()) {
		cout << "Error: cannot load segmentation file " << fileName << endl;
		return false;
	}

	int numSegments = 0;
	inFile >> numSegments;
	segmentData.clear();
	segmentData.reserve(numSegments);

	for (int segID = 0; segID < numSegments; segID++) {

		int numSamples;
		inFile >> numSamples;
		if (!inFile.good()) break;

		vector<int> segment;
		segment.resize(numSamples);
		for (int i = 0; i<numSamples; i++) {
			int sampleID;
			inFile >> sampleID;
			segment[i] = sampleID;
		}
		segmentData.push_back(segment);
	}

	inFile.close();

	return true;
}

bool SegmentUtil::savePatchData(string fileName, vector<vector<int>> &patchData, vector<vector<int>> &patchGraph) {

	if (patchData.size() != patchGraph.size()) {
		cout << "Error: inconsistent number of patches" << endl;
		return false;
	}

	ofstream outFile(fileName);

	outFile << patchData.size() << endl;

	for (auto &patch : patchData) {
		outFile << patch.size() << " ";
		for (auto &sampleID : patch) {
			outFile << sampleID << " ";
		}
		outFile << endl;
	}

	for (auto &node : patchGraph) {
		outFile << node.size() << " ";
		for (auto &neighborID : node) {
			outFile << neighborID << " ";
		}
		outFile << endl;
	}

	outFile.close();

	return true;
}

bool SegmentUtil::loadPatchData(string fileName, vector<vector<int>> &patchData, vector<vector<int>> &patchGraph) {

	ifstream inFile(fileName);
	if (!inFile.is_open()) {
		cout << "Error: cannot load patch file " << fileName << endl;
		return false;
	}

	int numPatches = 0;
	inFile >> numPatches;

	patchData.clear();
	patchData.resize(numPatches);
	for (auto &patch : patchData) {

		int numSamples;
		inFile >> numSamples;
		if (!inFile.good()) break;

		patch.clear();
		patch.resize(numSamples);
		for (int i = 0; i<numSamples; i++) {
			int sampleID;
			inFile >> sampleID;
			patch[i] = sampleID;
		}
	}

	patchGraph.clear();
	patchGraph.resize(numPatches);
	for (auto &node : patchGraph) {

		int numNeighbors;
		inFile >> numNeighbors;
		if (!inFile.good()) break;

		node.clear();
		node.resize(numNeighbors);
		for (int i = 0; i<numNeighbors; i++) {
			int neighborID;
			inFile >> neighborID;
			node[i] = neighborID;
		}
	}

	inFile.close();

	return true;
}

bool SegmentUtil::saveLabelData(string fileName, vector<int> &labelData) {

	ofstream outFile(fileName);

	outFile << labelData.size() << endl;
	for (int label : labelData) {
		outFile << label << endl;
	}

	outFile.close();

	return true;
}

bool SegmentUtil::loadLabelData(string fileName, vector<int> &labelData) {

	ifstream inFile(fileName);
	if (!inFile.is_open()) {
		cout << "Error: cannot load label file " << fileName << endl;
		return false;
	}

	int numLabels = 0;
	inFile >> numLabels;
	labelData.resize(numLabels);

	for (int labelID = 0; labelID < numLabels; labelID++) {		
		int label;
		inFile >> label;
		labelData[labelID] = label;
	}

	inFile.close();

	return true;
}

bool SegmentUtil::visualizeSegmentHierarchy(string fileName, TSampleSet &sample, vector<vector<int>> &patch, vector<vector<int>> &segment) {

	vec3 bbMin, bbMax;
	if (!SampleUtil::computeAABB(sample, bbMin, bbMax)) return false;
	vec3 space = (bbMax - bbMin)*1.2f;

	int numPatches = (int)patch.size();
	int numSegments = (int)segment.size();
	vector<bool> patchFlag(numPatches, false);

	PlyExporter pe;
	int shapeID = 0;
	for (int segmentID = 0; segmentID < numSegments; segmentID++) {

		vector<vec3> vp, vn;
		for (int patchID : segment[segmentID]) {
			if (patchFlag[patchID]) {
				patchFlag.assign(numPatches, false);
				patchFlag[patchID] = true;
				shapeID++;
			} else {
				patchFlag[patchID] = true;
			}
			for (int sampleID : patch[patchID]) {
				vp.push_back(sample.positions[sampleID]);
				vn.push_back(sample.normals[sampleID]);
			}
		}
		if (shapeID >= 27) break;
		vec3 offset(space[0] * ((shapeID % 9) % 3), space[1] * (shapeID / 9), space[2] * ((shapeID % 9) / 3));
		vec3i color = SegmentUtil::colorMapping(segmentID);
		if (!pe.addPoint(&vp, &vn, offset, color)) return false;
	}

	if (!pe.output(fileName)) return false;

	return true;
}

bool SegmentUtil::visualizeSegmentOnMesh(string fileName, TTriangleMesh &mesh, vector<vector<int>> &segment) {

	vector<vec3i> faceColors(mesh.indices.size(), vec3i(0, 0, 0));
	int numSegments = (int)segment.size();
	for (int segmentID = 0; segmentID < numSegments; segmentID++) {
		vec3i color = SegmentUtil::colorMapping(segmentID);
		for (int faceID : segment[segmentID]) {
			faceColors[faceID] = color;
		}
	}

	PlyExporter pe;
	if (!pe.addMesh(&mesh.positions, &mesh.normals, &mesh.indices, &faceColors)) return false;
	if (!pe.output(fileName)) return false;

	return true;
}

bool SegmentUtil::buildKNNGraph(
	TSampleSet &inSamples,
	vector<vector<int>> &outGraph,
	vector<bool> &outFlag,
	double optRadius,
	bool optNormalCheck)
{
	SKDTree tree;
	SKDTreeData treeData;
	if (!SampleUtil::buildKdTree(inSamples.positions, tree, treeData)) return false;

	vector<set<int>> sampleCloseNeighbors(inSamples.amount);

	// build graph
	const int numNeighbors = 7; // used for build graph
#pragma omp parallel for
	for (int sampleID = 0; sampleID < inSamples.amount; sampleID++) {
		vec3 sampleP = inSamples.positions[sampleID];
		vec3 sampleN = inSamples.normals[sampleID];
		SKDT::NamedPoint queryPoint(sampleP[0], sampleP[1], sampleP[2]);
		int capacity = (int)(3 * optRadius*(optRadius + 1) + 1);
		Thea::BoundedSortedArray<SKDTree::Neighbor> queryResult(capacity);
		tree.kClosestElements<Thea::MetricL2>(queryPoint, queryResult, inSamples.radius*optRadius);
		if (queryResult.size() < numNeighbors) {
			queryResult.clear();
			queryResult = Thea::BoundedSortedArray<SKDTree::Neighbor>(numNeighbors);
			tree.kClosestElements<Thea::MetricL2>(queryPoint, queryResult);
		}
		for (int id = 0; id < queryResult.size(); id++) {
			int neighborID = (int)tree.getElements()[queryResult[id].getIndex()].id;
			vec3 neighborN = inSamples.normals[neighborID];
			if (sampleID == neighborID) continue;
			if (optNormalCheck && cml::dot(sampleN, neighborN) < -0.01f) continue;
			sampleCloseNeighbors[sampleID].insert(neighborID);
		}
	}

	// detect outliers
	outFlag.assign(inSamples.amount, true);
#pragma omp parallel for
	for (int sampleID = 0; sampleID < inSamples.amount; sampleID++) {
		vec3 samplePos = inSamples.positions[sampleID];
		SKDT::NamedPoint queryPoint(samplePos[0], samplePos[1], samplePos[2]);
		Thea::BoundedSortedArray<SKDTree::Neighbor> queryResult(2);
		tree.kClosestElements<Thea::MetricL2>(queryPoint, queryResult);
		if (queryResult.size() > 1) {
			int neighborID = (int)tree.getElements()[queryResult[1].getIndex()].id;
			auto &topSet = sampleCloseNeighbors[neighborID];
			if (topSet.find(sampleID) == topSet.end()) outFlag[sampleID] = false;
		}
		else {
			outFlag[sampleID] = false;
		}
	}

	// make edges undirected
	for (int sampleID = 0; sampleID < inSamples.amount; sampleID++) {
		auto &nbSet = sampleCloseNeighbors[sampleID];
		if (outFlag[sampleID]) {			
			vector<int> nbList(nbSet.begin(), nbSet.end());
			nbSet.clear();
			for (int nbID : nbList) {
				if (outFlag[nbID]) {
					nbSet.insert(nbID);
					sampleCloseNeighbors[nbID].insert(sampleID);					
				}
			}
		} else {
			nbSet.clear();
		}
		if (nbSet.empty()) outFlag[sampleID] = false;
	}

	// prune small components

	vector<bool> visitFlag(inSamples.amount, false);
	for (int sampleID = 0; sampleID < inSamples.amount; sampleID++) {
		if (visitFlag[sampleID]) continue;
		visitFlag[sampleID] = true;

		// find connected component by BFS
		vector<int> queue(1, sampleID);
		int head = 0;
		while (head < (int)queue.size()) {
			int currentID = queue[head];
			for (int neighborID : sampleCloseNeighbors[currentID]) {
				if (!visitFlag[neighborID]) {
					queue.push_back(neighborID);
					visitFlag[neighborID] = true;
				}
			}
			head++;
		}

		if ((int)queue.size() < 50) { // UNDONE: param minimum size for connected component
			for (int currentID : queue) {
				outFlag[currentID] = false;
				sampleCloseNeighbors[currentID].clear();
			}
		}
	}

	// export graph
	outGraph.resize(inSamples.amount);
#pragma omp parallel for
	for (int sampleID = 0; sampleID < inSamples.amount; sampleID++) {
		if (outFlag[sampleID]) {
			auto &nbSet = sampleCloseNeighbors[sampleID];
			outGraph[sampleID].assign(nbSet.begin(), nbSet.end());
		} else {
			outGraph[sampleID].clear();
		}
	}

	return true;
}

vec3i SegmentUtil::colorMapping(int index) {

	//return vec3i(index%256, (index/256)%256, (index/256)/256);

	static vector<vec3i> colorMap(0);
	if (colorMap.size() == 0) {
		// random permutation of HSV color map
		// 0: red
		// 1: cyan
		// 2: green
		// 3: purple
		// 4: yellow
		// 5: blue
		// 6: dark green
		// 7: pink
		// 8: light cyan
		// 9: dark blue
		colorMap.push_back(vec3i(255, 0, 0));      // 0
		colorMap.push_back(vec3i(0, 128, 128));    // 1
		colorMap.push_back(vec3i(0, 255, 0));      // 2
		colorMap.push_back(vec3i(128, 0, 128));    // 3
		colorMap.push_back(vec3i(255, 255, 0));    // 4
		colorMap.push_back(vec3i(0, 0, 255));      // 5
		colorMap.push_back(vec3i(0, 128, 0));      // 6
		colorMap.push_back(vec3i(255, 0, 255));    // 7
		colorMap.push_back(vec3i(0, 255, 255));    // 8
		colorMap.push_back(vec3i(0, 0, 128));      // 9
		colorMap.push_back(vec3i(128, 0, 0));      // 10
		colorMap.push_back(vec3i(176, 255, 0));    // 11
		colorMap.push_back(vec3i(128, 128, 0));    // 12
		colorMap.push_back(vec3i(111, 0, 255));    // 13
		colorMap.push_back(vec3i(0, 120, 255));    // 14
		colorMap.push_back(vec3i(0, 255, 184));    // 15
		colorMap.push_back(vec3i(255, 0, 141));    // 16
		colorMap.push_back(vec3i(255, 103, 0));    // 17
		colorMap.push_back(vec3i(0, 107, 255));    // 18
		colorMap.push_back(vec3i(255, 231, 0));    // 19
		colorMap.push_back(vec3i(60, 0, 255));     // 20
		colorMap.push_back(vec3i(21, 255, 0));
		colorMap.push_back(vec3i(21, 0, 255));
		colorMap.push_back(vec3i(124, 255, 0));
		colorMap.push_back(vec3i(73, 255, 0));
		colorMap.push_back(vec3i(0, 255, 107));
		colorMap.push_back(vec3i(255, 0, 77));
		colorMap.push_back(vec3i(137, 255, 0));
		colorMap.push_back(vec3i(0, 255, 43));
		colorMap.push_back(vec3i(255, 0, 26));
		colorMap.push_back(vec3i(0, 255, 159));
		colorMap.push_back(vec3i(0, 236, 255));
		colorMap.push_back(vec3i(0, 210, 255));
		colorMap.push_back(vec3i(34, 0, 255));
		colorMap.push_back(vec3i(201, 255, 0));
		colorMap.push_back(vec3i(9, 255, 0));
		colorMap.push_back(vec3i(86, 0, 255));
		colorMap.push_back(vec3i(255, 13, 0));
		colorMap.push_back(vec3i(201, 0, 255));
		colorMap.push_back(vec3i(240, 0, 255));
		colorMap.push_back(vec3i(255, 0, 51));
		colorMap.push_back(vec3i(0, 255, 197));
		colorMap.push_back(vec3i(0, 255, 56));
		colorMap.push_back(vec3i(189, 255, 0));
		colorMap.push_back(vec3i(176, 0, 255));
		colorMap.push_back(vec3i(255, 90, 0));
		colorMap.push_back(vec3i(0, 133, 255));
		colorMap.push_back(vec3i(255, 0, 90));
		colorMap.push_back(vec3i(34, 255, 0));
		colorMap.push_back(vec3i(255, 0, 154));
		colorMap.push_back(vec3i(214, 255, 0));
		colorMap.push_back(vec3i(189, 0, 255));
		colorMap.push_back(vec3i(0, 30, 255));
		colorMap.push_back(vec3i(99, 0, 255));
		colorMap.push_back(vec3i(255, 64, 0));
		colorMap.push_back(vec3i(255, 0, 206));
		colorMap.push_back(vec3i(255, 0, 39));
		colorMap.push_back(vec3i(240, 255, 0));
		colorMap.push_back(vec3i(255, 0, 0));
		colorMap.push_back(vec3i(0, 4, 255));
		colorMap.push_back(vec3i(60, 255, 0));
		colorMap.push_back(vec3i(0, 69, 255));
		colorMap.push_back(vec3i(255, 116, 0));
		colorMap.push_back(vec3i(255, 0, 13));
		colorMap.push_back(vec3i(86, 255, 0));
		colorMap.push_back(vec3i(253, 0, 255));
		colorMap.push_back(vec3i(0, 223, 255));
		colorMap.push_back(vec3i(255, 180, 0));
		colorMap.push_back(vec3i(255, 77, 0));
		colorMap.push_back(vec3i(255, 0, 219));
		colorMap.push_back(vec3i(99, 255, 0));
		colorMap.push_back(vec3i(0, 255, 171));
		colorMap.push_back(vec3i(0, 255, 236));
		colorMap.push_back(vec3i(0, 255, 223));
		colorMap.push_back(vec3i(255, 39, 0));
		colorMap.push_back(vec3i(0, 255, 133));
		colorMap.push_back(vec3i(227, 0, 255));
		colorMap.push_back(vec3i(0, 255, 69));
		colorMap.push_back(vec3i(0, 159, 255));
		colorMap.push_back(vec3i(255, 219, 0));
		colorMap.push_back(vec3i(255, 0, 231));
		colorMap.push_back(vec3i(73, 0, 255));
		colorMap.push_back(vec3i(0, 184, 255));
		colorMap.push_back(vec3i(0, 255, 81));
		colorMap.push_back(vec3i(0, 146, 255));
		colorMap.push_back(vec3i(255, 0, 129));
		colorMap.push_back(vec3i(255, 0, 193));
		colorMap.push_back(vec3i(0, 94, 255));
		colorMap.push_back(vec3i(0, 56, 255));
		colorMap.push_back(vec3i(0, 171, 255));
		colorMap.push_back(vec3i(111, 255, 0));
		colorMap.push_back(vec3i(0, 255, 210));
		colorMap.push_back(vec3i(255, 167, 0));
		colorMap.push_back(vec3i(0, 255, 30));
		colorMap.push_back(vec3i(0, 255, 249));
		colorMap.push_back(vec3i(137, 0, 255));
		colorMap.push_back(vec3i(255, 206, 0));
		colorMap.push_back(vec3i(255, 193, 0));
		colorMap.push_back(vec3i(163, 255, 0));
		colorMap.push_back(vec3i(0, 255, 146));
		colorMap.push_back(vec3i(163, 0, 255));
		colorMap.push_back(vec3i(255, 26, 0));
		colorMap.push_back(vec3i(0, 255, 17));
		colorMap.push_back(vec3i(0, 197, 255));
		colorMap.push_back(vec3i(150, 255, 0));
		colorMap.push_back(vec3i(0, 43, 255));
		colorMap.push_back(vec3i(255, 129, 0));
		colorMap.push_back(vec3i(255, 0, 64));
		colorMap.push_back(vec3i(255, 141, 0));
		colorMap.push_back(vec3i(47, 0, 255));
		colorMap.push_back(vec3i(255, 154, 0));
		for (int k = 0; k < 1000; k++) {
			vec3i randColor(cml::random_integer(0, 255), cml::random_integer(0, 255), cml::random_integer(0, 255));
			colorMap.push_back(randColor);
		}
	}

	if (index < 0) return vec3i(0, 0, 0);
	if (index < (int)colorMap.size()) return colorMap[index];

	return vec3i(cml::random_integer(0, 255), cml::random_integer(0, 255), cml::random_integer(0, 255));
}