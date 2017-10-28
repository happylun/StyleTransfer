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

#include "SegmentGroupApxCvx.h"

#include <fstream>

#include "Mesh/MeshUtil.h"

#include "Segment/SegmentMeshApxCvx.h"
#include "Segment/SegmentMeshPriFit.h"
#include "Segment/SegmentUtil.h"

#include "Utility/PlyExporter.h"

#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

#define OUTPUT_PROGRESS

bool SegmentGroupApxCvx::runSegmentation(
	vector<TTriangleMesh> &inGroups,
	TTriangleMesh &outMesh,
	vector<vector<vector<int>>> &outSegments)
{

	int numGroups = (int)inGroups.size();

	auto &segmentVisibilities = StyleSynthesisConfig::mSegmentation_VisibilityThresholds;
	int numSegmentLevels = (int)segmentVisibilities.values.size() + 1;
	//if (true) {
	//	cout << numSegmentLevels << ": ";
	//	for (auto v : segmentVisibilities.values) cout << v << " ";
	//	cout << endl;
	//	system("pause");
	//}

	outSegments.resize(numSegmentLevels, vector<vector<int>>(0));

	int numAllFaces = 0;
	PlyExporter meshCollector;
	vector<TTriangleMesh> groupMesh(numGroups);

	for (int groupID = 0; groupID < numGroups; groupID++) {
#ifdef OUTPUT_PROGRESS
		cout << "\rSegmenting group " << (groupID + 1) << " / " << numGroups << "           ";
#endif

		SegmentMeshApxCvx smac;
		if (!smac.initSegmentation(inGroups[groupID])) return false;

		TTriangleMesh &segmentMesh = groupMesh[groupID];
		if (!smac.exportMesh(segmentMesh)) return false;
		if (!meshCollector.addMesh(&segmentMesh.indices, &segmentMesh.positions, &segmentMesh.normals)) return false;

		for (int segLevel = 1; segLevel < numSegmentLevels; segLevel++) {
			int levelID = numSegmentLevels - segLevel; // reverse level: fine to coarse -> coarse to fine
			vector<vector<int>> segments;
			if (!smac.runSegmentation(segmentVisibilities.values[segLevel - 1])) return false;
			if (!smac.exportSegmentation(segments)) return false;

			for (auto &segment : segments) {
				for (auto &faceID : segment) {
					faceID += numAllFaces;
				}
			}
			auto &levelSegments = outSegments[levelID];
			levelSegments.insert(levelSegments.end(), segments.begin(), segments.end());
		}

		int numFaces = (int)segmentMesh.indices.size();
		vector<int> segments(numFaces);
		for (int id = 0; id < numFaces; id++) {
			segments[id] = numAllFaces + id;
		}
		outSegments[0].push_back(segments);
		numAllFaces += numFaces;
	}
#ifdef OUTPUT_PROGRESS
	cout << endl;
#endif

	// finalize mesh

	outMesh.positions = meshCollector.mVertices;
	outMesh.normals = meshCollector.mNormals;
	outMesh.indices = meshCollector.mFaceIndices;
	outMesh.amount = (int)outMesh.positions.size();
	//if (!MeshUtil::reorientAnyMesh(outMesh)) return false;

	// process first level segments

	int minNumGroups = 2; // UNDONE: param good grouping criteria
	int maxNumGroups = 15;

#ifdef OUTPUT_PROGRESS
	cout << "Processing first level segments" << endl;
#endif
	if (numGroups >= minNumGroups && numGroups <= maxNumGroups) {
		// keep original grouping as first level segmentation
#ifdef OUTPUT_PROGRESS
		cout << "Using original grouping..." << endl;
#endif
	} else {
		// original grouping is bad -- use primitive fitting segmentation
#ifdef OUTPUT_PROGRESS
		cout << "Running primitive fitting segmentation..." << endl;
#endif
		SegmentMeshPriFit smpf;
		int levelID;
		if (numGroups < minNumGroups) {
			levelID = 1; // too few groups -- merge lower level patches
		} else {
			levelID = 0; // too many groups -- merge groups
		}
		if (!smpf.loadData(outMesh, outSegments[levelID])) return false;
		if (!smpf.runSegmentation()) return false;
		if (!smpf.exportSegmentation(outSegments[0])) return false;
	}

	return true;
}

bool SegmentGroupApxCvx::extractUniqueSegments(
	TTriangleMesh &inMesh,
	vector<vector<vector<int>>> &inSegments,
	vector<vector<int>> &outSegments)
{

#ifdef OUTPUT_PROGRESS
	cout << "Extracting unique segments..." << endl;
#endif

	int numLevels = (int)inSegments.size();
	int numGroups = (int)inSegments[0].size();
	int numFaces = (int)inMesh.indices.size();

	outSegments.clear();

	// process first level -- group level

	vector<int> lastLevelFaceMap(numFaces, -1);	
	for (int groupID = 0; groupID < numGroups; groupID++) {
		auto &group = inSegments[0][groupID];
		for (int faceID : group) lastLevelFaceMap[faceID] = groupID;
		outSegments.push_back(group);
	}

	// process remaining level

	vector<int> currentLevelFaceMap(numFaces, -1);
	for (int levelID = 1; levelID < numLevels; levelID++) {
		int numSegments = (int)inSegments[levelID].size();
		for (int segID = 0; segID < numSegments; segID++) {
			auto &segment = inSegments[levelID][segID];
			int segmentSize = (int)segment.size();

			int parentSegID = lastLevelFaceMap[segment[0]];
			int parentSegmentSize = (int)inSegments[levelID - 1][parentSegID].size();
			if (segmentSize != parentSegmentSize) outSegments.push_back(segment);

			for (int faceID : segment) currentLevelFaceMap[faceID] = segID;
		}

		lastLevelFaceMap.swap(currentLevelFaceMap);
		currentLevelFaceMap.assign(numFaces, -1);
	}

	// HACK: quickly prune flat segments

	if (true) {

		TTriangleMesh segmentMesh;
		segmentMesh.positions = inMesh.positions;
		segmentMesh.normals = inMesh.normals;
		segmentMesh.amount = (int)segmentMesh.positions.size();

		vector<vector<int>> prunedSegments(0);

		int numOutSegments = (int)outSegments.size();
		for (int segID = 0; segID < numOutSegments; segID++) {
			auto &segment = outSegments[segID];
			int numFaces = (int)segment.size();
			if (numFaces == 0) continue;
			segmentMesh.indices.resize(numFaces);
			for (int id = 0; id < numFaces; id++) {
				segmentMesh.indices[id] = inMesh.indices[segment[id]];
			}
			
			double volume;
			if (!MeshUtil::computeVolume(segmentMesh, volume)) return false;
			if (volume > 1e-6f) {
				prunedSegments.push_back(segment);
			}
		}

		outSegments.swap(prunedSegments);
	}

	return true;
}

bool SegmentGroupApxCvx::findElementMapping(
	vector<vector<vector<int>>> &inSegments,
	vector<vector<int>> &inElements,
	vector<vector<int>> &outMapping)
{
	// this can be done when extracting elements from segments but...

	int numElements = (int)inElements.size();
	vector<set<int>> elementFaceSets(numElements, set<int>());
	for (int elementID = 0; elementID < numElements; elementID++) {
		elementFaceSets[elementID].insert(inElements[elementID].begin(), inElements[elementID].end());
	}

	int numLevels = (int)inSegments.size();
	outMapping.resize(numLevels);
	for (int levelID = 0; levelID < numLevels; levelID++) {
		int numSegments = (int)inSegments[levelID].size();
		outMapping[levelID].assign(numSegments, -1);
		for (int segID = 0; segID < numSegments; segID++) {
			vector<int> &segment = inSegments[levelID][segID];
			for (int elementID = 0; elementID < numElements; elementID++) {
				auto &faceSet = elementFaceSets[elementID];
				if (segment.size() != faceSet.size()) continue;
				if (faceSet.find(segment[0]) == faceSet.end()) continue;
				outMapping[levelID][segID] = elementID;
				break;
			}
		}
	}

	return true;
}

bool SegmentGroupApxCvx::visualize(
	string fileName,
	TTriangleMesh &mesh,
	vector<vector<vector<int>>> &segments)
{

	vec3 bbMin, bbMax;
	if (!MeshUtil::computeAABB(mesh, bbMin, bbMax)) return false;

	PlyExporter pe;
	for (int levelID = 0; levelID < (int)segments.size(); levelID++) {
		auto &levelSegments = segments[levelID];
				
		vector<vec3i> faceColors(mesh.indices.size(), vec3i(0, 0, 0));		
		int numSegments = (int)levelSegments.size();
		for (int segmentID = 0; segmentID < numSegments; segmentID++) {
			vec3i color = SegmentUtil::colorMapping(segmentID);
			for (int faceID : levelSegments[segmentID]) {
				faceColors[faceID] = color;
			}
		}

		vec3 offset((bbMax[0] - bbMin[0])*levelID*1.5f, 0.0f, 0.0f);
		if (!pe.addMesh(&mesh.positions, &mesh.normals, &mesh.indices, &faceColors, offset)) return false;
	}
	if (!pe.output(fileName)) return false;

	return true;
}

bool SegmentGroupApxCvx::saveSegments(string fileName, vector<vector<vector<int>>> &segments) {

	ofstream file(fileName);

	int numLevels = (int)segments.size();
	file << numLevels << endl;

	for (auto &levelSegments : segments) {
		int numSegments = (int)levelSegments.size();
		file << numSegments << endl;
		for (auto &segment : levelSegments) {
			int numFaces = (int)segment.size();
			file << numFaces;
			for (int faceID : segment) file << " " << faceID;
			file << endl;
		}
	}

	file.close();

	return true;
}

bool SegmentGroupApxCvx::loadSegments(string fileName, vector<vector<vector<int>>> &segments) {

	ifstream file(fileName);

	int numLevels;
	file >> numLevels;
	segments.resize(numLevels);

	for (auto &levelSegments : segments) {
		int numSegments;
		file >> numSegments;
		levelSegments.clear();
		for (int segmentID = 0; segmentID < numSegments; segmentID++) {
			vector<int> segment;
			int numFaces;
			file >> numFaces;
			if (numFaces == 0) continue;
			segment.resize(numFaces);
			for (int id = 0; id < numFaces; id++) file >> segment[id];
			levelSegments.push_back(segment);
		}
	}

	file.close();

	return true;
}