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

#include "PipelineSegmentIO.h"

#include <iostream>
#include <fstream>

#include "Segment/SegmentGroupApxCvx.h"
#include "Segment/SegmentUtil.h"

#include "Mesh/MeshUtil.h"

#include "Data/StyleSynthesisConfig.h"

using namespace std;
using namespace StyleSynthesis;

bool PipelineSegmentIO::process() {

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string meshListFileName = datasetPrefix + "mesh/mesh-list-all.txt";
	vector<string> meshNameList;
	ifstream meshListFile(meshListFileName);
	while (!meshListFile.eof()) {
		string line;
		getline(meshListFile, line);
		if (line.empty()) break;
		meshNameList.push_back(StringUtil::trim(line));
	}
	meshListFile.close();

	int numMesh = (int)meshNameList.size();

	string meshLabelFileName = datasetPrefix + "mesh/mesh-labels-all.txt";
	vector<int> meshLabelList(0);
	ifstream meshLabelFile(meshLabelFileName);
	for (int meshID = 0; meshID < numMesh; meshID++) {
		int label;
		meshLabelFile >> label;
		meshLabelList.push_back(label);
	}
	meshLabelFile.close();

	for (int meshID = 0; meshID < numMesh; meshID++) {

		string meshName = meshNameList[meshID];
		cout << "Processing segmention " << meshName << endl;
		if (!runSingleModelSegmentation(meshName)) return error("single model segmentation error");

		//system("pause");
	}

	return true;
}

bool PipelineSegmentIO::runSingleModelSegmentation(string meshName) {

	// names...

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string partFolder = datasetPrefix + "scaled-part/" + meshName + "/";

	string segmentMeshName = datasetPrefix + "segment/" + meshName + ".ply";
	string segmentDataName = datasetPrefix + "segment/" + meshName + "-segment.txt";
	string segmentElementName = datasetPrefix + "segment/" + meshName + "-element.txt";
	string segmentVisName = datasetPrefix + "segment/" + meshName + "-segment.ply";
	string segmentVisPrefix = datasetPrefix + "segment/" + meshName + "-segment";
	if (!FileUtil::makedir(segmentDataName)) return false;

	if (FileUtil::existsfile(segmentElementName)) return true; // early quit

	// data

	vector<TTriangleMesh> groups;
	TTriangleMesh mesh;
	vector<vector<vector<int>>> segments;
	vector<vector<int>> elements;

	// load part groups

	int numGroups = 0;
	while (true) {
		string partName = partFolder + to_string(numGroups) + ".ply";
		if (!FileUtil::existsfile(partName)) break;
		groups.push_back(TTriangleMesh());
		if (!MeshUtil::loadMesh(partName, groups.back())) return false;
		numGroups++;
	}

	// algorithm

	if (FileUtil::existsfile(segmentDataName)) {
		if (!MeshUtil::loadMesh(segmentMeshName, mesh)) return false;
		if (!SegmentGroupApxCvx::loadSegments(segmentDataName, segments)) return false;
	} else {
		if (!SegmentGroupApxCvx::runSegmentation(groups, mesh, segments)) return false;
		if (!SegmentGroupApxCvx::visualize(segmentVisName, mesh, segments)) return false;
		if (!SegmentGroupApxCvx::saveSegments(segmentDataName, segments)) return false;
		if (!MeshUtil::saveMesh(segmentMeshName, mesh)) return false;
	}

	if (!FileUtil::existsfile(segmentElementName)) {
		// HACK: only use first 3 levels
		vector<vector<vector<int>>> segmentsOfInterest(segments.begin(), segments.begin() + 3);
		if (!SegmentGroupApxCvx::extractUniqueSegments(mesh, segmentsOfInterest, elements)) return false;
		if (!SegmentUtil::saveSegmentationData(segmentElementName, elements)) return false;
	}

	string anyVisualName = segmentVisPrefix + "-1.ply";
	if (!FileUtil::existsfile(anyVisualName)) {
		int visualLevel = min(3, (int)segments.size());
		if ((int)segments.size() >= visualLevel) {
			for (int level = 0; level < visualLevel; level++) {
				string visualName = segmentVisPrefix + "-" + to_string(level + 1) + ".ply";
				if (!SegmentUtil::visualizeSegmentOnMesh(visualName, mesh, segments[level])) return false;
			}
		}
	}

	return true;
}

bool PipelineSegmentIO::error(string s) {
	cout << "Error: " << s << endl;
	return false;
}