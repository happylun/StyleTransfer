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

#include "PipelineGraphIO.h"

#include <iostream>
#include <fstream>

#include "Segment/SegmentGroupApxCvx.h"
#include "Segment/SegmentUtil.h"

#include "Context/ContextPartGraph.h"
#include "Context/ContextPartGraphNodeGenerator.h"
#include "Context/ContextPartGraphMatch.h"

#include "Mesh/MeshUtil.h"

#include "Data/DataUtil.h"
#include "Data/StyleSynthesisConfig.h"

using namespace std;
using namespace StyleSynthesis;

bool PipelineGraphIO::process() {

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
		int meshLabel = meshLabelList[meshID];
		cout << "Processing context graph " << meshName << endl;
		if (!runSingleModelContextGraph(meshName, meshLabel)) return error("single model context graph error");

		//system("pause");
	}

	if (!computeSelfSimilarity(meshNameList)) return false;

	return true;
}

bool PipelineGraphIO::runSingleModelContextGraph(string meshName, int meshLabel) {

	// names...

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string partFolder = datasetPrefix + "scaled-part/" + meshName + "/";

	string segmentMeshName = datasetPrefix + "segment/" + meshName + ".ply";
	string segmentDataName = datasetPrefix + "segment/" + meshName + "-segment.txt";
	//string segmentMeshCopyName = datasetPrefix + "segment/" + meshName + "-copy.ply";

	string graphHierarchyName = datasetPrefix + "graph/" + meshName + "/graph-hierarchy.txt";
	string graphDescriptorName = datasetPrefix + "graph/" + meshName + "/graph-descriptor.txt";
	string graphContextName = datasetPrefix + "graph/" + meshName + "/graph-context.txt";
	if (!FileUtil::makedir(graphContextName)) return false;

	if (FileUtil::existsfile(graphContextName)) return true; // early quit

	// data

	TTriangleMesh mesh;
	vector<vector<vector<int>>> segments;
	ContextPartGraph graph;
	ContextPartGraphNodeGenerator graphNodeGenerator;

	// algorithm

	if (!FileUtil::existsfile(segmentDataName)) {
		return error("run single model segmentation first");
	} else {
		if (!MeshUtil::loadMesh(segmentMeshName, mesh)) return false;
		if (!SegmentGroupApxCvx::loadSegments(segmentDataName, segments)) return false;
	}

	int graphMode = 0;
	if (true) {
		// special handling for graph construction
		if (meshName[0] == 'L' && meshLabel == 4) {
			//// ceiling lamp
			//if (!FileUtil::existsfile(segmentMeshCopyName)) {
			//	if (!MeshUtil::saveMesh(segmentMeshCopyName, mesh)) return false;
			//	// flip along Y axis
			//	for (vec3 &v : mesh.positions) v[1] = -v[1];
			//	for (vec3 &n : mesh.normals) n[1] = -n[1];
			//	for (vec3i &i : mesh.indices) swap(i[1], i[2]);
			//	if (!MeshUtil::saveMesh(segmentMeshName, mesh)) return false;
			//}
			graphMode = 1;
		}
		if (meshName[0] == 'L' && meshLabel == 2) {
			//// wall lamp
			//if (!FileUtil::existsfile(segmentMeshCopyName)) {
			//	if (!MeshUtil::saveMesh(segmentMeshCopyName, mesh)) return false;
			//	// rotate to make the "wall" look like the "ground"
			//	for (vec3 &v : mesh.positions) {
			//		swap(v[1], v[2]);
			//		v[2] = -v[2];
			//	}
			//	for (vec3 &n : mesh.normals) {
			//		swap(n[1], n[2]);
			//		n[2] = -n[2];
			//	}
			//	if (!MeshUtil::saveMesh(segmentMeshName, mesh)) return false;
			//}
			//graphMode = 2;
		}
	}

	if (FileUtil::existsfile(graphHierarchyName)) {
		if (!graph.loadGraphHierarchy(graphHierarchyName, graphNodeGenerator, &mesh, &segments)) return false;
	} else {
		if (!graph.buildGraphHierarchy(graphNodeGenerator, &mesh, &segments)) return false;
		if (!graph.saveGraphHierarchy(graphHierarchyName)) return false;
	}
	if (FileUtil::existsfile(graphDescriptorName)) {
		if (!graph.loadGraphDescriptor(graphDescriptorName)) return false;
	} else {
		if (!graph.buildGraphDescriptor(graphMode)) return false;
		if (!graph.saveGraphDescriptor(graphDescriptorName)) return false;
	}
	if (FileUtil::existsfile(graphContextName)) {
		if (!graph.loadGraphContext(graphContextName)) return false;
	} else {
		if (!graph.buildGraphContext(graphMode)) return false;
		if (!graph.saveGraphContext(graphContextName)) return false;
		if (!graph.saveGraphDescriptor(graphDescriptorName)) return false; // NOTE: may be updated depending on context
	}

	return true;
}

bool PipelineGraphIO::computeSelfSimilarity(vector<string> &meshNameList) {

	// compute self similarity

	int numMeshes = (int)meshNameList.size();

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	for (int meshID = 0; meshID < numMeshes; meshID++) {

		string meshName = meshNameList[meshID];

		cout << "================ Computing self similarity " << meshName << " ================" << endl;

		// names

		string meshFileName = datasetPrefix + "segment/" + meshName + ".ply";
		string segmentName = datasetPrefix + "segment/" + meshName + "-segment.txt";

		string graphFolder = datasetPrefix + "graph/" + meshName + "/";
		string graphHName = graphFolder + "graph-hierarchy.txt";
		string graphDName = graphFolder + "graph-descriptor.txt";
		string graphCName = graphFolder + "graph-context.txt";

		string weightsFolder = datasetPrefix + "weights/";

		string outSimilarityName = graphFolder + "graph-self-similarity.txt";
		string outModeName = graphFolder + "graph-self-mode.txt";
		if (FileUtil::existsfile(outModeName)) continue; // early quit

		// data

		TTriangleMesh mesh;
		vector<vector<vector<int>>> segment;
		ContextPartGraphNodeGenerator cpgng;
		ContextPartGraph graph;

		if (!MeshUtil::loadMesh(meshFileName, mesh)) return false;
		if (!SegmentGroupApxCvx::loadSegments(segmentName, segment)) return false;

		if (!graph.loadGraphHierarchy(graphHName, cpgng, &mesh, &segment)) return false;
		if (!graph.loadGraphDescriptor(graphDName)) return false;
		if (!graph.loadGraphContext(graphCName)) return false;

		// matching

		Eigen::MatrixXd similarityMatrix;
		int matchingMode;
		ContextPartGraphMatch cpgm;
		if (!cpgm.loadWeights(weightsFolder)) return false;
		if (!cpgm.loadGraph(graph, graph)) return false;
		if (!cpgm.process()) return false;
		if (!cpgm.exportSimilarityMatrix(similarityMatrix)) return false;
		if (!cpgm.exportMatchingMode(matchingMode)) return false;

		if (!DataUtil::saveMatrixBinary(outSimilarityName, similarityMatrix)) return false;
		if (!DataUtil::saveIndexListASCII(outModeName, vector<int>(1, matchingMode))) return false;
	}

	return true;
}

bool PipelineGraphIO::error(string s) {
	cout << "Error: " << s << endl;
	return false;
}