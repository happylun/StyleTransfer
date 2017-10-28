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

#include "PipelineTrainValidationIO.h"

#include <iostream>
#include <fstream>
#include <set>

#include "Mesh/MeshUtil.h"

#include "Segment/SegmentGroupApxCvx.h"

#include "Context/ContextPartGraph.h"
#include "Context/ContextPartGraphNodeGenerator.h"
#include "Context/ContextPartGraphMatch.h"
#include "Context/ContextPartGraphCrossValidation.h"

#include "PipelineTrainPartIO.h"

#include "Data/DataUtil.h"
#include "Data/StyleSynthesisConfig.h"

using namespace std;
using namespace StyleSynthesis;

bool PipelineTrainValidationIO::process() {

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	vector<string> meshNameList(0);
	vector<string> meshSceneList(0);

	string meshListFileName = datasetPrefix + "mesh/mesh-list-all.txt";
	ifstream meshListFile(meshListFileName);
	while (!meshListFile.eof()) {
		string line;
		getline(meshListFile, line);
		if (line.empty()) break;

		string meshName = StringUtil::trim(line);
		vector<string> tokens;
		StringUtil::split(meshName, '/', tokens);
		string sceneName = tokens[0];

		meshNameList.push_back(meshName);
		meshSceneList.push_back(sceneName);
	}
	meshListFile.close();

	// functionality validation

	if (true) {

		int numCrossValidationSplits = StyleSynthesisConfig::mContext_FunctionalityLearningCrossValidationSplits;
		vector<string> allSceneList;
		for (string sceneName : meshSceneList) {
			if (allSceneList.empty() || sceneName != allSceneList.back()) {
				allSceneList.push_back(sceneName);
			}
		}
		int numAllScenes = (int)allSceneList.size();
		for (int split = 0; split < numCrossValidationSplits; split++) {
			vector<string> trainSceneList;
			for (int sceneID = (numAllScenes * split / numCrossValidationSplits);
				sceneID < (numAllScenes * (split + 1) / numCrossValidationSplits);
				sceneID++)
			{
				trainSceneList.push_back(allSceneList[sceneID]);
			}
			string learnFolder = "learn-CV" + to_string(split) + "/";

			if (!PipelineTrainPartIO::organizeTrainingData(meshNameList, meshSceneList, trainSceneList, learnFolder)) return false;
		}

		ContextPartGraphCrossValidation cpgcv;
		if (!cpgcv.train(numCrossValidationSplits)) return false;
		if (!cpgcv.test(numCrossValidationSplits)) return false;
		if (!cpgcv.evaluate()) return false;
		if (!cpgcv.compare()) return false;
	}

	// match all pairs of shapes (for auto target selection)
	//if (!matchAllPairShapes(meshNameList)) return false;

	return true;
}

bool PipelineTrainValidationIO::matchAllPairShapes(vector<string> &meshNameList) {

	// compute all pair similarity

	int numMeshes = (int)meshNameList.size();

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string outMatrixName = datasetPrefix + "train/all-pair-similarity.txt";
	string outFinalMatrixName = datasetPrefix + "train/all-pair-similarity-final.txt";
	string exemplarFileName = datasetPrefix + "train/exemplar-list.txt";

	vector<int> labelList;
	if (!DataUtil::loadIndexListASCII(datasetPrefix + "mesh/mesh-labels-all.txt", labelList)) return false;

	set<string> exemplarSet;
	ifstream exemplarListFile(exemplarFileName);
	while (!exemplarListFile.eof()) {
		string line;
		getline(exemplarListFile, line);
		if (line.empty()) break;

		string meshName = StringUtil::trim(line);
		exemplarSet.insert(meshName);
	}
	exemplarListFile.close();

	Eigen::MatrixXd outMatrix;
	if (!FileUtil::existsfile(outMatrixName)) {
		outMatrix.resize(numMeshes, numMeshes);
		outMatrix.setIdentity();
		if (!DataUtil::saveMatrixASCII(outMatrixName, outMatrix)) return false;
	} else {
		if (!DataUtil::loadMatrixASCII(outMatrixName, outMatrix)) return false;
	}

	for (int srcID = 0; srcID < numMeshes; srcID++) {
		if (exemplarSet.find(meshNameList[srcID]) == exemplarSet.end()) continue;
		for (int tgtID = 0; tgtID < numMeshes; tgtID++) {

			//if (outMatrix(srcID, tgtID) != 0) continue; // already computed
			//if (labelList[tgtID] != 3) continue;

			string srcName = meshNameList[srcID];
			string tgtName = meshNameList[tgtID];
			cout << "=============== Comparing mesh pair " << srcName << " -- " << tgtName << " ===============" << endl;

			// names

			string sourceMeshName = datasetPrefix + "segment/" + srcName + ".ply";
			string targetMeshName = datasetPrefix + "segment/" + tgtName+ ".ply";

			string sourceSegmentName = datasetPrefix + "segment/" + srcName+ "-segment.txt";
			string targetSegmentName = datasetPrefix + "segment/" + tgtName+ "-segment.txt";

			string sourceGraphFolder = datasetPrefix + "graph/" + srcName + "/";
			string targetGraphFolder = datasetPrefix + "graph/" + tgtName + "/";

			string sourceGraphHName = sourceGraphFolder + "graph-hierarchy.txt";
			string sourceGraphDName = sourceGraphFolder + "graph-descriptor.txt";
			string sourceGraphCName = sourceGraphFolder + "graph-context.txt";
			string sourceGraphSimName = sourceGraphFolder + "graph-self-similarity.txt";
			string sourceModeName = sourceGraphFolder + "graph-self-mode.txt";

			string targetGraphHName = targetGraphFolder + "graph-hierarchy.txt";
			string targetGraphDName = targetGraphFolder + "graph-descriptor.txt";
			string targetGraphCName = targetGraphFolder + "graph-context.txt";
			string targetGraphSimName = targetGraphFolder + "graph-self-similarity.txt";
			string targetModeName = targetGraphFolder + "graph-self-mode.txt";

			string weightsFolder = datasetPrefix + "weights/";

			// data

			TTriangleMesh sourceMesh;
			TTriangleMesh targetMesh;
			vector<vector<vector<int>>> sourceSegments;
			vector<vector<vector<int>>> targetSegments;

			if (!MeshUtil::loadMesh(sourceMeshName, sourceMesh)) return false;
			if (!MeshUtil::loadMesh(targetMeshName, targetMesh)) return false;

			if (!SegmentGroupApxCvx::loadSegments(sourceSegmentName, sourceSegments)) return false;
			if (!SegmentGroupApxCvx::loadSegments(targetSegmentName, targetSegments)) return false;

			ContextPartGraphNodeGenerator cpgng;

			ContextPartGraph sourceGraph;
			if (!sourceGraph.loadGraphHierarchy(sourceGraphHName, cpgng, &sourceMesh, &sourceSegments)) return false;
			if (!sourceGraph.loadGraphDescriptor(sourceGraphDName)) return false;
			if (!sourceGraph.loadGraphContext(sourceGraphCName)) return false;

			ContextPartGraph targetGraph;
			if (!targetGraph.loadGraphHierarchy(targetGraphHName, cpgng, &targetMesh, &targetSegments)) return false;
			if (!targetGraph.loadGraphDescriptor(targetGraphDName)) return false;
			if (!targetGraph.loadGraphContext(targetGraphCName)) return false;

			Eigen::MatrixXd sourceSimilarity;
			Eigen::MatrixXd targetSimilarity;
			if (!DataUtil::loadMatrixBinary(sourceGraphSimName, sourceSimilarity)) return false;
			if (!DataUtil::loadMatrixBinary(targetGraphSimName, targetSimilarity)) return false;

			int sourceMatchingMode;
			int targetMatchingMode;
			if (true) {
				vector<int> modes;
				if (!DataUtil::loadIndexListASCII(sourceModeName, modes)) return false;
				sourceMatchingMode = modes[0];
				if (!DataUtil::loadIndexListASCII(targetModeName, modes)) return false;
				targetMatchingMode = modes[0];
			}
			

			// matching

			Eigen::MatrixXd similarityMatrix;
			ContextPartGraphMatch cpgm;
			if (!cpgm.loadWeights(weightsFolder)) return false;
			if (!cpgm.loadGraph(sourceGraph, targetGraph)) return false;
			if (!cpgm.loadMatchingMode(max(sourceMatchingMode, targetMatchingMode))) return false;
			if (!cpgm.process()) return false;
			if (!cpgm.normalizeGraphSimilarity(sourceSimilarity, targetSimilarity)) return false;
			double similarity;
			if (!cpgm.exportGraphSimilarity(similarity)) return false;

			outMatrix(srcID, tgtID) = similarity;
			outMatrix(tgtID, srcID) = similarity;
		}
		//if (!DataUtil::saveMatrixASCII(outMatrixName, outMatrix)) return false;

		Eigen::VectorXd outVec = outMatrix.col(srcID);
		if (!DataUtil::saveVectorASCII(outMatrixName, outVec)) return false;
	}

	if (!FileUtil::copyfile(outMatrixName, outFinalMatrixName)) return false;

	return true;
}

bool PipelineTrainValidationIO::error(string s) {
	cout << "Error: " << s << endl;
	return false;
}