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

#include "PipelineTrainCurveIO.h"

#include <iostream>
#include <fstream>
#include <unordered_set>

#include "Context/ContextPartGraphTrainCurve.h"
#include "Context/ContextPartGraphLearningCurve.h"

#include "Mesh/MeshUtil.h"
#include "Curve/CurveUtil.h"

#include "Data/DataUtil.h"
#include "Data/StyleSynthesisConfig.h"

using namespace std;
using namespace StyleSynthesis;

bool PipelineTrainCurveIO::process() {

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	vector<string> meshNameList(0);
	vector<string> meshSceneList(0);
	vector<string> trainSceneList(0);
	map<string, int> meshLabelMap;
	int numLabels;

	if (true) {
		// load meshes / scenes
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
	}

	if (true) {
		// load training scenes
		string sceneListFileName = datasetPrefix + "mesh/train-scene-list.txt";
		ifstream sceneListFile(sceneListFileName);
		while (!sceneListFile.eof()) {
			string line;
			getline(sceneListFile, line);
			if (line.empty()) break;
			trainSceneList.push_back(line);
		}
		sceneListFile.close();
	}

	if (true) {
		// build label map
		vector<int> labelList;
		string meshLabelsFileName = datasetPrefix + "mesh/mesh-labels-all.txt";
		if (!DataUtil::loadIndexListASCII(meshLabelsFileName, labelList)) return false;
		if (labelList.size() != meshNameList.size()) return error("incorrect label list file");
		meshLabelMap.clear();
		unordered_set<int> labelSet;
		for (int meshID = 0; meshID < (int)meshNameList.size(); meshID++) {
			meshLabelMap[meshNameList[meshID]] = labelList[meshID] - 1; // labels are 1-indexed
			labelSet.insert(labelList[meshID]);
		}
		numLabels = (int)labelSet.size();
		cout << "Loaded " << numLabels << " labels" << endl;
	}

	// compute curve pairs for each scene...

	int numMeshes = (int)meshNameList.size();
	vector<string> meshNamesInScene(0);
	for (int meshID = 0; meshID < numMeshes; meshID++) {
		meshNamesInScene.push_back(meshNameList[meshID]);
		if (meshID + 1 == numMeshes || meshSceneList[meshID] != meshSceneList[meshID + 1]) {

			string sceneName = meshSceneList[meshID];
			cout << "Processing scene " << sceneName << endl;

			if (!processEachScene(sceneName, meshNamesInScene)) return false;

			meshNamesInScene.clear();
		}
	}

	// parameter learning for curve weights

	if (true) {		
		int numScenes = (int)trainSceneList.size();
		vector<vector<string>> trainSceneMeshes(numScenes);
		map<string, int> sceneMap;
		for (int sceneID = 0; sceneID < numScenes; sceneID++) {
			sceneMap[trainSceneList[sceneID]] = sceneID;
		}

		for (int meshID = 0; meshID < numMeshes; meshID++) {
			auto it = sceneMap.find(meshSceneList[meshID]);
			if (it != sceneMap.end()) {
				trainSceneMeshes[it->second].push_back(meshNameList[meshID]);
			}
		}

		if (!learnParameters(trainSceneList, trainSceneMeshes, meshLabelMap, numLabels)) return false;
	}

	return true;
}

bool PipelineTrainCurveIO::processEachScene(string sceneName, vector<string> &meshNames) {

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	// names

	string sceneFolder = datasetPrefix + "train-curve/" + sceneName + "/";
	if (!FileUtil::makedir(sceneFolder)) return false;

	string sceneCurvePairsName = sceneFolder + "curvePairs.txt";
	string sceneCurvePairsVName = sceneFolder + "curvePairs.ply";

	string curveViewPointName = datasetPrefix + "curve/viewpoint.txt";

	if (FileUtil::existsfile(sceneCurvePairsName)) return true; // early quit

	// prepare data

	int numMeshes = (int)meshNames.size();

	vector<string> sceneCurveFolders(numMeshes);
	vector<TTriangleMesh> sceneMeshes(numMeshes);

	for (int meshID = 0; meshID < numMeshes; meshID++) {
		string meshName = meshNames[meshID];
		string fullMeshName = datasetPrefix + "segment/" + meshName + ".ply";
		string curveFolder = datasetPrefix + "curve/" + meshName + "/";

		if (!MeshUtil::loadMesh(fullMeshName, sceneMeshes[meshID])) return false;
		sceneCurveFolders[meshID] = curveFolder;
	}

	vector<vec3> curveViewPoints;
	if (!CurveUtil::loadViewPoints(curveViewPointName, curveViewPoints)) return false;

	// algorithm

	ContextPartGraphTrainCurve cpgtc;
	if (!cpgtc.loadSceneCurves(sceneCurveFolders, curveViewPoints)) return false;
	if (!cpgtc.process()) return false;
	if (!cpgtc.visualizeCurvePairs(sceneCurvePairsVName, sceneMeshes)) return false;
	if (!cpgtc.outputCurvePairs(sceneCurvePairsName)) return false;

	return true;
}

bool PipelineTrainCurveIO::learnParameters(
	vector<string> &trainScenes,
	vector<vector<string>> &trainSceneMeshes,
	map<string, int> &meshLabelMap,
	int numLabels)
{

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;
	string trainFolder = datasetPrefix + "CurveLearning/";
	if (!FileUtil::makedir(trainFolder)) return false;

	for (int labelID1 = 0; labelID1 < numLabels; labelID1++) {
		for (int labelID2 = labelID1; labelID2 < numLabels; labelID2++) {

			string labelPairString = to_string(labelID1) + "--" + to_string(labelID2);
			string labelPairFolder = trainFolder + labelPairString + "/";
			cout << "||||||||||||||||| " << labelPairString << " |||||||||||||||||" << endl;

			string outSigmasName = labelPairFolder + "sigmas.txt";
			string outWeightsName = labelPairFolder + "weights.txt";
			string outCountsName = labelPairFolder + "counts.txt";

			if (FileUtil::existsfile(outWeightsName)) continue; // early quit

			vec2i labelPair(labelID1, labelID2);

			ContextPartGraphLearningCurve cpglc;
			int numMeshPairs, numCurvePairs;
			if (!cpglc.loadData(
				trainScenes, trainSceneMeshes,
				meshLabelMap, labelPair,
				numMeshPairs, numCurvePairs)) return false;
			if (numCurvePairs <= 0) continue;

			if (!FileUtil::makedir(labelPairFolder)) return false;
			if (!cpglc.process()) return false;
			if (!cpglc.outputSigmas(outSigmasName)) return false;
			if (!cpglc.outputWeights(outWeightsName)) return false;

			if (true) {
				vector<int> counts(0);
				counts.push_back(numMeshPairs);
				counts.push_back(numCurvePairs);
				if (!DataUtil::saveIndexListASCII(outCountsName, counts)) return false;
			}
		}
	}

	return true;
}

bool PipelineTrainCurveIO::error(string s) {
	cout << "Error: " << s << endl;
	return false;
}