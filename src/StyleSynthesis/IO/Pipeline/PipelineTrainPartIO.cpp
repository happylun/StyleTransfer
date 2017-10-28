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

#include "PipelineTrainPartIO.h"

#include <iostream>
#include <fstream>
#include <set>

#include "Mesh/MeshUtil.h"

#include "Segment/SegmentGroupApxCvx.h"

#include "Context/ContextPartGraph.h"
#include "Context/ContextPartGraphNodeGenerator.h"
#include "Context/ContextPartGraphTrain.h"

#include "Data/DataUtil.h"
#include "Data/StyleSynthesisConfig.h"

using namespace std;
using namespace StyleSynthesis;

bool PipelineTrainPartIO::process() {

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

	// compute element group, symmetry group, part pairs ...

	int numMeshes = (int)meshNameList.size();
	vector<string> meshNamesInScene(0);
	for (int meshID = 0; meshID < numMeshes; meshID++) {
		meshNamesInScene.push_back(meshNameList[meshID]);
		if (meshID + 1 == numMeshes || meshSceneList[meshID] != meshSceneList[meshID + 1]) {

			string sceneName = meshSceneList[meshID];
			cout << "Processing scene " << sceneName << endl;
			if (!processEachScene(sceneName, meshNamesInScene)) return false;

			meshNamesInScene.clear();

			//system("pause");
		}
	}

	// organize training data for functionality learning

	if (true) {
		vector<string> trainSceneList;
		string sceneListFileName = datasetPrefix + "mesh/train-scene-list.txt";
		ifstream sceneListFile(sceneListFileName);
		while (!sceneListFile.eof()) {
			string line;
			getline(sceneListFile, line);
			if (line.empty()) break;
			trainSceneList.push_back(line);
		}
		sceneListFile.close();
		string learnFolder = "learn/";
		if (!organizeTrainingData(meshNameList, meshSceneList, trainSceneList, learnFolder)) return false;
	}

	return true;
}

bool PipelineTrainPartIO::processEachScene(string sceneName, vector<string> &meshNames) {

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	// names

	string sceneFolder = datasetPrefix + "train/" + sceneName + "/";
	if (!FileUtil::makedir(sceneFolder)) return false;

	string sceneElementName = sceneFolder + "element.txt";
	string sceneElementVName = sceneFolder + "element.ply";
	string scenePartPairsName = sceneFolder + "partPairs.txt";
	string scenePartPairsVName = sceneFolder + "partPairs.ply";
	string sceneSymmetryVName = sceneFolder + "symmetry.ply";
	string sceneNodeSigmaName = sceneFolder + "nodeSigma.txt";
	string sceneEdgeSigmaName = sceneFolder + "edgeSigma.txt";

	if (FileUtil::existsfile(sceneNodeSigmaName)) return true; // early quit

	// load data

	int numMeshes = (int)meshNames.size();

	vector<TTriangleMesh> sceneMeshes(numMeshes);
	vector<vector<vector<vector<int>>>> sceneSegments(numMeshes); // 4-D vector... that's a bit scary...
	vector<ContextPartGraph> sceneGraphs(numMeshes);
	ContextPartGraphNodeGenerator sceneGraphNodeGenerator;

	for (int meshID = 0; meshID < numMeshes; meshID++) {

		//cout << "\rLoading data " << (meshID + 1) << " / " << numMeshes << "         ";

		// names

		string meshName = meshNames[meshID];
		string fullMeshName = datasetPrefix + "segment/" + meshName + ".ply";
		string segmentName = datasetPrefix + "segment/" + meshName + "-segment.txt";
		string graphFolder = datasetPrefix + "graph/" + meshName + "/";
		string graphHName = graphFolder + "graph-hierarchy.txt";
		string graphDName = graphFolder + "graph-descriptor.txt";
		string graphCName = graphFolder + "graph-context.txt";

		// data

		if (!MeshUtil::loadMesh(fullMeshName, sceneMeshes[meshID])) return false;
		if (!SegmentGroupApxCvx::loadSegments(segmentName, sceneSegments[meshID])) return false;

		if (!sceneGraphs[meshID].loadGraphHierarchy(graphHName, sceneGraphNodeGenerator, &sceneMeshes[meshID], &sceneSegments[meshID])) return false;
		if (!sceneGraphs[meshID].loadGraphDescriptor(graphDName)) return false;
		if (!sceneGraphs[meshID].loadGraphContext(graphCName)) return false;
	}
	//cout << endl;

	// algorithm

	ContextPartGraphTrain cpgt;
	if (!cpgt.loadSceneMesh(sceneMeshes)) return false;
	if (!cpgt.loadSceneSegments(sceneSegments)) return false;
	if (!cpgt.loadSceneGraph(sceneGraphs)) return false;
	if (!cpgt.process()) return false;
	if (!cpgt.visualizeElements(sceneElementVName)) return false;
	if (!cpgt.visualizePartPairs(scenePartPairsVName)) return false;
	if (!cpgt.visualizeSymmetries(sceneSymmetryVName)) return false;
	if (!cpgt.outputElements(sceneElementName)) return false;
	if (!cpgt.outputPartPairs(scenePartPairsName)) return false;
	if (!cpgt.outputSigma(sceneNodeSigmaName, sceneEdgeSigmaName)) return false;

	return true;
}

bool PipelineTrainPartIO::organizeTrainingData(
	vector<string> &meshNameList, vector<string> &meshSceneList,
	vector<string> &trainSceneList, string learnFolder)
{

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	//vector<string> trainSceneList;
	//string sceneListFileName = datasetPrefix + "mesh/train-scene-list.txt";
	//ifstream sceneListFile(sceneListFileName);
	//while (!sceneListFile.eof()) {
	//	string line;
	//	getline(sceneListFile, line);
	//	if (line.empty()) break;
	//	trainSceneList.push_back(line);
	//}
	//sceneListFile.close();

	string funcLearnFolder = datasetPrefix + learnFolder;
	if (!FileUtil::makedir(funcLearnFolder)) return false;

	string matchPairsName = funcLearnFolder + "matching_pairs.txt";
	if (FileUtil::existsfile(matchPairsName)) return true; // early quit

	string outNodeSigmaName = funcLearnFolder + "sigma-graph-node.txt";
	string outEdgeSigmaName = funcLearnFolder + "sigma-graph-edge.txt";

	string initWeightsFolder = datasetPrefix + "weights/real-init/";
	vector<string> initWeightsList;
	initWeightsList.push_back("weights-graph-node.txt");
	initWeightsList.push_back("weights-graph-edge.txt");
	initWeightsList.push_back("weights-graph-path.txt");
	initWeightsList.push_back("sigma-multipliers-graph-node.txt");
	initWeightsList.push_back("sigma-multipliers-graph-edge.txt");

	ofstream matchPairsFile(matchPairsName);

	int numMeshes = (int)meshNameList.size();
	int numTrainScene = (int)trainSceneList.size();

	vector<double> accumNodeSigma(0);
	vector<double> accumEdgeSigma(0);
	int dimNodeSigma = 0;
	int dimEdgeSigma = 0;

	int numEmptyScene = 0;

	for (int sceneID = 0; sceneID < numTrainScene; sceneID++) {
		string sceneName = trainSceneList[sceneID];
		vector<string> sceneMeshNames;
		for (int meshID = 0; meshID < numMeshes; meshID++) {
			if (meshSceneList[meshID] == sceneName) sceneMeshNames.push_back(meshNameList[meshID]);
		}

		cout << "Organizing scene " << sceneName << endl;

		// copy data

		string partPairsName = datasetPrefix + "train/" + sceneName + "/partPairs.txt";
		ifstream partPairsFile(partPairsName);

		int numMeshPairs;
		string line;
		getline(partPairsFile, line);
		numMeshPairs = stoi(line);
		for (int meshPairID = 0; meshPairID < numMeshPairs; meshPairID++) {

			getline(partPairsFile, line);
			vector<string> tokens;
			StringUtil::split(line, ' ', tokens);
			int srcMeshID = stoi(tokens[0]);
			int tgtMeshID = stoi(tokens[1]);
			string srcMeshName = sceneMeshNames[srcMeshID];
			string tgtMeshName = sceneMeshNames[tgtMeshID];

			string srcMeshShortName, tgtMeshShortName;
			if (true) {
				vector<string> srcTokens, tgtTokens;
				StringUtil::split(srcMeshName, '/', srcTokens);
				StringUtil::split(tgtMeshName, '/', tgtTokens);
				srcMeshShortName = srcTokens[1];
				tgtMeshShortName = tgtTokens[1];
			}

			matchPairsFile << srcMeshShortName << " " << tgtMeshShortName;
			for (int k = 2; k < (int)tokens.size(); k++) {
				matchPairsFile << " " << tokens[k];
			}
			matchPairsFile << endl;

			if (true) {
				string meshPairName = srcMeshShortName + "--" + tgtMeshShortName;

				string srcToFolder = funcLearnFolder + meshPairName + "/" + srcMeshShortName + "/";
				string tgtToFolder = funcLearnFolder + meshPairName + "/" + tgtMeshShortName + "/";
				if (!FileUtil::makedir(srcToFolder)) return false;
				if (!FileUtil::makedir(tgtToFolder)) return false;

				if (true) {
					// move parts

					string srcPartFromFolder = datasetPrefix + "scaled-part/" + srcMeshName + "/";
					string tgtPartFromFolder = datasetPrefix + "scaled-part/" + tgtMeshName + "/";


					int numSrcParts = 0;
					while (true) {
						string partName = to_string(numSrcParts) + ".ply";
						string fromFile = srcPartFromFolder + partName;
						string toFile = srcToFolder + partName;
						if (!FileUtil::existsfile(fromFile)) break;
						if (!FileUtil::copyfile(fromFile, toFile)) return false;
						numSrcParts++;
					}

					int numTgtParts = 0;
					while (true) {
						string partName = to_string(numTgtParts) + ".ply";
						string fromFile = tgtPartFromFolder + partName;
						string toFile = tgtToFolder + partName;
						if (!FileUtil::existsfile(fromFile)) break;
						if (!FileUtil::copyfile(fromFile, toFile)) return false;
						numTgtParts++;
					}
				}

				if (true) {
					// move raw mesh

					string srcRawMeshFromPrefix = datasetPrefix + "raw-mesh/" + srcMeshName;
					string tgtRawMeshFromPrefix = datasetPrefix + "raw-mesh/" + tgtMeshName;

					if (!FileUtil::copyfile(srcRawMeshFromPrefix + ".ply", srcToFolder + "sourceRawMesh.ply")) return false;
					if (!FileUtil::copyfile(tgtRawMeshFromPrefix + ".ply", tgtToFolder + "targetRawMesh.ply")) return false;
				}

				if (true) {
					// move segments

					string srcSegmentFromPrefix = datasetPrefix + "segment/" + srcMeshName;
					string tgtSegmentFromPrefix = datasetPrefix + "segment/" + tgtMeshName;

					if (!FileUtil::copyfile(srcSegmentFromPrefix + ".ply", srcToFolder + "sourceMesh.ply")) return false;
					if (!FileUtil::copyfile(srcSegmentFromPrefix + "-segment.txt", srcToFolder + "sourceSegment.txt")) return false;

					if (!FileUtil::copyfile(tgtSegmentFromPrefix + ".ply", tgtToFolder + "targetMesh.ply")) return false;
					if (!FileUtil::copyfile(tgtSegmentFromPrefix + "-segment.txt", tgtToFolder + "targetSegment.txt")) return false;
				}

				if (true) {
					// move graph

					string srcGraphFromFolder = datasetPrefix + "graph/" + srcMeshName + "/";
					string tgtGraphFromFolder = datasetPrefix + "graph/" + tgtMeshName + "/";

					if (!FileUtil::copyfile(srcGraphFromFolder + "graph-hierarchy.txt", srcToFolder + "sourceGraphHierarchy.txt")) return false;
					if (!FileUtil::copyfile(srcGraphFromFolder + "graph-descriptor.txt", srcToFolder + "sourceGraphDescriptor.txt")) return false;
					if (!FileUtil::copyfile(srcGraphFromFolder + "graph-context.txt", srcToFolder + "sourceGraphContext.txt")) return false;

					if (!FileUtil::copyfile(tgtGraphFromFolder + "graph-hierarchy.txt", tgtToFolder + "targetGraphHierarchy.txt")) return false;
					if (!FileUtil::copyfile(tgtGraphFromFolder + "graph-descriptor.txt", tgtToFolder + "targetGraphDescriptor.txt")) return false;
					if (!FileUtil::copyfile(tgtGraphFromFolder + "graph-context.txt", tgtToFolder + "targetGraphContext.txt")) return false;
				}
			}
		}

		partPairsFile.close();

		// accumulate sigma data

		if (true) {
			string nodeSigmaName = datasetPrefix + "train/" + sceneName + "/nodeSigma.txt";
			string edgeSigmaName = datasetPrefix + "train/" + sceneName + "/edgeSigma.txt";
			vector<double> nodeSigma;
			vector<double> edgeSigma;
			if (!DataUtil::loadValueListASCII(nodeSigmaName, nodeSigma)) return false;
			if (!DataUtil::loadValueListASCII(edgeSigmaName, edgeSigma)) return false;

			if (nodeSigma.empty() || edgeSigma.empty()) {
				numEmptyScene++;
			} else {
				if (accumNodeSigma.empty()) {
					accumNodeSigma = nodeSigma;
					dimNodeSigma = (int)nodeSigma.size();
				} else {
					for (int dim = 0; dim < dimNodeSigma; dim++) {
						accumNodeSigma[dim] += nodeSigma[dim];
					}
				}

				if (accumEdgeSigma.empty()) {
					accumEdgeSigma = edgeSigma;
					dimEdgeSigma = (int)edgeSigma.size();
				} else {
					for (int dim = 0; dim < dimEdgeSigma; dim++) {
						accumEdgeSigma[dim] += edgeSigma[dim];
					}
				}
			}
		}
	}

	matchPairsFile.close();

	cout << "Accumulating sigma..." << endl;

	if (true) {
		numTrainScene -= numEmptyScene;

		vector<double> outNodeSigma(dimNodeSigma, 0.0);
		vector<double> outEdgeSigma(dimEdgeSigma, 0.0);
		for (int dim = 0; dim < dimNodeSigma; dim++) {
			outNodeSigma[dim] = accumNodeSigma[dim] / numTrainScene;
		}
		for (int dim = 0; dim < dimEdgeSigma; dim++) {
			outEdgeSigma[dim] = accumEdgeSigma[dim] / numTrainScene;
		}
		if (!DataUtil::saveValueListASCII(outNodeSigmaName, outNodeSigma)) return false;
		if (!DataUtil::saveValueListASCII(outEdgeSigmaName, outEdgeSigma)) return false;
	}

	// copy initial weights

	for (string name : initWeightsList) {
		string fromFile = initWeightsFolder + name;
		string toFile = funcLearnFolder + name;
		if (!FileUtil::copyfile(fromFile, toFile)) return false;
	}

	return true;
}

bool PipelineTrainPartIO::error(string s) {
	cout << "Error: " << s << endl;
	return false;
}