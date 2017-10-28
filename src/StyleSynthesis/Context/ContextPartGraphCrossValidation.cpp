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

#include "ContextPartGraphCrossValidation.h"

#include <iostream>
#include <fstream>
#include <set>

#include "IO/Pipeline/PipelineTrainLearningIO.h"

#include "Mesh/MeshUtil.h"

#include "Segment/SegmentGroupApxCvx.h"

#include "ContextPartGraph.h"
#include "ContextPartGraphNodeGenerator.h"
#include "ContextPartGraphMatch.h"

#include "Utility/FileUtil.h"
#include "Utility/PlyExporter.h"

#include "Data/DataUtil.h"
#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

ContextPartGraphCrossValidation::ContextPartGraphCrossValidation() {
}

ContextPartGraphCrossValidation::~ContextPartGraphCrossValidation() {
}

bool ContextPartGraphCrossValidation::train(int numSplits) {
	
	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	//for (int splitID = numSplits-1; splitID >= 0; splitID--) { // do it reversely
	for (int splitID = 0; splitID < numSplits; splitID++) {

		string learnFolderName = "learn-CV" + to_string(splitID) + "/";
		string learnFolderPath = datasetPrefix + learnFolderName;

		string learnedWeightsFolder = datasetPrefix + "weights-CV" + to_string(splitID) + "/";
		if (FileUtil::existsfile(learnedWeightsFolder + "weights-graph-node.txt")) continue; // early quit
		if (!FileUtil::makedir(learnedWeightsFolder)) return false;

		string learnedWeightsName = learnFolderPath + "weights.txt";
		string initialNodeWeightName = learnFolderPath + "weights-graph-node.txt";
		string initialEdgeWeightName = learnFolderPath + "weights-graph-edge.txt";
		string initialPathWeightName = learnFolderPath + "weights-graph-path.txt";
		string initialNodeSigmaMultName = learnFolderPath + "sigma-multipliers-graph-node.txt";
		string initialEdgeSigmaMultName = learnFolderPath + "sigma-multipliers-graph-edge.txt";

		if (StyleSynthesisConfig::mContext_UseLagaDescriptors) {

			// copy weights
			if (!FileUtil::copyfile(initialNodeWeightName, learnedWeightsFolder + "weights-graph-node.txt")) return false;
			if (!FileUtil::copyfile(initialEdgeWeightName, learnedWeightsFolder + "weights-graph-edge.txt")) return false;
			if (!FileUtil::copyfile(initialPathWeightName, learnedWeightsFolder + "weights-graph-path.txt")) return false;
			if (!FileUtil::copyfile(initialNodeSigmaMultName, learnedWeightsFolder + "sigma-multipliers-graph-node.txt")) return false;
			if (!FileUtil::copyfile(initialEdgeSigmaMultName, learnedWeightsFolder + "sigma-multipliers-graph-edge.txt")) return false;
			if (!FileUtil::copyfile(learnFolderPath + "sigma-graph-node.txt", learnedWeightsFolder + "sigma-graph-node.txt")) return false;
			if (!FileUtil::copyfile(learnFolderPath + "sigma-graph-edge.txt", learnedWeightsFolder + "sigma-graph-edge.txt")) return false;

			continue; // no need to train
		}

		// do the learning

		if (!FileUtil::existsfile(learnedWeightsName)) {
			if (!PipelineTrainLearningIO::process(learnFolderName)) return false;
		}

		// determine dimension of weights

		vector<double> initNodeWeight, initEdgeWeight, initPathWeight;
		vector<double> initNodeSigmaMult, initEdgeSigmaMult;
		if (!DataUtil::loadValueListASCII(initialNodeWeightName, initNodeWeight)) return false;
		if (!DataUtil::loadValueListASCII(initialEdgeWeightName, initEdgeWeight)) return false;
		if (!DataUtil::loadValueListASCII(initialPathWeightName, initPathWeight)) return false;
		if (!DataUtil::loadValueListASCII(initialNodeSigmaMultName, initNodeSigmaMult)) return false;
		if (!DataUtil::loadValueListASCII(initialEdgeSigmaMultName, initEdgeSigmaMult)) return false;
		int dimNodeWeight = (int)initNodeWeight.size();
		int dimEdgeWeight = (int)initEdgeWeight.size();
		int dimPathWeight = (int)initPathWeight.size();
		int dimNodeSigmaMult = (int)initNodeSigmaMult.size();
		int dimEdgeSigmaMult = (int)initEdgeSigmaMult.size();

		// organize resulting weights

		vector<double> learnedWeights;
		if (!DataUtil::loadValueListASCII(learnedWeightsName, learnedWeights)) return false;
		if (dimNodeWeight + dimNodeSigmaMult + dimEdgeWeight + dimEdgeSigmaMult + dimPathWeight != (int)learnedWeights.size()) {
			cout << "Error: incorrect dimension of learned weights" << endl;
			return false;
		}
		int from, to;
		from = 0; to = dimNodeWeight;
		vector<double> learnedNodeWeight(learnedWeights.begin() + from, learnedWeights.begin() + to);
		from = to; to += dimNodeSigmaMult;
		vector<double> learnedNodeSigmaMult(learnedWeights.begin() + from, learnedWeights.begin() + to);
		from = to; to += dimEdgeWeight;
		vector<double> learnedEdgeWeight(learnedWeights.begin() + from, learnedWeights.begin() + to);
		from = to; to += dimEdgeSigmaMult;
		vector<double> learnedEdgeSigmaMult(learnedWeights.begin() + from, learnedWeights.begin() + to);
		from = to; to += dimPathWeight;
		vector<double> learnedPathWeight(learnedWeights.begin() + from, learnedWeights.begin() + to);

		if (!DataUtil::saveValueListASCII(learnedWeightsFolder + "weights-graph-node.txt", learnedNodeWeight)) return false;
		if (!DataUtil::saveValueListASCII(learnedWeightsFolder + "weights-graph-edge.txt", learnedEdgeWeight)) return false;
		if (!DataUtil::saveValueListASCII(learnedWeightsFolder + "weights-graph-path.txt", learnedPathWeight)) return false;
		if (!DataUtil::saveValueListASCII(learnedWeightsFolder + "sigma-multipliers-graph-node.txt", learnedNodeWeight)) return false;
		if (!DataUtil::saveValueListASCII(learnedWeightsFolder + "sigma-multipliers-graph-edge.txt", learnedEdgeWeight)) return false;
		if (!FileUtil::copyfile(learnFolderPath + "sigma-graph-node.txt", learnedWeightsFolder + "sigma-graph-node.txt")) return false;
		if (!FileUtil::copyfile(learnFolderPath + "sigma-graph-edge.txt", learnedWeightsFolder + "sigma-graph-edge.txt")) return false;
	}

	return true;
}

bool ContextPartGraphCrossValidation::test(int numSplits) {

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	vector<vec3i> splitResults;

	for (int splitID = 0; splitID < numSplits; splitID++) {
		string learnFolder = datasetPrefix + "learn-CV" + to_string(splitID) + "/";
		string weightsFolder = datasetPrefix + "weights-CV" + to_string(splitID) + "/";

		int numTrainPairs = 0;
		if (true) {
			vector<string> currentSourceList, currentTargetList;
			vector<vector<vec2i>> currentPairList;
			if (!loadMatchingPair(
				learnFolder + "matching_pairs.txt",
				currentSourceList,
				currentTargetList,
				currentPairList)) return false;
			numTrainPairs = (int)currentSourceList.size();
		}

		// collect all part pairs from other splits
		vector<string> sourceList, targetList, meshPairList;
		vector<vector<vec2i>> nodePairList;
		for (int otherID = 0; otherID < numSplits; otherID++) {
			if (otherID == splitID) continue;
			string otherFolder = datasetPrefix + "learn-CV" + to_string(otherID) + "/";

			vector<string> currentSourceList, currentTargetList;
			vector<vector<vec2i>> currentPairList;
			if (!loadMatchingPair(
				otherFolder + "matching_pairs.txt",
				currentSourceList,
				currentTargetList,
				currentPairList)) return false;
			int numMeshPairs = (int)currentSourceList.size();

			for (int meshPairID = 0; meshPairID < numMeshPairs; meshPairID++) {
				string sourceName = currentSourceList[meshPairID];
				string targetName = currentTargetList[meshPairID];
				string meshPairName = sourceName + "--" + targetName;
				string sourceFolder = otherFolder + meshPairName + "/" + sourceName + "/";
				string targetFolder = otherFolder + meshPairName + "/" + targetName + "/";

				sourceList.push_back(sourceFolder);
				targetList.push_back(targetFolder);
				meshPairList.push_back(meshPairName);
				nodePairList.push_back(currentPairList[meshPairID]);
			}
		}

		// test each pair
		int numMeshPairs = (int)sourceList.size();
		int numTestedMeshPairs = 0;
		int numAllNodePairs = 0;
		int numAllCorrectPairs = 0;
#pragma omp parallel for
		for (int pairID = 0; pairID < numMeshPairs; pairID++) {
			//if (meshPairList[pairID] != "L37-1--L37-5") continue; // HACK: debug
			int numNodePairs = (int)nodePairList[pairID].size();
			int numCorrectPairs = 0;
			if (!validatePair(
				sourceList[pairID],
				targetList[pairID],
				weightsFolder,
				nodePairList[pairID],
				numCorrectPairs))
			{
				cout << "Error: validate pairs" << endl;
				system("pause");
			}
#pragma omp critical
			{
				numTestedMeshPairs++;
				numAllNodePairs += numNodePairs;
				numAllCorrectPairs += numCorrectPairs;
				cout << "Testing pair " << meshPairList[pairID] << " : " << numTestedMeshPairs << " / " << numMeshPairs << " : ";
				cout << numAllCorrectPairs << " / " << numAllNodePairs << " = ";
				cout << (numAllCorrectPairs / (double)numAllNodePairs) << endl;
			}
		}

		vec3i result = vec3i(numTrainPairs, numAllCorrectPairs, numAllNodePairs);
		splitResults.push_back(result);
		cout << "Split result: " << result << " = " << (result[1] / (double)result[2]) << endl;
	}

	cout << "Overall results:" << endl;
	int overallCorrectPairs = 0;
	int overallAllPairs = 0;
	for (int splitID = 0; splitID < numSplits; splitID++) {
		vec3i result = splitResults[splitID];
		cout << "\t Split " << splitID << ": " << result << " = " << (result[1] / (double)result[2]) << endl;
		overallCorrectPairs += result[1];
		overallAllPairs += result[2];
	}
	cout << "Overall: " << overallCorrectPairs << " / " << overallAllPairs << " = " << (overallCorrectPairs / (double)overallAllPairs) << endl;

	return true;
}

bool ContextPartGraphCrossValidation::evaluate() {

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string matchPairFileName = datasetPrefix + "validation/matching-pairs.txt";
	if (!FileUtil::existsfile(matchPairFileName)) return true;

	string weightsFolder = datasetPrefix + "weights/" +
		(StyleSynthesisConfig::mContext_UseLagaDescriptors ? "laga/" : "original/");

	vector<string> sourceList, targetList;
	vector<vector<vec2i>> pairList;
	if (!loadMatchingPair(matchPairFileName, sourceList, targetList, pairList)) return false;

	int numAllNodePairs = 0;
	int numAllCorrectPairs = 0;

	ContextPartGraphNodeGenerator nodeGen;

	int numMeshPairs = (int)sourceList.size();
	for (int meshPairID = 0; meshPairID < numMeshPairs; meshPairID++) {
		string sourceName = sourceList[meshPairID];
		string targetName = targetList[meshPairID];
		vector<vec2i> &nodePairList = pairList[meshPairID];

		TTriangleMesh sourceMesh, targetMesh;
		vector<vector<vector<int>>> sourceSegment, targetSegment;
		ContextPartGraph sourceGraph, targetGraph;
		if (!loadMatchingData(
			sourceName, targetName, nodeGen, sourceMesh, targetMesh,
			sourceSegment, targetSegment, sourceGraph, targetGraph)) return false;

		// validation...

		int numNodePairs = (int)nodePairList.size();
		int numCorrectPairs = 0;
		if (!validatePair(&sourceGraph, &targetGraph, weightsFolder, nodePairList, numCorrectPairs)) return false;

		numAllNodePairs += numNodePairs;
		numAllCorrectPairs += numCorrectPairs;
		cout << "Checking mesh pair " << (meshPairID + 1) << " / " << numMeshPairs << ": ";
		cout << numAllCorrectPairs << " / " << numAllNodePairs;
		cout << " = " << (numAllCorrectPairs / (double)numAllNodePairs) << endl;
	}

	cout << "Final accuracy = " << numAllCorrectPairs << " / " << numAllNodePairs;
	cout << " = " << (numAllCorrectPairs / (double)numAllNodePairs) << endl;

	return true;
}

bool ContextPartGraphCrossValidation::compare() {

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;
	string outTripletFileName = datasetPrefix + "validation/functionality-triplets.txt";
	string outAllTripletFileName = datasetPrefix + "validation/all-functionality-triplets.txt";

	bool filterFromUserData = false;

	vector<string> sourceList, targetList;
	vector<vector<vec2i>> dataPairList;
	if (filterFromUserData) {
		string matchPairFileName = datasetPrefix + "validation/matching-pairs.txt";
		// dataPairList is match pair list
		if (!loadMatchingPair(matchPairFileName, sourceList, targetList, dataPairList)) return false;
	} else {
		string meshPairFileName = datasetPrefix + "mesh/functionality-evaluation-list.txt";
		// dataPairList is mesh pair list
		if (!loadMatchingPair(meshPairFileName, sourceList, targetList, dataPairList)) return false;
	}

	int numTotalTriplets = 0;
	int numTriplets = 0;
	vec3i overallStats(0, 0, 0);
	ofstream outTripletFile(outTripletFileName);
	ofstream outAllTripletFile(outAllTripletFileName);

	ContextPartGraphNodeGenerator nodeGen;

	int numMeshPairs = (int)sourceList.size();
	for (int meshPairID = 0; meshPairID < numMeshPairs; meshPairID++) {

		cout << "Checking mesh pair " << (meshPairID + 1) << " / " << numMeshPairs << endl;

		string sourceName = sourceList[meshPairID];
		string targetName = targetList[meshPairID];

		TTriangleMesh sourceMesh, targetMesh;
		vector<vector<vector<int>>> sourceSegment, targetSegment;
		ContextPartGraph sourceGraph, targetGraph;
		if (!loadMatchingData(
			sourceName, targetName, nodeGen, sourceMesh, targetMesh,
			sourceSegment, targetSegment, sourceGraph, targetGraph)) return false;

		vector<int> validSourceNodes(0);
		if (filterFromUserData) {
			vector<vec2i> &nodePairList = dataPairList[meshPairID];
			for (vec2i &nodePair : nodePairList) validSourceNodes.push_back(nodePair[0]);
		} else {
			int numSourceNodes = (int)sourceGraph.mRootNode->mChildren.size();
			for (int nodeID = 0; nodeID < numSourceNodes; nodeID++) validSourceNodes.push_back(nodeID);
		}

		if (true) {
			bool oldFlag = StyleSynthesisConfig::mContext_UseLagaDescriptors;
			string weightsFolder;

			cout << "Matching using our method..." << endl;
			StyleSynthesisConfig::mContext_UseLagaDescriptors = false;
			weightsFolder = datasetPrefix + "weights/original/";
			Eigen::MatrixXd simMatOurs;
			if (!computeSimilarity(&sourceGraph, &targetGraph, weightsFolder, simMatOurs)) return false;

			cout << "Matching using Laga et al.'s method..." << endl;
			StyleSynthesisConfig::mContext_UseLagaDescriptors = true;
			weightsFolder = datasetPrefix + "weights/laga/";
			Eigen::MatrixXd simMatLaga;
			if (!computeSimilarity(&sourceGraph, &targetGraph, weightsFolder, simMatLaga)) return false;

			StyleSynthesisConfig::mContext_UseLagaDescriptors = oldFlag;

			cout << "Finding triplets..." << endl;
			vector<vec3i> validTriplets, diffTriplets;
			vec3i stats;
			if (!getDifferentMatching(&sourceGraph, &targetGraph, simMatOurs, simMatLaga,
				validSourceNodes, validTriplets, diffTriplets, stats)) return false;
			overallStats += stats;

			numTotalTriplets += (int)validTriplets.size();
			for (vec3i &triplet : validTriplets) {
				outAllTripletFile << sourceName << " " << targetName << " " << triplet << endl;
			}

			if (!diffTriplets.empty()) {
				for (vec3i &triplet : diffTriplets) {
					outTripletFile << sourceName << " " << targetName << " " << triplet << endl;
					numTriplets++;

					if (true) {
						// visualization

						string visFolder = datasetPrefix + "validation/triplet-" + to_string(numTriplets) + "/";
						if (!FileUtil::makedir(visFolder)) return false;

						auto srcNode = sourceGraph.mAllNodes[triplet[0]];
						vector<vec3i> srcColors(sourceMesh.indices.size(), vec3i(127, 127, 127));
						for (int faceID : sourceSegment[srcNode->mPartLevelID][srcNode->mPartSegmentID]) {
							srcColors[faceID] = vec3i(255, 0, 0);
						}
						auto tgtNodeOurs = targetGraph.mAllNodes[triplet[1]];
						vector<vec3i> tgtColorsOurs(targetMesh.indices.size(), vec3i(127, 127, 127));
						for (int faceID : targetSegment[tgtNodeOurs->mPartLevelID][tgtNodeOurs->mPartSegmentID]) {
							tgtColorsOurs[faceID] = vec3i(255, 0, 0);
						}
						auto tgtNodeLaga = targetGraph.mAllNodes[triplet[2]];
						vector<vec3i> tgtColorsLaga(targetMesh.indices.size(), vec3i(127, 127, 127));
						for (int faceID : targetSegment[tgtNodeLaga->mPartLevelID][tgtNodeLaga->mPartSegmentID]) {
							tgtColorsLaga[faceID] = vec3i(255, 0, 0);
						}
						PlyExporter peA, peB, peC, peAll;

						if (!peA.addMesh(&sourceMesh.positions, &sourceMesh.normals, &sourceMesh.indices, &srcColors)) return false;
						if (!peAll.addMesh(&sourceMesh.positions, &sourceMesh.normals, &sourceMesh.indices, &srcColors, vec3(0.0f, 3.0f, 0.0f))) return false;

						if (!peB.addMesh(&targetMesh.positions, &targetMesh.normals, &targetMesh.indices, &tgtColorsOurs)) return false;
						if (!peAll.addMesh(&targetMesh.positions, &targetMesh.normals, &targetMesh.indices, &tgtColorsOurs, vec3(-2.0f, 0.0f, 0.0f))) return false;

						if (!peC.addMesh(&targetMesh.positions, &targetMesh.normals, &targetMesh.indices, &tgtColorsLaga)) return false;
						if (!peAll.addMesh(&targetMesh.positions, &targetMesh.normals, &targetMesh.indices, &tgtColorsLaga, vec3(2.0f, 0.0f, 0.0f))) return false;

						if (!peA.output(visFolder + "a.ply")) return false;
						if (!peB.output(visFolder + "b.ply")) return false;
						if (!peC.output(visFolder + "c.ply")) return false;
						if (!peAll.output(visFolder + "triplet.ply")) return false;
					}
				}
			}
		}
	}

	outTripletFile.close();
	outAllTripletFile.close();

	cout << "Number of triplets: " << numTriplets << " / " << numTotalTriplets << endl;
	cout << "Overall statistics: " << overallStats << endl;

	return true;
}

bool ContextPartGraphCrossValidation::validatePair(
	string sourceFolder,
	string targetFolder,
	string weightsFolder,
	vector<vec2i> &nodePairs,
	int &numCorrectPairs)
{

	// names...

	string sourceMeshName = sourceFolder + "sourceMesh.ply";
	string sourceSegmentName = sourceFolder + "sourceSegment.txt";
	string sourceGraphHName = sourceFolder + "sourceGraphHierarchy.txt";
	string sourceGraphDName = sourceFolder + "sourceGraphDescriptor.txt";
	string sourceGraphCName = sourceFolder + "sourceGraphContext.txt";

	string targetMeshName = targetFolder + "targetMesh.ply";
	string targetSegmentName = targetFolder + "targetSegment.txt";
	string targetGraphHName = targetFolder + "targetGraphHierarchy.txt";
	string targetGraphDName = targetFolder + "targetGraphDescriptor.txt";
	string targetGraphCName = targetFolder + "targetGraphContext.txt";

	// data...

	TTriangleMesh sourceMesh, targetMesh;
	vector<vector<vector<int>>> sourceSegment, targetSegment;
	ContextPartGraph sourceGraph, targetGraph;
	ContextPartGraphNodeGenerator nodeGen;

	if (!MeshUtil::loadMesh(sourceMeshName, sourceMesh)) return false;
	if (!SegmentGroupApxCvx::loadSegments(sourceSegmentName, sourceSegment)) return false;
	if (!sourceGraph.loadGraphHierarchy(sourceGraphHName, nodeGen, &sourceMesh, &sourceSegment)) return false;
	if (!sourceGraph.loadGraphDescriptor(sourceGraphDName)) return false;
	if (!sourceGraph.loadGraphContext(sourceGraphCName)) return false;

	if (!MeshUtil::loadMesh(targetMeshName, targetMesh)) return false;
	if (!SegmentGroupApxCvx::loadSegments(targetSegmentName, targetSegment)) return false;
	if (!targetGraph.loadGraphHierarchy(targetGraphHName, nodeGen, &targetMesh, &targetSegment)) return false;
	if (!targetGraph.loadGraphDescriptor(targetGraphDName)) return false;
	if (!targetGraph.loadGraphContext(targetGraphCName)) return false;

	// matching...

	return validatePair(&sourceGraph, &targetGraph, weightsFolder, nodePairs, numCorrectPairs);
}

bool ContextPartGraphCrossValidation::validatePair(
	ContextPartGraph *sourceGraph,
	ContextPartGraph *targetGraph,
	string weightsFolder,
	vector<vec2i> &nodePairs,
	int &numCorrectPairs)
{
	Eigen::MatrixXd fullSimilarityMatrix;
	if (!computeSimilarity(sourceGraph, targetGraph, weightsFolder, fullSimilarityMatrix)) return false;

	// only consider first level nodes

	int numSourceValidNodes = (int)sourceGraph->mRootNode->mChildren.size();
	int numTargetValidNodes = (int)targetGraph->mRootNode->mChildren.size();
	Eigen::MatrixXd similarityMatrix = fullSimilarityMatrix.topLeftCorner(numSourceValidNodes, numTargetValidNodes);

	// validation..

	int numNodePairs = (int)nodePairs.size();
	numCorrectPairs = 0;
	for (int pairID = 0; pairID < numNodePairs; pairID++) {
		int sourceID = nodePairs[pairID][0];
		int targetID = nodePairs[pairID][1];

		// all symmetric nodes are valid
		set<int> sourceSet, targetSet;
		sourceSet.insert(sourceID);
		targetSet.insert(targetID);
		for (auto &node : sourceGraph->mAllNodes[sourceID]->mSymmetry) {
			sourceSet.insert(node->mUID);
		}
		for (auto &node : targetGraph->mAllNodes[targetID]->mSymmetry) {
			targetSet.insert(node->mUID);
		}

		int row, col;
		similarityMatrix.row(sourceID).maxCoeff(&row, &col);
		int matchTargetID = col;
		similarityMatrix.col(targetID).maxCoeff(&row, &col);
		int matchSourceID = row;
		if (targetSet.find(matchTargetID) != targetSet.end() ||
			sourceSet.find(matchSourceID) != sourceSet.end())
		{
			numCorrectPairs++;
		}
	}

	return true;
}

bool ContextPartGraphCrossValidation::computeSimilarity(
	ContextPartGraph *sourceGraph,
	ContextPartGraph *targetGraph,
	string weightsFolder,
	Eigen::MatrixXd &similarityMatrix)
{
	Eigen::MatrixXd sourceSelfSimilarity;
	if (!StyleSynthesisConfig::mContext_UseLagaDescriptors) {
		ContextPartGraphMatch cpgm;
		if (!cpgm.loadWeights(weightsFolder)) return false;
		if (!cpgm.loadGraph(*sourceGraph, *sourceGraph)) return false;
		if (!cpgm.process()) return false;
		if (!cpgm.exportSimilarityMatrix(sourceSelfSimilarity)) return false;
	}

	Eigen::MatrixXd targetSelfSimilarity;
	if (!StyleSynthesisConfig::mContext_UseLagaDescriptors) {
		ContextPartGraphMatch cpgm;
		if (!cpgm.loadWeights(weightsFolder)) return false;
		if (!cpgm.loadGraph(*targetGraph, *targetGraph)) return false;
		if (!cpgm.process()) return false;
		if (!cpgm.exportSimilarityMatrix(targetSelfSimilarity)) return false;
	}

	if (true) {
		ContextPartGraphMatch cpgm;

		if (!cpgm.loadWeights(weightsFolder)) return false;
		if (!cpgm.loadGraph(*sourceGraph, *targetGraph)) return false;
		if (!cpgm.process()) return false;
		if (!StyleSynthesisConfig::mContext_UseLagaDescriptors) {
			if (!cpgm.normalizeGraphSimilarity(sourceSelfSimilarity, targetSelfSimilarity)) return false;
		}
		if (!cpgm.exportSimilarityMatrix(similarityMatrix)) return false;

		//// DEBUG: visualization
		//if (!cpgm.matchGraphNodes()) return false;
		//if (!cpgm.visualize("match.ply")) return false;
		//if (!DataUtil::saveMatrixASCII("match.txt", similarityMatrix)) return false;
		//system("pause");
	}

	return true;
}

bool ContextPartGraphCrossValidation::getDifferentMatching(
	ContextPartGraph *sourceGraph,
	ContextPartGraph *targetGraph,
	Eigen::MatrixXd &simMatOurs,
	Eigen::MatrixXd &simMatLaga,
	vector<int> &validSourceNodes,
	vector<vec3i> &validTriplets,
	vector<vec3i> &diffTriplets,
	vec3i &stats)
{

	// only consider first level nodes
	int numSourceNodes = (int)sourceGraph->mRootNode->mChildren.size();
	int numTargetNodes = (int)targetGraph->mRootNode->mChildren.size();
	Eigen::MatrixXd matOurs = simMatOurs.topLeftCorner(numSourceNodes, numTargetNodes);
	Eigen::MatrixXd matLaga = simMatLaga.topLeftCorner(numSourceNodes, numTargetNodes);
	double maxScore = max(matOurs.maxCoeff(), matLaga.maxCoeff());

	diffTriplets.clear();
	validTriplets.clear();
	stats = vec3i(0, 0, 0);
	for (int sourceID : validSourceNodes) {
		int targetIDOurs, targetIDLaga;
		double scoreOurs, scoreLaga;
		if (true) {
			int row, col;
			scoreOurs = matOurs.row(sourceID).maxCoeff(&row, &col);
			targetIDOurs = col;
			for (auto node : sourceGraph->mAllNodes[sourceID]->mSymmetry) {
				double symScore = matOurs.row(node->mUID).maxCoeff(&row, &col);
				if (symScore > scoreOurs) {
					scoreOurs = symScore;
					targetIDOurs = col;
				}
			}
		}
		if (true) {
			int row, col;
			scoreLaga = matLaga.row(sourceID).maxCoeff(&row, &col);
			targetIDLaga = col;
		}
		if (scoreOurs < StyleSynthesisConfig::mData_CustomNumber2) continue; // UNDONE: param matched node pair score threshold

		validTriplets.push_back(vec3i(sourceID, targetIDOurs, targetIDLaga));
		stats[2]++;
		if (targetIDOurs != targetIDLaga) stats[1]++;
		for (auto node : targetGraph->mAllNodes[targetIDOurs]->mSymmetry) {
			if (node->mUID == targetIDLaga) {
				targetIDLaga = targetIDOurs;
				break;
			}
		}
		if (targetIDOurs == targetIDLaga) continue;
		diffTriplets.push_back(vec3i(sourceID, targetIDOurs, targetIDLaga));
		stats[0]++;
	}

	return true;
}

bool ContextPartGraphCrossValidation::loadMatchingData(
	string sourceName, string targetName,
	ContextPartGraphNodeGenerator &nodeGen,
	TTriangleMesh &sourceMesh,
	TTriangleMesh &targetMesh,
	vector<vector<vector<int>>> &sourceSegment,
	vector<vector<vector<int>>> &targetSegment,
	ContextPartGraph &sourceGraph,
	ContextPartGraph &targetGraph)
{
	// names...

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string sourceMeshName = datasetPrefix + "segment/" + sourceName + ".ply";
	string sourceSegmentName = datasetPrefix + "segment/" + sourceName + "-segment.txt";
	string sourceGraphHName = datasetPrefix + "graph/" + sourceName + "/graph-hierarchy.txt";
	string sourceGraphDName = datasetPrefix + "graph/" + sourceName + "/graph-descriptor.txt";
	string sourceGraphCName = datasetPrefix + "graph/" + sourceName + "/graph-context.txt";

	string targetMeshName = datasetPrefix + "segment/" + targetName + ".ply";
	string targetSegmentName = datasetPrefix + "segment/" + targetName + "-segment.txt";
	string targetGraphHName = datasetPrefix + "graph/" + targetName + "/graph-hierarchy.txt";
	string targetGraphDName = datasetPrefix + "graph/" + targetName + "/graph-descriptor.txt";
	string targetGraphCName = datasetPrefix + "graph/" + targetName + "/graph-context.txt";

	// data...

	if (!MeshUtil::loadMesh(sourceMeshName, sourceMesh)) return false;
	if (!SegmentGroupApxCvx::loadSegments(sourceSegmentName, sourceSegment)) return false;
	if (!sourceGraph.loadGraphHierarchy(sourceGraphHName, nodeGen, &sourceMesh, &sourceSegment)) return false;
	if (!sourceGraph.loadGraphDescriptor(sourceGraphDName)) return false;
	if (!sourceGraph.loadGraphContext(sourceGraphCName)) return false;

	if (!MeshUtil::loadMesh(targetMeshName, targetMesh)) return false;
	if (!SegmentGroupApxCvx::loadSegments(targetSegmentName, targetSegment)) return false;
	if (!targetGraph.loadGraphHierarchy(targetGraphHName, nodeGen, &targetMesh, &targetSegment)) return false;
	if (!targetGraph.loadGraphDescriptor(targetGraphDName)) return false;
	if (!targetGraph.loadGraphContext(targetGraphCName)) return false;

	return true;
}

bool ContextPartGraphCrossValidation::loadMatchingPair(
	string fileName,
	vector<string> &sourceList,
	vector<string> &targetList,
	vector<vector<vec2i>> &pairList)
{

	ifstream file(fileName);
	while (!file.eof()) {
		string line;
		getline(file, line);
		line = StringUtil::trim(line);
		if (file.fail() || line.empty()) break;
		vector<string> tokens;
		StringUtil::split(line, ' ', tokens);
		string sourceName = tokens[0];
		string targetName = tokens[1];
		int numPairs = ((int)tokens.size() - 2) / 2;
		vector<vec2i> currentPairList(numPairs);
		for (int k = 0; k < numPairs; k++) {
			currentPairList[k] = vec2i(stoi(tokens[k * 2 + 2]), stoi(tokens[k * 2 + 3]));
		}

		sourceList.push_back(sourceName);
		targetList.push_back(targetName);
		pairList.push_back(currentPairList);
	}
	file.close();

	return true;
}