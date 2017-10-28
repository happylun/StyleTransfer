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

#include "PipelineMatchPartIO.h"

#include <iostream>
#include <fstream>
#include <set>

#include "Mesh/MeshUtil.h"

#include "Segment/SegmentGroupApxCvx.h"
#include "Segment/SegmentUtil.h"

#include "Similarity/SimilarityMetric.h"
#include "Similarity/SimilarityDistance.h"

#include "Context/ContextPartGraph.h"
#include "Context/ContextPartGraphNodeGenerator.h"
#include "Context/ContextPartGraphMatch.h"
#include "Context/ContextPartGraphAssemble.h"
#include "Context/ContextPartGraphTabuSearchPart.h"

#include "Data/StyleSynthesisConfig.h"
#include "Data/DataUtil.h"

using namespace std;
using namespace StyleSynthesis;

bool PipelineMatchPartIO::process() {

	if (false) {
		string sourceName = StyleSynthesisConfig::mData_CustomString1;
		string targetName = StyleSynthesisConfig::mData_CustomString2;
		string sourceShortName = sourceName.substr(sourceName.find_last_of("/\\") + 1);
		string targetShortName = targetName.substr(targetName.find_last_of("/\\") + 1);
		string pairName = sourceShortName + "--" + targetShortName;

		if (!runModelPairs(sourceName, targetName, pairName)) return false;
		if (!organizeResults(sourceName, targetName, pairName)) return false;
		return true;
	}

	if (true) {

		bool useLaga = StyleSynthesisConfig::mContext_UseLagaDescriptors;
		double epsilonScaleFactor = 1.0; // threshold scaling
		string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

		vector<string> sourceList(0), targetList(0), pairList(0);

		string pairListName = datasetPrefix + "assemble/pair-list.txt";
		ifstream pairListFile(pairListName);
		if (!pairListFile.is_open()) return error("cannot open pair list file");
		while (!pairListFile.eof()) {
			string line;
			getline(pairListFile, line);
			line = StringUtil::trim(line);
			if (pairListFile.fail() || line.length() == 0) break;
			stringstream ss(line);
			string sourceName, targetName, pairName;
			ss >> sourceName >> targetName >> pairName;
			if (ss.fail()) return error("incorrect pair list token '" + line + "'");
			sourceList.push_back(sourceName);
			targetList.push_back(targetName);
			pairList.push_back(pairName);
		}
		pairListFile.close();

		string baseConfigName = "base.cfg";
		if (!StyleSynthesisConfig::saveConfig(baseConfigName)) return false;

		int numPairs = (int)pairList.size();
		for (int pairID = 0; pairID < numPairs; pairID++) {
			cout << "=========== Processing pair " << (pairID + 1) << " / " << numPairs << " : " << pairList[pairID] << " ===========" << endl;

			
			if (!StyleSynthesisConfig::loadConfig(baseConfigName)) return false;

			cout << "Loading parameters..." << endl;
			string caseConfigName = datasetPrefix + "assemble/" + pairList[pairID] + "/params.cfg";
			if (!FileUtil::existsfile(caseConfigName)) {
				if (!FileUtil::makedir(caseConfigName)) return false;
				string goodConfigName = datasetPrefix + "assemble/good-params.cfg";
				if (!FileUtil::existsfile(goodConfigName)) {
					goodConfigName = baseConfigName;
				}
				if (!FileUtil::copyfile(goodConfigName, caseConfigName)) return false;
			}

			if (FileUtil::existsfile(caseConfigName)) {
				if (!StyleSynthesisConfig::loadConfig(caseConfigName)) return false;
				if (useLaga) {
					StyleSynthesisConfig::mAssemble_DebugNoBranching = true;
					StyleSynthesisConfig::mContext_UseLagaDescriptors = true;
				}
				StyleSynthesisConfig::mContext_MatchNodeSimilarityThreshold *= epsilonScaleFactor;
			} else return error("missing param data file");

			if (!runModelPairs(sourceList[pairID], targetList[pairID], pairList[pairID])) return false;
			if (!organizeResults(sourceList[pairID], targetList[pairID], pairList[pairID])) return false;
		}
	}

	////////////////////////////////////////////////////////////////////////////////

	if (false) {
		string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

		string matchPairFileName = datasetPrefix + "train/match-shape-pairs.txt";
		if (!generatePairs(matchPairFileName)) return false;

		vector<pair<string, string>> matchPairs;
		if (true) {
			ifstream matchPairFile(matchPairFileName);
			while (!matchPairFile.eof()) {
				string line;
				getline(matchPairFile, line);
				if (matchPairFile.fail()) break;
				line = StringUtil::trim(line);
				if (line.empty()) break;
				vector<string> tokens;
				StringUtil::split(line, ' ', tokens);
				if ((int)tokens.size() != 2) return error("incorrect match pair file");
				matchPairs.push_back(make_pair(tokens[0], tokens[1]));
			}
			matchPairFile.close();
		}

		int numMatchPairs = (int)matchPairs.size();

		for (int pairID = 0; pairID < numMatchPairs; pairID++) {
			string sourceName = matchPairs[pairID].first;
			string targetName = matchPairs[pairID].second;
			string sourceShortName = sourceName.substr(sourceName.find_last_of("/\\") + 1);
			string targetShortName = targetName.substr(targetName.find_last_of("/\\") + 1);
			string pairName = sourceShortName + "--" + targetShortName;

			cout << "============== Matching pairs " << pairName << " ==============" << endl;
			if (!runModelPairs(sourceName, targetName, pairName)) return false;
			if (!organizeResults(sourceName, targetName, pairName)) return false;
		}
	}

	return true;
}

bool PipelineMatchPartIO::generatePairs(string pairFileName) {

	if (FileUtil::existsfile(pairFileName)) return true; // already generated

	// names

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;
	string meshListFileName = datasetPrefix + "mesh/mesh-list-all.txt";
	string meshLabelFileName = datasetPrefix + "mesh/mesh-labels-all.txt";
	string similarityFileName = datasetPrefix + "train/all-pair-similarity-final.txt";
	string exemplarFileName = datasetPrefix + "train/exemplar-list.txt";

	// data

	vector<string> meshNameList(0);
	vector<string> meshSceneList(0);
	vector<int> meshLabelList(0);
	Eigen::MatrixXd allPairSimilarity;
	map<string, int> meshNameMap;
	vector<int> exemplarList(0);

	ifstream meshListFile(meshListFileName);
	while (!meshListFile.eof()) {
		string line;
		getline(meshListFile, line);
		if (line.empty()) break;

		string meshName = StringUtil::trim(line);
		vector<string> tokens;
		StringUtil::split(meshName, '/', tokens);
		string sceneName = tokens[0];

		meshNameMap[meshName] = (int)meshNameList.size();
		meshNameList.push_back(meshName);
		meshSceneList.push_back(sceneName);
	}
	meshListFile.close();

	int numMeshes = (int)meshNameList.size();

	ifstream meshLabelFile(meshLabelFileName);
	for (int meshID = 0; meshID < numMeshes; meshID++) {
		int label;
		meshLabelFile >> label;
		meshLabelList.push_back(label);
	}
	meshLabelFile.close();

	if (!DataUtil::loadMatrixASCII(similarityFileName, allPairSimilarity)) return false;

	ifstream exemplarListFile(exemplarFileName);
	while (!exemplarListFile.eof()) {
		string line;
		getline(exemplarListFile, line);
		if (line.empty()) break;

		string meshName = StringUtil::trim(line);
		int meshID = meshNameMap[meshName];
		exemplarList.push_back(meshID);
	}
	exemplarListFile.close();

	// generate match pairs

	vector<pair<string, string>> matchPairs(0);
	vector<double> matchScores(0);

	//for (int meshID = 0; meshID < numMeshes; meshID++) {
	//	int matchID = meshID;
	//	double matchScore = -1;
	//	for (int otherID = 0; otherID < numMeshes; otherID++) {
	//		if (meshSceneList[meshID] == meshSceneList[otherID]) continue;
	//		if (meshLabelList[meshID] == meshLabelList[otherID]) continue;
	//		double score = allPairSimilarity(meshID, otherID);
	//		if (score > matchScore) {
	//			matchScore = score;
	//			matchID = otherID;
	//		}
	//	}
	//	matchPairs.push_back(make_pair(meshNameList[meshID], meshNameList[matchID]));
	//	matchScores.push_back(matchScore);
	//}

	if (true) {
		int numLabels = 6;
		for (int meshID : exemplarList) {
			int label = meshLabelList[meshID];
			vector<double> maxScore(numLabels, -1);
			vector<int> maxID(numLabels, -1);
			for (int otherID = 0; otherID < numMeshes; otherID++) {
				if (meshSceneList[meshID] == meshSceneList[otherID]) continue;
				int otherLabel = meshLabelList[otherID];
				if (otherLabel == label) continue;
				double score = allPairSimilarity(meshID, otherID);
				if (score > maxScore[otherLabel]) {
					maxScore[otherLabel] = score;
					maxID[otherLabel] = otherID;
				}
			}
			for (int otherLabel = 0; otherLabel < numLabels; otherLabel++) {
				int otherID = maxID[otherLabel];
				if (otherID < 0) continue;
				matchPairs.push_back(make_pair(meshNameList[meshID], meshNameList[otherID]));
				matchScores.push_back(maxScore[otherLabel]);
			}
		}
	}


	// sort match pairs

	int numMatches = (int)matchPairs.size();
	vector<int> matchOrder(numMatches);
	for (int k = 0; k < numMatches; k++) matchOrder[k] = k;
	//sort(matchOrder.begin(), matchOrder.end(),
	//	[&matchScores](int lhs, int rhs) {return matchScores[lhs] > matchScores[rhs]; });

	// output pairs

	ofstream pairFile(pairFileName);
	for (int matchID = 0; matchID < numMatches; matchID++) {
		pairFile << matchPairs[matchOrder[matchID]].first << " ";
		pairFile << matchPairs[matchOrder[matchID]].second << endl;
	}
	pairFile.close();

	return true;
}

bool PipelineMatchPartIO::runModelPairs(string sourceName, string targetName, string pairName) {

	// names

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string sourceMeshName = datasetPrefix + "segment/" + sourceName + ".ply";
	string targetMeshName = datasetPrefix + "segment/" + targetName + ".ply";

	string sourceSegmentName = datasetPrefix + "segment/" + sourceName + "-segment.txt";
	string targetSegmentName = datasetPrefix + "segment/" + targetName + "-segment.txt";

	string sourceGraphFolder = datasetPrefix + "graph/" + sourceName + "/";
	string targetGraphFolder = datasetPrefix + "graph/" + targetName + "/";

	string sourceGraphHName = sourceGraphFolder + "graph-hierarchy.txt";
	string sourceGraphDName = sourceGraphFolder + "graph-descriptor.txt";
	string sourceGraphCName = sourceGraphFolder + "graph-context.txt";

	string targetGraphHName = targetGraphFolder + "graph-hierarchy.txt";
	string targetGraphDName = targetGraphFolder + "graph-descriptor.txt";
	string targetGraphCName = targetGraphFolder + "graph-context.txt";

	string weightsFolder = datasetPrefix + "weights/";

	string assembleFolder = datasetPrefix + "assemble/" + pairName + "/";
	if (!FileUtil::makedir(assembleFolder)) return false;

	string matchDataName = assembleFolder + "match.txt";
	string matchSimilarityName = assembleFolder + "similarity.txt";
	string matchGroupSimilarityName = assembleFolder + "group-similarity.txt";
	string matchVisualName = assembleFolder + "match.ply";

	string shapeStyleDistanceName = assembleFolder + "shape-style-distance.txt";
	string sourceNodeContributionName = assembleFolder + "source-contribution.txt";
	string targetNodeContributionName = assembleFolder + "target-contribution.txt";
	string sourceNodeSaliencyName = assembleFolder + "source-saliency.txt";
	string targetNodeSaliencyName = assembleFolder + "target-saliency.txt";
	string matchContributionName = assembleFolder + "match-contribution.txt";
	string matchOrderName = assembleFolder + "match-order.txt";
	string sourceKeyGroupName = assembleFolder + "source-key-group.txt";
	string targetKeyGroupName = assembleFolder + "target-key-group.txt";

	string assembleResultName = assembleFolder + "result.txt";

	if (FileUtil::existsfile(assembleResultName)) return true; // early quit
	if (FileUtil::existsfile(assembleFolder + "search-result.txt")) return true; // early quit
	
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

	// algorithm

	cout << "Computing node contribution & saliency..." << endl;

	double shapeStyleDistance;
	vector<double> sourceNodeContribution;
	vector<double> targetNodeContribution;
	vector<double> sourceNodeSaliency;
	vector<double> targetNodeSaliency;
	if (FileUtil::existsfile(sourceNodeContributionName) && FileUtil::existsfile(targetNodeContributionName)) {
		if (!DataUtil::loadValueListASCII(sourceNodeContributionName, sourceNodeContribution)) return false;
		if (!DataUtil::loadValueListASCII(targetNodeContributionName, targetNodeContribution)) return false;
		if (!DataUtil::loadValueListASCII(sourceNodeSaliencyName, sourceNodeSaliency)) return false;
		if (!DataUtil::loadValueListASCII(targetNodeSaliencyName, targetNodeSaliency)) return false;
	} else {
		if (!computeNodeContribution(
			sourceName, targetName,
			&sourceGraph, &targetGraph,
			&sourceSegments, &targetSegments,
			shapeStyleDistance,
			sourceNodeContribution, targetNodeContribution,
			sourceNodeSaliency, targetNodeSaliency)) return false;
		if (!DataUtil::saveValueListASCII(sourceNodeContributionName, sourceNodeContribution)) return false;
		if (!DataUtil::saveValueListASCII(targetNodeContributionName, targetNodeContribution)) return false;
		if (!DataUtil::saveValueListASCII(sourceNodeSaliencyName, sourceNodeSaliency)) return false;
		if (!DataUtil::saveValueListASCII(targetNodeSaliencyName, targetNodeSaliency)) return false;
	}

	cout << "Graph matching..." << endl;

	if (!FileUtil::existsfile(matchDataName)) {

		// compute self similarity
		Eigen::MatrixXd sourceSelfSimilarity;
		if (true) {
			cout << "Computing source self-similarity" << endl;
			ContextPartGraphMatch cpgm;
			if (!cpgm.loadWeights(weightsFolder)) return false;
			if (!cpgm.loadGraph(sourceGraph, sourceGraph)) return false;
			if (!cpgm.process()) return false;
			if (!cpgm.exportSimilarityMatrix(sourceSelfSimilarity)) return false;
		}
		Eigen::MatrixXd targetSelfSimilarity;
		if (true) {
			cout << "Computing target self-similarity" << endl;
			ContextPartGraphMatch cpgm;
			if (!cpgm.loadWeights(weightsFolder)) return false;
			if (!cpgm.loadGraph(targetGraph, targetGraph)) return false;
			if (!cpgm.process()) return false;
			if (!cpgm.exportSimilarityMatrix(targetSelfSimilarity)) return false;
		}

		// do matching
		vector<vector<vec2i>> matchings;
		vector<double> matchingsSimilarity;
		ContextPartGraphMatch cpgm;
		if (!cpgm.loadWeights(weightsFolder)) return false;
		if (!cpgm.loadGraph(sourceGraph, targetGraph)) return false;
		if (!cpgm.process()) return false;
		if (!cpgm.normalizeGraphSimilarity(sourceSelfSimilarity, targetSelfSimilarity)) return false;
		if (!cpgm.matchGraphNodes()) return false;
		if (!cpgm.exportMatchings(matchings)) return false;
		if (!cpgm.exportMatchingsSimilarity(matchingsSimilarity)) return false;
		if (!cpgm.visualize(matchVisualName)) return false;

		// save similarity matrix
		Eigen::MatrixXd similarityMatrix;
		if (!cpgm.exportSimilarityMatrix(similarityMatrix)) return false;
		if (!DataUtil::saveMatrixASCII(matchSimilarityName, similarityMatrix)) return false;
		if (!DataUtil::saveValueListASCII(matchGroupSimilarityName, matchingsSimilarity)) return false;

		// compute contribution
		int numMatchings = (int)matchings.size();
		vector<double> matchingsScore(numMatchings, 0);
		for (int groupID = 0; groupID < numMatchings; groupID++) {
			auto &matchingGroup = matchings[groupID];
			for (vec2i &matchingPair : matchingGroup) {
				matchingsScore[groupID] += sourceNodeContribution[matchingPair[0]]
					+ targetNodeContribution[matchingPair[1]];
			}
		}
		if (!DataUtil::saveValueListASCII(matchContributionName, matchingsScore)) return false;

		// sort by contribution
		vector<int> matchingsOrder(numMatchings);
		if (true) {
			for (int k = 0; k < numMatchings; k++) matchingsOrder[k] = k;
			//sort(matchingsOrder.begin(), matchingsOrder.end(),
			//	[&matchingsScore](int lhs, int rhs){return matchingsScore[lhs] > matchingsScore[rhs]; });
		}
		if (!DataUtil::saveIndexListASCII(matchOrderName, matchingsOrder)) return false;

		// skip non-replaceable matchings
		if (true) {
			double threshold = StyleSynthesisConfig::mStyle_ReplaceableContributionThreshold;
			vector<vector<vec2i>> newMatchings(0);
			vector<double> newMatchingsScore(0);
			for (int k = 0; k < numMatchings; k++) {
				double score = matchingsScore[matchingsOrder[k]];
				if (score < threshold) continue;
				newMatchings.push_back(matchings[matchingsOrder[k]]);
				newMatchingsScore.push_back(score);
			}
			matchings.swap(newMatchings);
			matchingsScore.swap(newMatchingsScore);
		}

		if (!ContextPartGraphMatch::saveMatchings(matchDataName, matchings)) return false;		
	}

	if (!FileUtil::existsfile(sourceKeyGroupName) || !FileUtil::existsfile(targetKeyGroupName)) {

		vector<vector<int>> srcGroups, tgtGroups;
		vector<double> srcGroupScore, tgtGroupScore;
		if (!getKeyGroups(&sourceGraph, sourceNodeContribution, srcGroups, srcGroupScore)) return false;
		if (!getKeyGroups(&targetGraph, targetNodeContribution, tgtGroups, tgtGroupScore)) return false;

		// prune groups by contributions
		double threshold = StyleSynthesisConfig::mStyle_TransferrableContributionThreshold;
		vector<vector<int>> newSrcGroups, newTgtGroups;
		for (int groupID = 0; groupID < (int)srcGroups.size(); groupID++) {
			cout << "Source key group: " << srcGroupScore[groupID] << ": ";
			for (int id : srcGroups[groupID]) cout << id << " ";
			cout << endl;
			if (srcGroupScore[groupID] >= threshold) {
				newSrcGroups.push_back(srcGroups[groupID]);
			}
		}
		for (int groupID = 0; groupID < (int)tgtGroups.size(); groupID++) {
			cout << "Target key group: " << tgtGroupScore[groupID] << ": ";
			for (int id : tgtGroups[groupID]) cout << id << " ";
			cout << endl;
			if (tgtGroupScore[groupID] >= threshold) {
				newTgtGroups.push_back(tgtGroups[groupID]);
			}
		}

		if (!DataUtil::saveGroupListASCII(sourceKeyGroupName, newSrcGroups)) return false;
		if (!DataUtil::saveGroupListASCII(targetKeyGroupName, newTgtGroups)) return false;
	}

	//return true;
	//system("pause");

	vector<vector<vec2i>> matchings;
	if (!ContextPartGraphMatch::loadMatchings(matchDataName, matchings)) return false;

	vector<vector<int>> sourceKeyGroups, targetKeyGroups;
	if (!DataUtil::loadGroupListASCII(sourceKeyGroupName, sourceKeyGroups)) return false;
	if (!DataUtil::loadGroupListASCII(targetKeyGroupName, targetKeyGroups)) return false;

	cout << "Running tabu search for part replacement..." << endl;

	if (true) {

		ContextPartGraphTabuSearchPart cpgtsp;
		if (!cpgtsp.loadGraph(&cpgng, &sourceGraph, &targetGraph)) return false;
		if (!cpgtsp.loadMatchings(matchings)) return false;
		if (!cpgtsp.loadContributions(
			sourceNodeContribution, targetNodeContribution,
			sourceNodeSaliency, targetNodeSaliency)) return false;
		if (!cpgtsp.loadKeyGroups(sourceKeyGroups, targetKeyGroups)) return false;
		if (!cpgtsp.loadNames(weightsFolder, assembleFolder)) return false;
		if (!cpgtsp.process()) return false;

		int numSolutions = cpgtsp.getNumSolutions();
		if (!DataUtil::saveIndexListASCII(assembleResultName, vector<int>(1, numSolutions))) return false;
	}

	return true;
}


bool PipelineMatchPartIO::organizeResults(string sourceName, string targetName, string pairName) {

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;
	string assembleFolder = datasetPrefix + "assemble/" + pairName + "/";

	string searchResultName = assembleFolder + "result.txt";
	if (!FileUtil::existsfile(searchResultName)) return error("run batch search first");

	string allOutFolder = datasetPrefix + "output/" + pairName + "/";
	if (!FileUtil::makedir(allOutFolder)) return false;

	string outFolder = assembleFolder + "solution-search/";
	if (!FileUtil::makedir(outFolder)) return false;

	if (true) {
		string inExemplarName = datasetPrefix + "segment/" + sourceName + ".ply";
		string inCandidateName = datasetPrefix + "segment/" + targetName + ".ply";

		string allOutExemplarName = allOutFolder + "exemplar.ply";
		string allOutCandidateName = allOutFolder + "candidate.ply";
		//if (FileUtil::existsfile(allOutExemplarName)) return true; // early quit
		if (!FileUtil::copyfile(inExemplarName, allOutExemplarName)) return false;
		if (!FileUtil::copyfile(inCandidateName, allOutCandidateName)) return false;

		string outExemplarName = outFolder + "exemplar.ply";
		string outCandidateName = outFolder + "candidate.ply";
		if (!FileUtil::copyfile(inExemplarName, outExemplarName)) return false;
		if (!FileUtil::copyfile(inCandidateName, outCandidateName)) return false;
	}

	if (true) {
		string solutionListName = assembleFolder + "solution-sorted-list.txt";
		ifstream solutionListFile(solutionListName);
		if (!solutionListFile.is_open()) return error("cannot open solution list file");
		int numSolutions = 0;
		while (!solutionListFile.eof()) {
			int solutionID;
			double solutionDistance;
			solutionListFile >> solutionID >> solutionDistance;
			if (solutionID == 0 || solutionListFile.fail()) break;

			string inResultName = assembleFolder + "solution-" + to_string(solutionID) + "/mesh.ply";

			string allOutResultName = allOutFolder + "mesh-" + to_string(numSolutions) + ".ply";
			if (!FileUtil::copyfile(inResultName, allOutResultName)) return false;

			string outResultName = outFolder + "mesh-" + to_string(numSolutions) + ".ply";
			if (!FileUtil::copyfile(inResultName, outResultName)) return false;

			numSolutions++;
		}
		solutionListFile.close();
	}

	return true;
}

bool PipelineMatchPartIO::computeNodeContribution(
	string sourceName, string targetName,
	ContextPartGraph *sourceGraph,
	ContextPartGraph *targetGraph,
	vector<vector<vector<int>>> *sourceSegments,
	vector<vector<vector<int>>> *targetSegments,
	double &shapeStyleDistance,
	vector<double> &sourceNodeContribution,
	vector<double> &targetNodeContribution,
	vector<double> &sourceNodeSaliency,
	vector<double> &targetNodeSaliency)
{

	// names

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string weightsFolder = datasetPrefix + "weights/";

	string sourceShortName = sourceName.substr(sourceName.find_last_of("/\\") + 1);
	string targetShortName = targetName.substr(targetName.find_last_of("/\\") + 1);

	string elementFolder = datasetPrefix + "element/" + sourceShortName + "--" + targetShortName + "/";
	if (!FileUtil::makedir(elementFolder)) return false;

	string elementVisualName = elementFolder + "vis-element.ply";

	// do style matching

	cout << "Doing style matching" << endl;

	if (!SimilarityMetric::initMetric(weightsFolder)) return false;

	SimilarityData data;
	if (!data.loadData(sourceName, targetName)) return false;

	SimilarityDistance sd(&data);
	if (!sd.process()) return false;
	if (!sd.output(elementFolder)) return false;
	if (!sd.visualize(elementVisualName)) return false;

	// compute node contribution

	cout << "Computing node contribution" << endl;

	vector<double> sourceFaceContribution;
	vector<double> targetFaceContribution;
	vector<double> sourceFaceSaliency;
	vector<double> targetFaceSaliency;
	if (!sd.evaluateFaceContribution(
		shapeStyleDistance,
		sourceFaceContribution,
		targetFaceContribution,
		sourceFaceSaliency,
		targetFaceSaliency)) return false;

	int numSourceNodes = (int)sourceGraph->mAllNodes.size();
	int numTargetNodes = (int)targetGraph->mAllNodes.size();
	sourceNodeContribution.resize(numSourceNodes);
	targetNodeContribution.resize(numTargetNodes);
	sourceNodeSaliency.resize(numSourceNodes);
	targetNodeSaliency.resize(numTargetNodes);

	for (int nodeID = 0; nodeID < numSourceNodes; nodeID++) {
		ContextPartGraphNode *node = sourceGraph->mAllNodes[nodeID];
		auto &patch = (*sourceSegments)[node->mPartLevelID][node->mPartSegmentID];

		double contribution = 0;
		for (int faceID : patch) contribution += sourceFaceContribution[faceID];
		sourceNodeContribution[nodeID] = contribution;

		double saliency = 0;
		for (int faceID : patch) saliency += sourceFaceSaliency[faceID];
		sourceNodeSaliency[nodeID] = saliency;

		if (nodeID < (int)sourceGraph->mRootNode->mChildren.size()) {
			cout << "Source node " << nodeID << " contribution = " << contribution << "; saliency = " << saliency << endl;
		}
	}
	cout << "==========================" << endl;

	for (int nodeID = 0; nodeID < numTargetNodes; nodeID++) {
		ContextPartGraphNode *node = targetGraph->mAllNodes[nodeID];
		auto &patch = (*targetSegments)[node->mPartLevelID][node->mPartSegmentID];

		double contribution = 0;
		for (int faceID : patch) contribution += targetFaceContribution[faceID];
		targetNodeContribution[nodeID] = contribution;

		double saliency = 0;
		for (int faceID : patch) saliency += targetFaceSaliency[faceID];
		targetNodeSaliency[nodeID] = saliency;

		if (nodeID < (int)targetGraph->mRootNode->mChildren.size()) {
			cout << "Target node " << nodeID << " contribution = " << contribution << "; saliency = " << saliency << endl;
		}
	}
	cout << "==========================" << endl;

	//system("pause");

	return true;
}

bool PipelineMatchPartIO::getKeyGroups(
	ContextPartGraph *graph,
	vector<double> &nodeContributions,
	vector<vector<int>> &outGroups,
	vector<double> &outGroupContributions)
{

	int numNodes = (int)graph->mAllNodes.size();

	vector<bool> validFlags(numNodes, false);
	for (int nodeID = 0; nodeID < numNodes; nodeID++) {
		if (graph->mAllNodes[nodeID]->mParent->mUID < 0) { // HACK: only use nodes in first level
			validFlags[nodeID] = true;
		}
	}

	outGroups.clear();
	outGroupContributions.clear();

	vector<bool> visitedFlags(numNodes, false);
	for (int nodeID = 0; nodeID < numNodes; nodeID++) {
		if (!validFlags[nodeID]) continue;
		if (visitedFlags[nodeID]) continue;
		visitedFlags[nodeID] = true;

		vector<int> clique(1, nodeID);
		for (auto nbNode : graph->mAllNodes[nodeID]->mSymmetry) {
			int nbNodeID = nbNode->mUID;
			if (!validFlags[nbNodeID]) continue;
			if (visitedFlags[nbNodeID]) continue; // it won't happen
			visitedFlags[nbNodeID] = true;
			clique.push_back(nbNodeID);
		}

		double contribution = 0;
		for (int id : clique) {
			contribution += nodeContributions[id];
		}

		outGroups.push_back(clique);
		outGroupContributions.push_back(contribution);
	}

	return true;
}

bool PipelineMatchPartIO::error(string s) {
	cout << "Error: " << s << endl;
	return false;
}