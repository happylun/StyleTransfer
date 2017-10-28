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

#include "ContextPartGraphTabuSearchPart.h"

#include <iostream>
#include <fstream>
#include <set>

#include "Context/ContextPartGraphMatch.h"
#include "Context/ContextPartGraphAssemble.h"
#include "Context/ContextPartGraphAssembleResult.h"

#include "Utility/FileUtil.h"

#include "Data/DataUtil.h"
#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

#define OUTPUT_PROGRESS

ContextPartGraphTabuSearchPart::ContextPartGraphTabuSearchPart() {
}

ContextPartGraphTabuSearchPart::~ContextPartGraphTabuSearchPart() {
	if (!cleanUp()) error("failed to cleanup memory for Tabu search");
}

bool ContextPartGraphTabuSearchPart::loadGraph(TNodeGen *generator, TGraph *source, TGraph *target) {
	mpNodeGenerator = generator;
	mpSourceGraph = source;
	mpTargetGraph = target;
	return true;
}

bool ContextPartGraphTabuSearchPart::loadMatchings(vector<vector<vec2i>> &matchings) {
	mMatchings = matchings;
	return true;
}

bool ContextPartGraphTabuSearchPart::loadContributions(
	vector<double> &sourceContribution,
	vector<double> &targetContribution,
	vector<double> &sourceSaliency,
	vector<double> &targetSaliency)
{
	mSourceContribution = sourceContribution;
	mTargetContribution = targetContribution;
	mSourceSaliency = sourceSaliency;
	mTargetSaliency = targetSaliency;
	return true;
}

bool ContextPartGraphTabuSearchPart::loadKeyGroups(vector<vector<int>> &sourceGroups, vector<vector<int>> &targetGroups) {
	mSourceKeyGroups = sourceGroups;
	mTargetKeyGroups = targetGroups;
	return true;
}

bool ContextPartGraphTabuSearchPart::loadNames(string weightsFolder, string resultFolder) {
	mWeightsFolder = weightsFolder;
	mResultFolder = resultFolder;
	return true;
}

bool ContextPartGraphTabuSearchPart::process() {

	if (!initialize()) return false;
	if (!runTabuSearch()) return false;
	if (!analyzeSolutions()) return false;

	return true;
}

bool ContextPartGraphTabuSearchPart::initialize() {

	// put identity solution (no replacement) into solution pool

	TSolution *identitySolution = new TSolution();
	if (!identitySolution->identity(mpTargetGraph)) return false;

	string solutionFolder = mResultFolder + "solution-0/";
	if (!FileUtil::makedir(solutionFolder)) return false;
	string similarityName = solutionFolder + "similarity.txt";
	string thresholdName = mResultFolder + "threshold.txt";

	if (!identitySolution->saveData(solutionFolder)) return false;
	if (!DataUtil::saveValueListASCII(similarityName, vector<double>(1, 1.0))) return false;
	if (!DataUtil::saveValueListASCII(thresholdName, vector<double>(1, StyleSynthesisConfig::mStyle_ReplaceableContributionThreshold))) return false;

	mSolutionPool.push_back(identitySolution);
	mSolutionSourceContributors.push_back(vector<int>(0));
	mSolutionTargetContributors.push_back(vector<int>(0));

	// compute identity similarity

	if (true) {
		ContextPartGraphMatch cpgm;
		if (!cpgm.loadWeights(mWeightsFolder)) return false;
		if (!cpgm.loadGraph(*mpTargetGraph, *mpTargetGraph)) return false;
		if (!cpgm.process()) return false;
		if (!cpgm.exportSimilarityMatrix(mOriginalIdentitySimilarity)) return false;
	}

	return true;
}

bool ContextPartGraphTabuSearchPart::runTabuSearch() {

	if (mSolutionPool.empty()) return error("solution pool is not initialized");

	int maxNumGroups = 20; // UNDONE: param maximum number of groups to be considered in tabu search

	// test apply / not apply replacing for each matching group one by one on all previous solutions

	int numMatchingGroups = min(maxNumGroups, (int)mMatchings.size());

	for (int matchGroupID = 0; matchGroupID < numMatchingGroups; matchGroupID++) {
#ifdef OUTPUT_PROGRESS
		cout << "============== Testing matching group " << (matchGroupID + 1) << " / " << numMatchingGroups << " ==============" << endl;
#endif

		int numCurrentSolutions = (int)mSolutionPool.size();
		int startSolution = StyleSynthesisConfig::mAssemble_DebugNoBranching ? numCurrentSolutions - 1 : 0;
		for (int solID = startSolution; solID < numCurrentSolutions; solID++) {
#ifdef OUTPUT_PROGRESS
			cout << "============== Testing solution " << (solID + 1) << " / " << numCurrentSolutions << " ==============" << endl;
#endif

			TSolution *currentSolution = mSolutionPool[solID];
			vector<int> &currentSourceContributors = mSolutionSourceContributors[solID];
			vector<int> &currentTargetContributors = mSolutionTargetContributors[solID];

			// map node

			vector<vec2i> currentMatching;
#ifdef OUTPUT_PROGRESS
			cout << "Matching: ";
#endif
			for (vec2i matchPair : mMatchings[matchGroupID]) {
#ifdef OUTPUT_PROGRESS
				cout << matchPair[0] << "-" << matchPair[1] << "  ";
#endif
				vec2i mappedMatchPair = matchPair;
				mappedMatchPair[1] = currentSolution->mNodeMapping[matchPair[1]];
				if (mappedMatchPair[1] < 0 || mappedMatchPair[1] >= (int)currentSolution->mGraph.mAllNodes.size()) {
					return error("incorrect node mapping");
				}
				currentMatching.push_back(mappedMatchPair);
			}
#ifdef OUTPUT_PROGRESS
			cout << endl;
#endif

			// part assembly

			ContextPartGraphAssemble cpga;
			bool sucessFlag;
			if (!cpga.loadGraph(mpNodeGenerator, mpSourceGraph, &currentSolution->mGraph)) return false;
			if (!cpga.loadOldSolution(mSolutionPool[0], currentSolution)) return false;
			if (!cpga.loadMatching(currentMatching)) return false;
			if (!cpga.process(sucessFlag)) return false;

			// generate solution

			int solutionID = (int)mSolutionPool.size();
			string solutionFolder = mResultFolder + "solution-" + to_string(solutionID) + "/";
			if (!FileUtil::makedir(solutionFolder)) return false; // if there is anything inside this folder
			string assembleVName = solutionFolder + "assemble.ply";
			string similarityName = solutionFolder + "similarity.txt";

			TSolution *newSolution = new TSolution();
			if (!newSolution->generate(&cpga)) return false;
			if (!newSolution->chain(currentSolution)) return false;

			// calculate style score

			vector<int> newSourceContributors = currentSourceContributors;
			vector<int> newTargetContributors = currentTargetContributors;
			for (vec2i matchPair : mMatchings[matchGroupID]) {
				newSourceContributors.push_back(matchPair[0]);
				newTargetContributors.push_back(matchPair[1]);
			}
			double newStyleScore;
			if (!computeStyleScore(newSourceContributors, newTargetContributors, newStyleScore)) return false;

			// check functionality compatibility

			double graphSimilarity;
			if (!computeFunctionalSimilarity(newSolution, graphSimilarity)) return false;

#ifdef OUTPUT_PROGRESS
			cout << "Graph score = " << graphSimilarity << endl;
			cout << "Style score = " << newStyleScore << endl;
#endif
			if (graphSimilarity < 0) continue; // UNDONE: param functionality compatibility threshold

			if (!cpga.visualize(solutionFolder)) return false;
			if (!newSolution->visualize(assembleVName, currentSolution)) return false;
			if (!newSolution->saveData(solutionFolder)) return false;
			if (!DataUtil::saveValueListASCII(similarityName, vector<double>(1, graphSimilarity))) return false;

			if (sucessFlag) {
				mSolutionPool.push_back(newSolution);
				mSolutionSourceContributors.push_back(newSourceContributors);
				mSolutionTargetContributors.push_back(newTargetContributors);
			} else {
#ifdef OUTPUT_PROGRESS
				cout << "========= Failed to go any further from this solution... =========" << endl;
#endif
				//system("pause");
			}
			//system("pause");
		}
	}

	if (true) {
		// HACK: only do removing / post-alignment / adding in last solution
		int lastSolutionID = (int)mSolutionPool.size() - 1;
		if (!finalizeSolution(lastSolutionID)) return false;
	}

	if (false) {
		// do removing / post-alignment / adding in each valid solution
		int numValidSolutions = (int)mSolutionPool.size();
		for (int solutionID = 0; solutionID < numValidSolutions; solutionID++) {
			if (!finalizeSolution(solutionID)) return false;
		}
	}

	return true;
}

bool ContextPartGraphTabuSearchPart::analyzeSolutions() {

	int numSolutions = (int)mSolutionPool.size();

	vector<double> scoreList(0);
	for (int solID = 0; solID < numSolutions; solID++) {
		double score;
		if (!computeStyleScore(
			mSolutionSourceContributors[solID],
			mSolutionTargetContributors[solID],
			score)) return false;
		scoreList.push_back(score);
	}

	ofstream listFile(mResultFolder + "solution-list.txt");
	cout << "Solution list:" << endl;
	for (int solID = 0; solID < numSolutions; solID++) {
		listFile << solID << "   " << scoreList[solID] << endl;
		cout << solID << "   " << scoreList[solID] << endl;
	}
	listFile.close();

	vector<int> orderList(numSolutions);
	for (int k = 0; k < numSolutions; k++) orderList[k] = k;
	
	sort(orderList.begin(), orderList.end(),
		[&scoreList](int &lhs, int &rhs) {
			return scoreList[lhs] < scoreList[rhs] || // most stylistically similar ones come first
				(scoreList[lhs] == scoreList[rhs] && lhs > rhs); }); // if equally similar, later solution come first

	ofstream sortedListFile(mResultFolder + "solution-sorted-list.txt");
	cout << "Sorted solution list:" << endl;
	for (int orderID = 0; orderID < numSolutions; orderID++) {
		int solID = orderList[orderID];
		//if (orderID > 0 && scoreList[orderList[orderID - 1]] == scoreList[solID]) continue; // skip solution before post-alignment
		sortedListFile << solID << "   " << scoreList[solID] << endl;
		cout << solID << "   " << scoreList[solID] << endl;
	}
	sortedListFile.close();

	return true;
}

bool ContextPartGraphTabuSearchPart::finalizeSolution(int solutionID) {

#ifdef OUTPUT_PROGRESS
	cout << "================= Finalizing solution =================" << endl;
#endif

	TSolution *currentSolution = mSolutionPool[solutionID];
	vector<int> currentSourceContributors = mSolutionSourceContributors[solutionID];
	vector<int> currentTargetContributors = mSolutionTargetContributors[solutionID];

	// get set of replaced nodes
	set<int> replacedSourceSet, replacedTargetSet;
	for (vec2i replacedPair : currentSolution->mReplaceMapping) {
		replacedSourceSet.insert(replacedPair[0]);
	}
	for (int targetID = 0; targetID < (int)currentSolution->mNodeMapping.size(); targetID++) {
		if (currentSolution->mNodeMapping[targetID] < 0) {
			replacedTargetSet.insert(targetID);
		}
	}

	int maxRemoving = 3; // UNDONE: param maximum number of removing/adding operations
	int maxAdding = 3;
	double styleScoreThreshold = 0.0; // UNDONE: param style threshold

	int numSourceKeyGroups = (int)mSourceKeyGroups.size();
	int numTargetKeyGroups = (int)mTargetKeyGroups.size();
	if (!StyleSynthesisConfig::mAssemble_AllowRemovingParts) numTargetKeyGroups = 0;
	if (!StyleSynthesisConfig::mAssemble_AllowAddingParts || solutionID == 0) numSourceKeyGroups = 0;

	int numRemoving = 0;
	for (int groupID = 0; groupID < numTargetKeyGroups; groupID++) {

		vector<int> &group = mTargetKeyGroups[groupID];
		bool isValid = true;
		for (int nodeID : group) {
			if (replacedTargetSet.find(nodeID) != replacedTargetSet.end()) {
				isValid = false;
				break;
			}
		}
		if (!isValid) continue;

#ifdef OUTPUT_PROGRESS
		cout << "================= Removing group " << (groupID + 1) << " =================" << endl;
#endif

		vector<int> addingList(0);
		vector<int> removingList;

		for (int nodeID : group) {
			int mappedNodeID = currentSolution->mNodeMapping[nodeID];
			if (mappedNodeID >= 0) removingList.push_back(mappedNodeID);
		}

		TSolution *newSolution;
		if (!runAddingRemoving(addingList, removingList, currentSolution, newSolution)) return false;
		if (currentSolution != newSolution) {

			vector<int> newSourceContributors = currentSourceContributors;
			vector<int> newTargetContributors = currentTargetContributors;
			for (int nodeID : group) {
				newTargetContributors.push_back(nodeID);
			}
			double newStyleScore;
			if (!computeStyleScore(newSourceContributors, newTargetContributors, newStyleScore)) return false;
#ifdef OUTPUT_PROGRESS
			cout << "Style score = " << newStyleScore << endl;
#endif
			if (newStyleScore < styleScoreThreshold) {
				if(newSolution) delete newSolution;
				newSolution = currentSolution;
			} else {
				mSolutionPool.push_back(newSolution);
				mSolutionSourceContributors.push_back(newSourceContributors);
				mSolutionTargetContributors.push_back(newTargetContributors);

				currentSolution = newSolution;
				currentSourceContributors = newSourceContributors;
				currentTargetContributors = newTargetContributors;

				numRemoving++;
				if (numRemoving >= maxRemoving) break;
			}
		}
	}

	if (StyleSynthesisConfig::mAssemble_AllowPostAlignment && solutionID > 0) {
#ifdef OUTPUT_PROGRESS
		cout << "================= Re-aligning solution =================" << endl;
#endif

		TSolution *newSolution;
		if (!runPostProcessing(currentSolution, newSolution)) return false;
		if (currentSolution != newSolution) {

			mSolutionPool.push_back(newSolution);
			mSolutionSourceContributors.push_back(currentSourceContributors); // unchanged
			mSolutionTargetContributors.push_back(currentTargetContributors);

			currentSolution = newSolution;
		}
	}

	int numAdding = 0;
	for (int groupID = 0; groupID < numSourceKeyGroups; groupID++) {

		vector<int> &group = mSourceKeyGroups[groupID];
		bool isValid = true;
		for (int nodeID : group) {
			if (replacedSourceSet.find(nodeID) != replacedSourceSet.end()) {
				isValid = false;
				break;
			}
		}
		if (!isValid) continue;

#ifdef OUTPUT_PROGRESS
		cout << "================= Adding group " << (groupID + 1) << " =================" << endl;
#endif

		vector<int> addingList = group;
		vector<int> removingList(0);

		TSolution *newSolution;
		if (!runAddingRemoving(addingList, removingList, currentSolution, newSolution)) return false;
		if (currentSolution != newSolution) {

			vector<int> newSourceContributors = currentSourceContributors;
			vector<int> newTargetContributors = currentTargetContributors;
			for (int nodeID : group) {
				newSourceContributors.push_back(nodeID);
			}
			double newStyleScore;
			if (!computeStyleScore(newSourceContributors, newTargetContributors, newStyleScore)) return false;
#ifdef OUTPUT_PROGRESS
			cout << "Style score = " << newStyleScore << endl;
#endif
			if (newStyleScore < styleScoreThreshold) {
				if (newSolution) delete newSolution;
				newSolution = currentSolution;
			} else {
				mSolutionPool.push_back(newSolution);
				mSolutionSourceContributors.push_back(newSourceContributors);
				mSolutionTargetContributors.push_back(newTargetContributors);

				currentSolution = newSolution;
				currentSourceContributors = newSourceContributors;
				currentTargetContributors = newTargetContributors;

				numAdding++;
				if (numAdding >= maxAdding) break;
			}
		}
	}

	return true;
}

bool ContextPartGraphTabuSearchPart::runAddingRemoving(
	vector<int> &addingList, vector<int> &removingList, 
	TSolution* currentSolution, TSolution* &newSolution)
{

#ifdef OUTPUT_PROGRESS
	cout << "Doing assembly" << endl;
#endif

	// do assembly

	ContextPartGraphAssemble cpga;
	bool sucessFlag;
	if (!cpga.loadGraph(mpNodeGenerator, mpSourceGraph, &currentSolution->mGraph)) return false;
	if (!cpga.loadOldSolution(mSolutionPool[0], currentSolution)) return false;
	if (!cpga.loadAdding(addingList)) return false;
	if (!cpga.loadRemoving(removingList)) return false;
	if (!cpga.process(sucessFlag)) return false;

	// generate solution

	int solutionID = (int)mSolutionPool.size();
	string solutionFolder = mResultFolder + "solution-" + to_string(solutionID) + "/";
	if (!FileUtil::makedir(solutionFolder)) return false; // if there is anything inside this folder
	string assembleVName = solutionFolder + "assemble.ply";
	string similarityName = solutionFolder + "similarity.txt";

	newSolution = new TSolution();
	if (!newSolution->generate(&cpga)) return false;
	if (!newSolution->chain(currentSolution)) return false;

	if (newSolution->mGraph.mRootNode->mChildren.empty()) { // empty graph -- removed all nodes
		if(newSolution) delete newSolution;
		newSolution = currentSolution;
		return true;
	}

	// check functionality compatibility

	double graphSimilarity;
	if (!computeFunctionalSimilarity(newSolution, graphSimilarity)) return false;

#ifdef OUTPUT_PROGRESS
	cout << "Graph score = " << graphSimilarity << endl;
#endif
	if (graphSimilarity < 2.0) { // UNDONE: param functionality compatibility threshold
		if (newSolution) delete newSolution;
		newSolution = currentSolution;
		return true;
	}

	if (!cpga.visualize(solutionFolder)) return false;
	if (!newSolution->visualize(assembleVName, currentSolution)) return false;
	if (!newSolution->saveData(solutionFolder)) return false;
	if (!DataUtil::saveValueListASCII(similarityName, vector<double>(1, graphSimilarity))) return false;

	if (!sucessFlag) {
		if (newSolution) delete newSolution;
		newSolution = currentSolution;
#ifdef OUTPUT_PROGRESS
		cout << "========= Failed to go any further from this solution... =========" << endl;
#endif
		//system("pause");
	}
	//system("pause");

	return true;
}

bool ContextPartGraphTabuSearchPart::runPostProcessing(TSolution* currentSolution, TSolution* &newSolution)
{

	// almost identical to the codes above...

#ifdef OUTPUT_PROGRESS
	cout << "Doing assembly" << endl;
#endif

	// do assembly

	ContextPartGraphAssemble cpga;
	bool sucessFlag;
	if (!cpga.loadGraph(mpNodeGenerator, mpSourceGraph, &currentSolution->mGraph)) return false;
	if (!cpga.loadOldSolution(mSolutionPool[0], currentSolution)) return false;
	if (!cpga.postProcess(sucessFlag)) return false;

	// generate solution

	int solutionID = (int)mSolutionPool.size();
	string solutionFolder = mResultFolder + "solution-" + to_string(solutionID) + "/";
	if (!FileUtil::makedir(solutionFolder)) return false; // if there is anything inside this folder
	string assembleVName = solutionFolder + "assemble.ply";
	string similarityName = solutionFolder + "similarity.txt";

	newSolution = new TSolution();
	if (!newSolution->generate(&cpga)) return false;
	if (!newSolution->chain(currentSolution)) return false;

	if (newSolution->mGraph.mRootNode->mChildren.empty()) { // empty graph -- removed all nodes
		if (newSolution) delete newSolution;
		newSolution = currentSolution;
		return true;
	}

	// check functionality compatibility

	double graphSimilarity;
	if (!computeFunctionalSimilarity(newSolution, graphSimilarity)) return false;

#ifdef OUTPUT_PROGRESS
	cout << "Solution score = " << graphSimilarity << endl;
#endif
	if (graphSimilarity < 0) { // UNDONE: param functionality compatibility threshold
		if (newSolution) delete newSolution;
		newSolution = currentSolution;
		return true;
	}

	if (!cpga.visualize(solutionFolder)) return false;
	if (!newSolution->visualize(assembleVName, currentSolution)) return false;
	if (!newSolution->saveData(solutionFolder)) return false;
	if (!DataUtil::saveValueListASCII(similarityName, vector<double>(1, graphSimilarity))) return false;

	if (!sucessFlag) {
		if (newSolution) delete newSolution;
		newSolution = currentSolution;
#ifdef OUTPUT_PROGRESS
		cout << "========= Failed to go any further from this solution... =========" << endl;
#endif
		//system("pause");
	}
	//system("pause");

	return true;
}

bool ContextPartGraphTabuSearchPart::computeFunctionalSimilarity(TSolution *solution, double &similarity) {

	Eigen::MatrixXd assembleIdentitySimilarity;
	if (true) {
		ContextPartGraphMatch cpgm;
		if (!cpgm.loadWeights(mWeightsFolder)) return false;
		if (!cpgm.loadGraph(solution->mGraph, solution->mGraph)) return false;
		if (!cpgm.process()) return false;
		if (!cpgm.exportSimilarityMatrix(assembleIdentitySimilarity)) return false;
	}
	if (true) {
		ContextPartGraphMatch cpgm;
		if (!cpgm.loadWeights(mWeightsFolder)) return false;
		if (!cpgm.loadGraph(*mpTargetGraph, solution->mGraph)) return false;
		if (!cpgm.process()) return false;
		if (!cpgm.normalizeGraphSimilarity(mOriginalIdentitySimilarity, assembleIdentitySimilarity)) return false;
		if (!cpgm.exportGraphSimilarity(similarity)) return false;
	}

	return true;
}

bool ContextPartGraphTabuSearchPart::computeStyleScore(
	vector<int> &sourceContributors,
	vector<int> &targetContributors,
	double &score)
{
	set<int> sourceSet(sourceContributors.begin(), sourceContributors.end());
	set<int> targetSet(targetContributors.begin(), targetContributors.end());

	score = 2.0;
	for (int sourceID : sourceSet) score -= mSourceContribution[sourceID];
	for (int targetID : targetSet) score -= mTargetContribution[targetID];
	score *= 0.5;

	// re-normalize saliency
	double oldSaliency = 2.0;
	double newSaliency = oldSaliency;
	for (int sourceID : sourceContributors) newSaliency += mSourceSaliency[sourceID];
	for (int targetID : targetContributors) newSaliency -= mTargetSaliency[targetID];
	if(newSaliency > 0) score *= oldSaliency / newSaliency;

	//cout << "Source contributors: ";
	//for (int id : sourceSet) cout << " " << id;
	//cout << endl;
	//cout << "Target contributors: ";
	//for (int id : targetSet) cout << " " << id;
	//cout << endl;

	return true;
}

bool ContextPartGraphTabuSearchPart::cleanUp() {

	for (TSolution *solution : mSolutionPool) {
		delete solution;
	}
	mSolutionPool.clear();

	return true;
}
