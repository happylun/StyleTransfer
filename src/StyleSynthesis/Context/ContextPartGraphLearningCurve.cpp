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

#include "ContextPartGraphLearningCurve.h"

#include <fstream>

#include "Mesh/MeshUtil.h"

#include "Segment/SegmentGroupApxCvx.h"

#include "Context/ContextPartGraph.h"
#include "Context/ContextPartGraphMatch.h"
#include "Context/ContextPartGraphMatchCurve.h"

#include "Curve/CurveUtil.h"

#include "Data/StyleSynthesisConfig.h"
#include "Data/DataUtil.h"

#include "Utility/FileUtil.h"
#include "Utility/Timer.h"

using namespace StyleSynthesis;

#define DEBUG_OUTPUT

ContextPartGraphLearningCurve::ContextPartGraphLearningCurve() {
}

ContextPartGraphLearningCurve::~ContextPartGraphLearningCurve() {
}

bool ContextPartGraphLearningCurve::loadData(
	vector<string> &trainScenes,
	vector<vector<string>> &trainMeshes,
	map<string, int> &meshLabelMap,
	vec2i labelPair,
	int &outNumMeshPairs,
	int &outNumCurvePairs)
{

	// load curve pairs data

	cout << "Loading data..." << endl;

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;
	string trainCurveFolder = datasetPrefix + "train-curve/";

	int numScenes = (int)trainScenes.size();
	outNumMeshPairs = 0;
	outNumCurvePairs = 0;
	mMeshPairs.clear();
	mCurvePairs.clear();
	for (int sceneID = 0; sceneID < numScenes; sceneID++) {
		cout << "Loading scene " << trainScenes[sceneID] << "..." << endl;

		string sceneCurveFolder = trainCurveFolder + trainScenes[sceneID] + "/";

		ifstream curvePairsFile(sceneCurveFolder + "curvePairs.txt");

		int numMeshPairs;
		curvePairsFile >> numMeshPairs;
		for (int meshPairID = 0; meshPairID < numMeshPairs; meshPairID++) {
			int meshID1, meshID2;
			int numCurvePairs;
			curvePairsFile >> meshID1 >> meshID2 >> numCurvePairs;

			pair<string, string> meshPair;
			meshPair.first = trainMeshes[sceneID][meshID1];
			meshPair.second = trainMeshes[sceneID][meshID2];

			// load matched curve pairs

			vector<vec2i> matchedCurvePairs(0);
			for (int curvePairID = 0; curvePairID < numCurvePairs; curvePairID++) {
				int curveID1, curveID2;
				curvePairsFile >> curveID1 >> curveID2;
				matchedCurvePairs.push_back(vec2i(curveID1, curveID2));
			}
			if (numCurvePairs <= 0) continue;

			// check label pair

			if (true) {
				auto it1 = meshLabelMap.find(meshPair.first);
				auto it2 = meshLabelMap.find(meshPair.second);
				if (it1 == meshLabelMap.end() || it2 == meshLabelMap.end()) {
					return error("incorrect mesh " + meshPair.first + " & " + meshPair.second);
				}
				int label1 = it1->second;
				int label2 = it2->second;
				if ((label1 != labelPair[0] || label2 != labelPair[1]) &&
					(label2 != labelPair[0] || label1 != labelPair[1]))
				{
					continue; // not targeting label pair -- SKIP
				}
			}

			// retain only a subset of curves if there are too many

			int maxNumCurvePairs = 100;
			if (numCurvePairs > maxNumCurvePairs) {
				random_shuffle(matchedCurvePairs.begin(), matchedCurvePairs.end());
				matchedCurvePairs.resize(maxNumCurvePairs);
				numCurvePairs = maxNumCurvePairs;
			}

			// find third curves

			vector<vec2i> extraCurvePairs(0);
			vec3i numTypeCurves; // only consider curves on second shape
			if (!getNumCurves(meshPair.second, numTypeCurves)) return false;
			for (int curvePairID = 0; curvePairID < numCurvePairs; curvePairID++) {
				int curveID = matchedCurvePairs[curvePairID][1];
				int otherID = -1;
				if (curveID < numTypeCurves[0]) {
					otherID = cml::random_integer(0, numTypeCurves[0]);
				} else if (curveID < numTypeCurves[0] + numTypeCurves[1]) {
					otherID = numTypeCurves[0] + cml::random_integer(0, numTypeCurves[1]);
				} else {
					otherID = numTypeCurves[0] + numTypeCurves[1] + cml::random_integer(0, numTypeCurves[2]);
				}
				extraCurvePairs.push_back(vec2i(matchedCurvePairs[curvePairID][0], otherID));
			}

			vector<vec2i> curvePairs(0);
			curvePairs.reserve(matchedCurvePairs.size() + extraCurvePairs.size());
			curvePairs.insert(curvePairs.end(), matchedCurvePairs.begin(), matchedCurvePairs.end());
			curvePairs.insert(curvePairs.end(), extraCurvePairs.begin(), extraCurvePairs.end());

			mMeshPairs.push_back(meshPair);
			mCurvePairs.push_back(curvePairs);
			outNumMeshPairs++;
			outNumCurvePairs += numCurvePairs;
		}

		curvePairsFile.close();
	}

	return true;
}

bool ContextPartGraphLearningCurve::process() {

	if (!initData()) return false;
	if (!initLearning()) return false;
	if (!runSolver()) return false;

	return true;
}

bool ContextPartGraphLearningCurve::initData() {

	// names

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;
	string weightsFolder = datasetPrefix + "weights/";
	string curveFolder = datasetPrefix + "curve/";

	// init matching

	vector<vec3> curveViewPoints;
	if (!CurveUtil::loadViewPoints(curveFolder + "viewpoint.txt", curveViewPoints)) return false;

	int numMeshPairs = (int)mMeshPairs.size();
	mCurveMatchings.resize(numMeshPairs);

	for (int meshPairID = 0; meshPairID < numMeshPairs; meshPairID++) {
		cout << "====== Initializing data " << (meshPairID + 1) << " / " << numMeshPairs << " ======" << endl;

		string meshName1 = mMeshPairs[meshPairID].first;
		string meshName2 = mMeshPairs[meshPairID].second;
		string curveFolder1 = datasetPrefix + "curve/" + meshName1 + "/";
		string curveFolder2 = datasetPrefix + "curve/" + meshName2 + "/";

		string shortName1 = meshName1.substr(meshName1.find_last_of("/\\") + 1);
		string shortName2 = meshName2.substr(meshName2.find_last_of("/\\") + 1);
		string pairName = shortName1 + "--" + shortName2;
		string curveDistName = datasetPrefix + "CurveLearning/" + pairName + "/curveDist.txt";
		if (!FileUtil::makedir(curveDistName)) return false;

		ContextPartGraphMatchCurve *cpgmc = new ContextPartGraphMatchCurve();
		if (!cpgmc->loadWeights(weightsFolder)) return false;

		if (FileUtil::existsfile(curveDistName)) {
			if (!cpgmc->loadCurveDistances(curveDistName)) return false;
		} else {

			// I guess those data do not require persistent storage...
			TTriangleMesh mesh1, mesh2;
			vector<vector<vector<int>>> segment1, segment2;
			ContextPartGraph graph1, graph2;
			ContextPartGraphNodeGenerator nodeGen;
			Eigen::MatrixXd simMat;
			if (!initGraph(
				meshName1, meshName2,
				mesh1, mesh2, segment1, segment2,
				&graph1, &graph2, &nodeGen, simMat)) return false;

			if (!cpgmc->loadGraph(&graph1, &graph2, 0, simMat)) return false;
			if (!cpgmc->loadCurve(curveFolder1, curveFolder2, curveViewPoints)) return false;
			if (!cpgmc->loadCandidates(mCurvePairs[meshPairID])) return false;
			if (!cpgmc->initialize()) return false;
			if (!cpgmc->computeDistance(false)) return false;
			if (!cpgmc->saveCurveDistances(curveDistName)) return false;
		}

		mCurveMatchings[meshPairID] = cpgmc;

		cout << endl;
	}

	return true;
}

bool ContextPartGraphLearningCurve::initLearning() {

	int numMeshPairs = (int)mMeshPairs.size();
	int dimWeights = (int)mCurveMatchings[0]->mCurveWeights.size();

	// compute sigmas

	if (true) {
		Eigen::MatrixXd allSigmas(numMeshPairs, dimWeights);
		for (int meshPairID = 0; meshPairID < numMeshPairs; meshPairID++) {
			cout << "\rComputing sigma " << (meshPairID + 1) << " / " << numMeshPairs << "      ";
			ContextPartGraphMatchCurve *cpgmc = mCurveMatchings[meshPairID];
			if (!cpgmc->computeSigma()) return false;
			for (int dim = 0; dim < dimWeights; dim++) {
				allSigmas(meshPairID, dim) = cpgmc->mCurveSigma[dim];
			}
		}
		cout << endl;

		mSigmas = allSigmas.colwise().mean();
	}

	// initialize weights
	cout << "Initializing weights" << endl;

	double initialWeight = 1.0 / dimWeights;
	mWeights.resize(dimWeights * 2);
	mWeights << Eigen::VectorXd::Ones(dimWeights) * initialWeight, // curve weights
		Eigen::VectorXd::Ones(dimWeights); // sigma multipliers

	return true;
}

bool ContextPartGraphLearningCurve::runSolver() {

	Eigen::VectorXd lowerBound = 0 * mWeights.array() + 1e-7;
	setLowerBound(lowerBound);
	cppoptlib::LbfgsbSolver<double> solver;
	solver.settings_.maxIter = StyleSynthesisConfig::mContext_FunctionalityLearningIteration;
	solver.settings_.rate = StyleSynthesisConfig::mContext_FunctionalityLearningStepSize;

	cout << "Optimizing..." << endl;
	solver.minimize(*this, mWeights);

	cout << "Solved. loss = " << value(mWeights) << endl;
	cout << "#errors: " << numErrors() << endl;
	cout << "weights = " << mWeights.transpose() << endl;

	for (auto &cm : mCurveMatchings) delete cm; // free up

	return true;
}

bool ContextPartGraphLearningCurve::getNumCurves(string inMeshName, vec3i &outNumCurves) {

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;
	string curveFolder = datasetPrefix + "curve/" + inMeshName + "/";

	string curveRVName = curveFolder + "data-snap-rv.txt";
	string curveBName = curveFolder + "data-snap-b.txt";
	string curveCName = curveFolder + "data-snap-c.txt";

	int numRV, numB, numC;

	ifstream fileRV(curveRVName, ios::binary);
	fileRV.read((char*)(&numRV), sizeof(numRV));
	fileRV.close();

	ifstream fileB(curveBName, ios::binary);
	fileB.read((char*)(&numB), sizeof(numB));
	fileB.close();

	numC = 0;
	ifstream fileC(curveCName, ios::binary);
	int numViews;
	fileC.read((char*)(&numViews), sizeof(numViews));
	for (int viewID = 0; viewID < numViews; viewID++) {
		int numCurves;
		fileC.read((char*)(&numCurves), sizeof(numCurves));
		numC += numCurves;
		for (int curveID = 0; curveID < numCurves; curveID++) {
			int numPoints;
			fileC.read((char*)(&numPoints), sizeof(numPoints));
			fileC.seekg(sizeof(vec3)*numPoints, ios::cur);
		}
	}
	fileC.close();

	outNumCurves = vec3i(numRV, numB, numC);

	return true;
}

bool ContextPartGraphLearningCurve::initGraph(
	string sourceName, string targetName,
	TTriangleMesh &sourceMesh, TTriangleMesh &targetMesh,
	vector<vector<vector<int>>> &sourceSegments,
	vector<vector<vector<int>>> &targetSegments,
	ContextPartGraph *sourceGraph, ContextPartGraph *targetGraph,
	ContextPartGraphNodeGenerator *nodeGen,
	Eigen::MatrixXd &simMat)
{
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

	string sourceShortName = sourceName.substr(sourceName.find_last_of("/\\") + 1);
	string targetShortName = targetName.substr(targetName.find_last_of("/\\") + 1);
	string pairName = sourceShortName + "--" + targetShortName;
	string similarityName = datasetPrefix + "CurveLearning/" + pairName + "/simMat.txt";
	if (!FileUtil::makedir(similarityName)) return false;

	// load graph

	if (!MeshUtil::loadMesh(sourceMeshName, sourceMesh)) return false;
	if (!MeshUtil::loadMesh(targetMeshName, targetMesh)) return false;

	if (!SegmentGroupApxCvx::loadSegments(sourceSegmentName, sourceSegments)) return false;
	if (!SegmentGroupApxCvx::loadSegments(targetSegmentName, targetSegments)) return false;

	if (!sourceGraph->loadGraphHierarchy(sourceGraphHName, *nodeGen, &sourceMesh, &sourceSegments)) return false;
	if (!sourceGraph->loadGraphDescriptor(sourceGraphDName)) return false;
	if (!sourceGraph->loadGraphContext(sourceGraphCName)) return false;

	if (!targetGraph->loadGraphHierarchy(targetGraphHName, *nodeGen, &targetMesh, &targetSegments)) return false;
	if (!targetGraph->loadGraphDescriptor(targetGraphDName)) return false;
	if (!targetGraph->loadGraphContext(targetGraphCName)) return false;

	// match graph

	if (FileUtil::existsfile(similarityName)) {
		if (!DataUtil::loadMatrixASCII(similarityName, simMat)) return false;
	} else {
		// compute self similarity
		Eigen::MatrixXd sourceSelfSimilarity;
		if (true) {
			ContextPartGraphMatch cpgm;
			if (!cpgm.loadWeights(weightsFolder)) return false;
			if (!cpgm.loadGraph(*sourceGraph, *sourceGraph)) return false;
			if (!cpgm.process()) return false;
			if (!cpgm.exportSimilarityMatrix(sourceSelfSimilarity)) return false;
		}
		Eigen::MatrixXd targetSelfSimilarity;
		if (true) {
			ContextPartGraphMatch cpgm;
			if (!cpgm.loadWeights(weightsFolder)) return false;
			if (!cpgm.loadGraph(*targetGraph, *targetGraph)) return false;
			if (!cpgm.process()) return false;
			if (!cpgm.exportSimilarityMatrix(targetSelfSimilarity)) return false;
		}
		// do matching
		vector<vector<vec2i>> matchings;
		vector<double> matchingsSimilarity;
		ContextPartGraphMatch cpgm;
		if (!cpgm.loadWeights(weightsFolder)) return false;
		if (!cpgm.loadGraph(*sourceGraph, *targetGraph)) return false;
		if (!cpgm.process()) return false;
		if (!cpgm.normalizeGraphSimilarity(sourceSelfSimilarity, targetSelfSimilarity)) return false;
		if (!cpgm.exportSimilarityMatrix(simMat)) return false;
		if (!DataUtil::saveMatrixASCII(similarityName, simMat)) return false;
	}

	return true;
}

double ContextPartGraphLearningCurve::value(const Eigen::VectorXd &x) {

	mWeights = x;
	updateWeights();

	double regFactor = StyleSynthesisConfig::mContext_FunctionalityLearningRegularization;
	double objective = 0.0;
	int numErrors = 0;
	int numTriplets = 0;

	int numMeshPairs = (int)mMeshPairs.size();
	for (int meshPairID = 0; meshPairID < numMeshPairs; meshPairID++) {

		int numCurvePairs = (int)mCurvePairs[meshPairID].size() / 2; // first half: matched; second half: extra
		numTriplets += numCurvePairs;

		if (!mCurveMatchings[meshPairID]->computeScore()) error("compute scores");
		Eigen::VectorXd &score = mCurveMatchings[meshPairID]->mAllPairScores;

		Eigen::ArrayXd diff = (score.topRows(numCurvePairs) - score.bottomRows(numCurvePairs)).array();
		diff = diff.cwiseMax(-10.0).cwiseMin(10.0); // clamp outlier value
		numErrors += (int)(diff < 0.0).count();

		objective += ((-diff).exp() + 1.0).log().sum(); // be careful about this...
	}
	objective /= numTriplets;
	objective += x.array().abs().sum() * regFactor; // L1 regularization

#ifdef DEBUG_OUTPUT
	cout << "Objective function called (#loss=" << objective << ", #err=" << numErrors << ")" << endl;
#endif

	return objective;
}

void ContextPartGraphLearningCurve::gradient(Eigen::MatrixXd &dx) {

	double regFactor = StyleSynthesisConfig::mContext_FunctionalityLearningRegularization;

	dx = Eigen::MatrixXd::Zero(mWeights.size(), 1);

	int numTriplets = 0;

	int numMeshPairs = (int)mMeshPairs.size();
	for (int meshPairID = 0; meshPairID < numMeshPairs; meshPairID++) {

		int numCurvePairs = (int)mCurvePairs[meshPairID].size() / 2; // first half: matched; second half: extra
		numTriplets += numCurvePairs;

		if (!mCurveMatchings[meshPairID]->computeScoreDerivative()) error("compute score derivatives");
		Eigen::VectorXd &score = mCurveMatchings[meshPairID]->mAllPairScores;
		Eigen::MatrixXd &deriv = mCurveMatchings[meshPairID]->mAllPairScoreDerivatives;

		Eigen::ArrayXd diffScore = (score.topRows(numCurvePairs) - score.bottomRows(numCurvePairs)).array();
		diffScore = diffScore.cwiseMax(-10.0).cwiseMin(10.0); // clamp outlier value
		Eigen::VectorXd dscore = ((-diffScore).exp() + 1.0).inverse() - 1.0; // be careful about this...

		Eigen::MatrixXd diffDeriv = deriv.topRows(numCurvePairs) - deriv.bottomRows(numCurvePairs);
		dx += diffDeriv.transpose() * dscore;
	}
	dx /= numTriplets;

	Eigen::VectorXd reg = ((mWeights.array() > 0.0).cast<double>() - (mWeights.array() < 0.0).cast<double>()) * regFactor; // L1 regularization
	dx += reg;
}

bool ContextPartGraphLearningCurve::updateWeights() {

	int numMeshPairs = (int)mCurveMatchings.size();
	int dimWeights = (int)mSigmas.size();

	for (int meshPairID = 0; meshPairID < numMeshPairs; meshPairID++) {
		ContextPartGraphMatchCurve *cpgmc = mCurveMatchings[meshPairID];
		for (int dim = 0; dim < dimWeights; dim++) {
			cpgmc->mCurveSigma[dim] = mSigmas[dim];
			cpgmc->mCurveWeights[dim] = mWeights[dim];
			cpgmc->mCurveSigmaMultipliers[dim] = mWeights[dim + dimWeights];
		}
	}

	return true;
}


int ContextPartGraphLearningCurve::numErrors() {
	
	int errorCount = 0;
	int numMeshPairs = (int)mMeshPairs.size();
	for (int meshPairID = 0; meshPairID < numMeshPairs; meshPairID++) {

		int numCurvePairs = (int)mCurvePairs[meshPairID].size() / 2;
		Eigen::VectorXd &score = mCurveMatchings[meshPairID]->mAllPairScores;

		Eigen::ArrayXd diff = (score.topRows(numCurvePairs) - score.bottomRows(numCurvePairs)).array();
		diff = diff.cwiseMax(-10.0).cwiseMin(10.0); // clamp outlier value
		errorCount += (int)(diff < 0.0).count();
	}

	return errorCount;
}

bool ContextPartGraphLearningCurve::computeThreshold(double &threshold) {

	cout << "Computing threshold" << endl;

	vector<double> allScoreList(0);
	int numMeshPairs = (int)mMeshPairs.size();
	for (int meshPairID = 0; meshPairID < numMeshPairs; meshPairID++) {
		Eigen::VectorXd &score = mCurveMatchings[meshPairID]->mAllPairScores;
		vector<double> scoreList(score.data(), score.data() + score.size());
		allScoreList.insert(allScoreList.end(), scoreList.begin(), scoreList.end());
	}

	int n = (int)(allScoreList.size() * 0.5);
	nth_element(allScoreList.begin(), allScoreList.begin() + n, allScoreList.end());
	threshold = allScoreList[n];

	return true;
}

bool ContextPartGraphLearningCurve::outputSigmas(string fileName) {

	vector<double> sigmas(mSigmas.data(), mSigmas.data() + mSigmas.size());
	if (!DataUtil::saveValueListASCII(fileName, sigmas)) return false;

	return true;
}

bool ContextPartGraphLearningCurve::outputWeights(string fileName) {

	vector<double> weights(mWeights.data(), mWeights.data() + mWeights.size());
	if (!DataUtil::saveValueListASCII(fileName, weights)) return false;

	return true;
}