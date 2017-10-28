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

#include "ContextPartGraphMatchCurve.h"

#include <fstream>
#include <sstream>
#include <set>

#include "Mesh/MeshUtil.h"
#include "Sample/SampleUtil.h"
#include "Segment/SegmentUtil.h"
#include "Curve/CurveUtil.h"

#include "Context/ContextPartGraph.h"
#include "ContextPartGraphAssembleResult.h"

#include "Match/MatchCurveICP.h"
#include "Match/MatchRigidICP.h"

#include "Data/StyleSynthesisConfig.h"
#include "Data/DataUtil.h"

#include "Utility/PlyExporter.h"

using namespace StyleSynthesis;

ContextPartGraphMatchCurve::ContextPartGraphMatchCurve() {

}

ContextPartGraphMatchCurve::~ContextPartGraphMatchCurve() {

}

bool ContextPartGraphMatchCurve::loadWeights(string weightsFolder) {

	string graphCurveWeightsName = weightsFolder + "weights-graph-curve.txt";
	string graphCurveSigmaName = weightsFolder + "sigma-graph-curve.txt";
	string graphCurveSigmaMultipliersName = weightsFolder + "sigma-multipliers-graph-curve.txt";

	if (!DataUtil::loadValueListASCII(graphCurveWeightsName, mCurveWeights)) return false;
	if (!DataUtil::loadValueListASCII(graphCurveSigmaName, mCurveSigma)) return false;
	if (!DataUtil::loadValueListASCII(graphCurveSigmaMultipliersName, mCurveSigmaMultipliers)) return false;

	return true;
}

bool ContextPartGraphMatchCurve::loadGraph(ContextPartGraph *srcGraph, ContextPartGraph *tgtGraph, ContextPartGraphAssembleResult *solution, Eigen::MatrixXd &similarity) {

	mpSourceGraph = srcGraph;
	mpTargetGraph = tgtGraph;
	mpSolution = solution;
	mNodeSimilarity = similarity;

	// HACK: match nodes only within specific levels
	if (true) {
		Eigen::MatrixXd matSim = mNodeSimilarity;
		matSim.setZero();

		int numSourceNodes = (int)srcGraph->mAllNodes.size();
		int numTargetNodes = (int)tgtGraph->mAllNodes.size();

		vector<int> sourceNodeLevel(numSourceNodes, -1);
		vector<int> targetNodeLevel(numTargetNodes, -1);
		for (int nodeID = 0; nodeID < numSourceNodes; nodeID++) {
			auto node = mpSourceGraph->mAllNodes[nodeID];
			if (node->mParent == mpSourceGraph->mRootNode) sourceNodeLevel[nodeID] = 0;
			else sourceNodeLevel[nodeID] = sourceNodeLevel[node->mParent->mUID] + 1;
		}
		for (int nodeID = 0; nodeID < numTargetNodes; nodeID++) {
			auto node = mpTargetGraph->mAllNodes[nodeID];
			if (node->mParent == mpTargetGraph->mRootNode) targetNodeLevel[nodeID] = 0;
			else targetNodeLevel[nodeID] = targetNodeLevel[node->mParent->mUID] + 1;
		}

		int matchLevels = StyleSynthesisConfig::mContext_MatchNodeLevels;
		for (int sourceID = 0; sourceID < numSourceNodes; sourceID++) {
			if (sourceNodeLevel[sourceID] >= matchLevels) continue;
			for (int targetID = 0; targetID < numTargetNodes; targetID++) {
				if (targetNodeLevel[targetID] >= matchLevels) continue;
				matSim(sourceID, targetID) = mNodeSimilarity(sourceID, targetID);
			}
		}

		mNodeSimilarity.swap(matSim);
	}

	return true;
}

bool ContextPartGraphMatchCurve::process() {

	if (!initialize()) return false;
	if (!gatherCandidates()) return false;
	if (!computeDistance()) return false;
	if (!computeScore()) return false;
	if (!extractPairs()) return false;

	return true;
}

bool ContextPartGraphMatchCurve::initialize() {

	if (!initContext()) return false;
	if (!initCurves()) return false;
	if (!findCliques()) return false;

	return true;
}

bool ContextPartGraphMatchCurve::initContext() {

	// find curve nodes

	cout << "Finding curve nodes" << endl;

	if (!findCurveNodes(mSourceAllCurves, *mpSourceGraph, mSourceAllCurvesNodes)) return false;
	if (!findCurveNodes(mTargetAllCurves, *mpTargetGraph, mTargetAllCurvesNodes)) return false;

	// skip replaced nodes

	if (mpSolution) {
		for (int curveID = 0; curveID < (int)mTargetAllCurves.size(); curveID++) {
			vector<int> updatedNodes(0);
			for (int nodeID : mTargetAllCurvesNodes[curveID]) {
				if (mpSolution->mNodeMapping[nodeID] >= 0) updatedNodes.push_back(nodeID);
			}
			mTargetAllCurvesNodes[curveID].swap(updatedNodes);
		}
	}

	return true;
}

bool ContextPartGraphMatchCurve::initCurves() {

	cout << "Initializing curves data" << endl;

	int numSourceCurves = (int)mSourceAllCurves.size();
	int numTargetCurves = (int)mTargetAllCurves.size();

	// check whether curve is valid

	if (mAllPairIndices.empty()) {
		// process all curves

		mSourceAllCurvesValidFlags.assign(numSourceCurves, true);
		mTargetAllCurvesValidFlags.assign(numTargetCurves, true);

		for (int curveID = 0; curveID < numTargetCurves; curveID++) {
			if (mTargetAllCurvesNodes[curveID].empty()) {
				mTargetAllCurvesValidFlags[curveID] = false; // corresponding part has been replaced
			}
		}

		double minCurveLength = StyleSynthesisConfig::mDeform_MinimumMatchedCurveLength;
		for (int curveID = 0; curveID < numSourceCurves; curveID++) {
			double length;
			if (!CurveUtil::computeCurveLength(mSourceAllCurves[curveID], length)) return false;
			if (length < minCurveLength) mSourceAllCurvesValidFlags[curveID] = false; // curve is too short
		}
		for (int curveID = 0; curveID < numTargetCurves; curveID++) {
			double length;
			if (!CurveUtil::computeCurveLength(mTargetAllCurves[curveID], length)) return false;
			if (length < minCurveLength) mTargetAllCurvesValidFlags[curveID] = false; // curve is too short
		}

	} else {
		// keep all candidate curves

		mSourceAllCurvesValidFlags.assign(numSourceCurves, false);
		mTargetAllCurvesValidFlags.assign(numTargetCurves, false);
		for (vec2i idx : mAllPairIndices) {
			mSourceAllCurvesValidFlags[idx[0]] = true;
			mTargetAllCurvesValidFlags[idx[1]] = true;
		}

	}

	// check whether curve is straight

	mSourceAllCurvesStraightFlags.resize(numSourceCurves, true);
	mTargetAllCurvesStraightFlags.resize(numTargetCurves, true);

	for (int curveID = 0; curveID < numSourceCurves; curveID++) {
		bool straightFlag;
		if (!CurveUtil::checkStraightCurve(mSourceAllCurves[curveID], straightFlag)) return false;
		mSourceAllCurvesStraightFlags[curveID] = straightFlag;
	}
	for (int curveID = 0; curveID < numTargetCurves; curveID++) {
		bool straightFlag;
		if (!CurveUtil::checkStraightCurve(mTargetAllCurves[curveID], straightFlag)) return false;
		mTargetAllCurvesStraightFlags[curveID] = straightFlag;
	}

	// sampling on curves

	double sampleRadius = StyleSynthesisConfig::mCurve_SamplingRadius;

	mCurveMatSP.resize(numSourceCurves);
	mCurveMatSN.resize(numSourceCurves);
	mCurveMatTP.resize(numTargetCurves);
	mCurveMatTN.resize(numTargetCurves);
#pragma omp parallel for
	for (int curveID = 0; curveID < numSourceCurves; curveID++) {
		vector<vec3> curveSamples;
		vector<vec3> curveDirections;
		if (!CurveUtil::sampleLine(mSourceAllCurves[curveID], curveSamples, (float)sampleRadius)) error("sampling curves");
		if (!CurveUtil::computeCurveDirections(curveSamples, curveDirections)) error("curve directions");
		if (!CurveUtil::buildMatrix(curveSamples, mCurveMatSP[curveID])) error("curve matrix");
		if (!CurveUtil::buildMatrix(curveDirections, mCurveMatSN[curveID])) error("curve matrix");
	}
#pragma omp parallel for
	for (int curveID = 0; curveID < numTargetCurves; curveID++) {
		vector<vec3> curveSamples;
		vector<vec3> curveDirections;
		if (!CurveUtil::sampleLine(mTargetAllCurves[curveID], curveSamples, (float)sampleRadius)) error("sampling curves");
		if (!CurveUtil::computeCurveDirections(curveSamples, curveDirections)) error("curve directions");
		if (!CurveUtil::buildMatrix(curveSamples, mCurveMatTP[curveID])) error("curve matrix");
		if (!CurveUtil::buildMatrix(curveDirections, mCurveMatTN[curveID])) error("curve matrix");
	}

	// compute curve relative position

	if (true) {
		vec3 bbMin, bbMax;
		if (!MeshUtil::computeAABB(*mpSourceGraph->mRootNode->mpGraphMesh, bbMin, bbMax)) return false;
		vec3 bbSize = bbMax - bbMin;

		mSourceAllCurveRelPos.resize(numSourceCurves);
#pragma omp parallel for
		for (int curveID = 0; curveID < numSourceCurves; curveID++) {
			Eigen::Vector3d center = mCurveMatSP[curveID].rowwise().mean();
			Eigen::Vector3d &relPos = mSourceAllCurveRelPos[curveID];
			for (int k = 0; k < 3; k++) {
				relPos[k] = (center[k] - bbMin[k]) / bbSize[k];
			}
			relPos = relPos.cwiseMin(1.0);
			relPos = relPos.cwiseMax(0.0);
		}
	}

	if (true) {
		vec3 bbMin, bbMax;
		if (!MeshUtil::computeAABB(*mpTargetGraph->mRootNode->mpGraphMesh, bbMin, bbMax)) return false;
		vec3 bbSize = bbMax - bbMin;

		mTargetAllCurveRelPos.resize(numTargetCurves);
#pragma omp parallel for
		for (int curveID = 0; curveID < numTargetCurves; curveID++) {
			Eigen::Vector3d center = mCurveMatTP[curveID].rowwise().mean();
			Eigen::Vector3d &relPos = mTargetAllCurveRelPos[curveID];
			for (int k = 0; k < 3; k++) {
				relPos[k] = (center[k] - bbMin[k]) / bbSize[k];
			}
			relPos = relPos.cwiseMin(1.0);
			relPos = relPos.cwiseMax(0.0);
		}
	}

	return true;
}

bool ContextPartGraphMatchCurve::findCliques() {

	cout << "Finding symmetric cliques" << endl;

	int numSourceCurves = (int)mSourceAllCurves.size();
	int numTargetCurves = (int)mTargetAllCurves.size();

	// find symmetric source curves

	mSourceAllCurvesCliqueID.assign(numSourceCurves, -1);
	mSourceCliques.clear();
	for (int curveID = 0; curveID < numSourceCurves; curveID++) {
		if ((curveID + 1) % 10 == 0) cout << "\rProcessed " << (curveID + 1) << " / " << numSourceCurves << "          ";
		if (!mSourceAllCurvesValidFlags[curveID]) continue;
		if (mSourceAllCurvesCliqueID[curveID] >= 0) continue;

		vector<int> clique(1, curveID);
		for (int otherID = curveID + 1; otherID < numSourceCurves; otherID++) {
			if (!mSourceAllCurvesValidFlags[otherID]) continue;
			if (mSourceAllCurvesCliqueID[otherID] >= 0) continue;
			bool isSymmetric;
			if (!detectSymmetry(mCurveMatSP[curveID], mCurveMatSN[curveID],
				mCurveMatSP[otherID], mCurveMatSN[otherID], isSymmetric)) return false;
			if (isSymmetric) {
				clique.push_back(otherID);
			}
		}

		if ((int)clique.size() > 1) {
			int cliqueID = (int)mSourceCliques.size();
			mSourceCliques.push_back(clique);
			for (int id : clique) mSourceAllCurvesCliqueID[id] = cliqueID;
			//cout << "Clique: " << clique.size() << endl;
		}
	}
	cout << endl;

	// find symmetric target curves

	// yes those are duplicated codes...
	mTargetAllCurvesCliqueID.assign(numTargetCurves, -1);
	mTargetCliques.clear();
	for (int curveID = 0; curveID < numTargetCurves; curveID++) {
		if ((curveID + 1) % 10 == 0) cout << "\rProcessed " << (curveID + 1) << " / " << numTargetCurves << "          ";
		if (!mTargetAllCurvesValidFlags[curveID]) continue;
		if (mTargetAllCurvesCliqueID[curveID] >= 0) continue;

		vector<int> clique(1, curveID);
		for (int otherID = curveID + 1; otherID < numTargetCurves; otherID++) {
			if (!mTargetAllCurvesValidFlags[otherID]) continue;
			if (mTargetAllCurvesCliqueID[otherID] >= 0) continue;
			bool isSymmetric;
			if (!detectSymmetry(mCurveMatTP[curveID], mCurveMatTN[curveID],
				mCurveMatTP[otherID], mCurveMatTN[otherID], isSymmetric)) return false;
			if (isSymmetric) {
				clique.push_back(otherID);
			}
		}

		if ((int)clique.size() > 1) {
			int cliqueID = (int)mTargetCliques.size();
			mTargetCliques.push_back(clique);
			for (int id : clique) mTargetAllCurvesCliqueID[id] = cliqueID;
			//cout << "Clique: " << clique.size() << endl;
		}
	}
	cout << endl;

	return true;
}

bool ContextPartGraphMatchCurve::gatherCandidates() {

	// find all potential matching pairs

	cout << "Gathering candidates for matching curve pairs" << endl;
	
	mAllPairIndices.clear();
	mAllPairViewDistances.clear();

	// match non-contour to non-contour
	for (int srcID = mSourceCurveIndex[0][0]; srcID < mSourceCurveIndex[1][1]; srcID++) {
		if (!mSourceAllCurvesValidFlags[srcID]) continue;
		for (int tgtID = mTargetCurveIndex[0][0]; tgtID < mTargetCurveIndex[1][1]; tgtID++) {
			if (!mTargetAllCurvesValidFlags[tgtID]) continue;
			// delayed pruning straight curve pairs. same below
			if (mSourceAllCurvesStraightFlags[srcID] && mTargetAllCurvesStraightFlags[tgtID]) continue;
			mAllPairIndices.push_back(vec2i(srcID, tgtID));
			mAllPairViewDistances.push_back(0.0);
		}
	}

	//// match non-contour to contour
	//for (int srcID = mSourceCurveIndex[0][0]; srcID < mSourceCurveIndex[1][1]; srcID++) {
	//	if (!mSourceAllCurvesValidFlags[srcID]) continue;
	//	for (int tgtID = mTargetCurveIndex[2][0]; tgtID < (int)mTargetAllCurves.size(); tgtID++) {
	//		if (!mTargetAllCurvesValidFlags[tgtID]) continue;
	//		//if (mSourceAllCurvesStraightFlags[srcID] && mTargetAllCurvesStraightFlags[tgtID]) continue;
	//		mAllPairIndices.push_back(vec2i(srcID, tgtID));
	//		mAllPairViewDistances.push_back(0.0);
	//	}
	//}

	//// match contour to non-contour
	//for (int srcID = mSourceCurveIndex[2][0]; srcID < (int)mSourceAllCurves.size(); srcID++) {
	//	if (!mSourceAllCurvesValidFlags[srcID]) continue;
	//	for (int tgtID = mTargetCurveIndex[0][0]; tgtID < mTargetCurveIndex[1][1]; tgtID++) {
	//		if (!mTargetAllCurvesValidFlags[tgtID]) continue;
	//		//if (mSourceAllCurvesStraightFlags[srcID] && mTargetAllCurvesStraightFlags[tgtID]) continue;
	//		mAllPairIndices.push_back(vec2i(srcID, tgtID));
	//		mAllPairViewDistances.push_back(0.0);
	//	}
	//}

	// match contour to contour (view dependent)
	int numViews = (int)mCurveViewPoints.size();
	double maxViewDist = 0.6; // neighboring view distance ~= 0.518 for 12x3 viewing configuration
	for (int srcViewID = 0; srcViewID < numViews; srcViewID++) {
		int srcOffID = srcViewID + 2;
		for (int tgtViewID = 0; tgtViewID < numViews; tgtViewID++) {
			int tgtOffID = tgtViewID + 2;
			vec3d viewDiff = mCurveViewPoints[srcViewID] - mCurveViewPoints[tgtViewID];
			double viewDist = viewDiff.length();
			if (viewDist > maxViewDist) continue; // early pruning

			for (int srcID = mSourceCurveIndex[srcOffID][0]; srcID < mSourceCurveIndex[srcOffID][1]; srcID++) {
				if (!mSourceAllCurvesValidFlags[srcID]) continue;
				for (int tgtID = mTargetCurveIndex[tgtOffID][0]; tgtID < mTargetCurveIndex[tgtOffID][1]; tgtID++) {
					if (!mTargetAllCurvesValidFlags[tgtID]) continue;
					if (mSourceAllCurvesStraightFlags[srcID] && mTargetAllCurvesStraightFlags[tgtID]) continue;
					mAllPairIndices.push_back(vec2i(srcID, tgtID));
					mAllPairViewDistances.push_back(viewDist);
				}
			}
		}
	}

	cout << "Number of matching pairs: " << mAllPairIndices.size() << endl;

	return true;
}

bool ContextPartGraphMatchCurve::computeDistance(bool doPruning) {

	// compair curves

	int numMatchPairs = (int)mAllPairIndices.size();

	mAllPairDistances.resize(numMatchPairs);
	mAllPairNodes.resize(numMatchPairs);
	mAllPairNodeScores.resize(numMatchPairs);

	int count = 0;
#pragma omp parallel for
	for (int pairID = 0; pairID < numMatchPairs; pairID++) {
#pragma omp atomic
		count++;
#pragma omp critical
		if (count % 1000 == 0) cout << "\rComparing curve pair " << count << " / " << numMatchPairs << "      ";

		vec2i pairIdx = mAllPairIndices[pairID];
		vector<double> distances;
		int sourceNodeID = -1;
		int targetNodeID = -1;
		double nodeScore = 0;
		if (!computeCurveDistance(pairIdx[0], pairIdx[1], distances, sourceNodeID, targetNodeID, nodeScore, doPruning)) error("compare curve");

		if (!distances.empty()) distances.push_back(mAllPairViewDistances[pairID]); // don't forget to add view distance

		mAllPairDistances[pairID].swap(distances);
		mAllPairNodes[pairID] = vec2i(sourceNodeID, targetNodeID);
		mAllPairNodeScores[pairID] = nodeScore;
	}
	cout << endl;

	return true;
}

bool ContextPartGraphMatchCurve::computeSigma() {

	int numMatchPairs = (int)mAllPairDistances.size();
	int dimDistances = (int)mCurveWeights.size();
	mCurveSigma.resize(dimDistances);

	for (int dim = 0; dim < dimDistances; dim++) {
		vector<double> allDistances(0);
		allDistances.reserve(numMatchPairs);
		for (int pairID = 0; pairID < numMatchPairs; pairID++) {
			if (!mAllPairDistances[pairID].empty()) {
				allDistances.push_back(mAllPairDistances[pairID][dim]);
			}
		}
		int n = (int)(allDistances.size() * 0.5);
		nth_element(allDistances.begin(), allDistances.begin() + n, allDistances.end());
		mCurveSigma[dim] = allDistances[n];
		if (mCurveSigma[dim] == 0) mCurveSigma[dim] = 1.0; // NOTE: be careful!!
	}

	return true;
}

bool ContextPartGraphMatchCurve::computeScore() {

	int numMatchPairs = (int)mAllPairDistances.size();
	int dimDistances = (int)mCurveWeights.size();

	mAllPairScores.resize(numMatchPairs);
#pragma omp parallel for
	for (int pairID = 0; pairID < numMatchPairs; pairID++) {
		double score = mAllPairNodeScores[pairID];
		if (!mAllPairDistances[pairID].empty()) {
			for (int dim = 0; dim < dimDistances; dim++) {
				score += mCurveWeights[dim] * exp(-cml::sqr(mAllPairDistances[pairID][dim] / mCurveSigma[dim]) * mCurveSigmaMultipliers[dim]);
			}
		}
		mAllPairScores[pairID] = score;
	}

	return true;
}

bool ContextPartGraphMatchCurve::computeScoreDerivative() {

	int numMatchPairs = (int)mAllPairDistances.size();
	int dimDistances = (int)mCurveWeights.size();

	mAllPairScoreDerivatives.resize(numMatchPairs, dimDistances*2);
#pragma omp parallel for
	for (int pairID = 0; pairID < numMatchPairs; pairID++) {
		for (int dim = 0; dim < dimDistances; dim++) {
			mAllPairScoreDerivatives(pairID, dim) = exp(-cml::sqr(mAllPairDistances[pairID][dim] / mCurveSigma[dim]) * mCurveSigmaMultipliers[dim]);
			mAllPairScoreDerivatives(pairID, dim + dimDistances) = mAllPairScoreDerivatives(pairID, dim) * mCurveWeights[dim]
				* (-cml::sqr(mAllPairDistances[pairID][dim] / mCurveSigma[dim]));
		}
	}

	return true;
}

bool ContextPartGraphMatchCurve::extractPairs() {

	int numMaximumGroups = 30; // UNDONE: param maximum number of matching groups

	int numSourceCurves = (int)mSourceAllCurves.size();
	int numTargetCurves = (int)mTargetAllCurves.size();
	int numSourceNodes = (int)mpSourceGraph->mAllNodes.size();
	int numTargetNodes = (int)mpTargetGraph->mAllNodes.size();
	int numMatchPairs = (int)mAllPairIndices.size();

	// build matrix

	Eigen::MatrixXd matScore(numSourceCurves, numTargetCurves); // match score : # of all source curves X # of all target curves
	Eigen::MatrixXi matSourceNode(numSourceCurves, numTargetCurves); // contextual source node ID (-1 if curves not matched) : # of all source curves X # of all target curves
	Eigen::MatrixXi matTargetNode(numSourceCurves, numTargetCurves); // contextual target node ID (-1 if curves not matched) : # of all source curves X # of all target curves
	matScore.setZero();
	matSourceNode.setConstant(-1);
	matTargetNode.setConstant(-1);

#pragma omp parallel for
	for (int pairID = 0; pairID < numMatchPairs; pairID++) {
		int sourceCurveID = mAllPairIndices[pairID][0];
		int targetCurveID = mAllPairIndices[pairID][1];
		matScore(sourceCurveID, targetCurveID) = mAllPairScores[pairID];
		matSourceNode(sourceCurveID, targetCurveID) = mAllPairNodes[pairID][0];
		matTargetNode(sourceCurveID, targetCurveID) = mAllPairNodes[pairID][1];
	}

	// extract pairs

	mMatchingsSourceCurves.clear();
	mMatchingsTargetCurves.clear();
	mMatchingsSourceNodes.clear();
	mMatchingsTargetNodes.clear();

	double minMatchingSimilarity = matScore.maxCoeff() * 0.0; // UNDONE: param curve matching pruning threshold

	for (int groupID = 0; groupID < numMaximumGroups; groupID++) {

		// find top matched pair

		int maxSourceID;
		int maxTargetID;
		double maxSimilarity = matScore.maxCoeff(&maxSourceID, &maxTargetID);
		if (maxSimilarity == 0) {
			cout << "No more matched pairs" << endl;
			break; // no more matched pairs
		}
		if (maxSimilarity < minMatchingSimilarity) {
			cout << "No more good candidates" << endl;
			break; // low-quality matching pair
		}

		cout << "Best curve: " << maxSourceID << ", " << maxTargetID << ", " << maxSimilarity;
		cout << " " << (mSourceAllCurvesStraightFlags[maxSourceID] ? "T" : "F");
		cout << (mTargetAllCurvesStraightFlags[maxTargetID] ? "T" : "F") << endl;

		cout << "Extracting matching group " << groupID << endl;

		set<int> usedSourceCurves;
		set<int> usedTargetCurves;

		int maxSourceNode = matSourceNode(maxSourceID, maxTargetID);

		//int maxTargetNode = mAllPairTargetNode(maxSourceID, maxTargetID);
		//for (bool done = false; !done; done = true) {
		
		for (int maxTargetNode : mTargetAllCurvesNodes[maxTargetID]) {

			// find valid nodes to be considered in a symmetry group

			vector<bool> sourceNodesValidFlag(numSourceNodes, false);
			vector<bool> targetNodesValidFlag(numTargetNodes, false);

			sourceNodesValidFlag[maxSourceNode] = true;
			for (auto &node : mpSourceGraph->mAllNodes[maxSourceNode]->mSymmetry) {
				sourceNodesValidFlag[node->mUID] = true;
			}

			targetNodesValidFlag[maxTargetNode] = true;
			for (auto &node : mpTargetGraph->mAllNodes[maxTargetNode]->mSymmetry) {
				targetNodesValidFlag[node->mUID] = true;
			}

			// find all symmetric curves (curves on symmetric nodes only)

			int srcCliqueID = mSourceAllCurvesCliqueID[maxSourceID];
			int tgtCliqueID = mTargetAllCurvesCliqueID[maxTargetID];

			vector<int> srcCurveList;
			vector<int> tgtCurveList;
			vector<int> srcCurveNodes;
			vector<int> tgtCurveNodes;

			if (srcCliqueID >= 0) {
				for (int curveID : mSourceCliques[srcCliqueID]) {
					int curveNode = -1;
					for (int nodeID : mSourceAllCurvesNodes[curveID]) {
						if (sourceNodesValidFlag[nodeID]) {
							curveNode = nodeID;
							break;
						}
					}
					if (curveNode >= 0) {
						srcCurveList.push_back(curveID);
						srcCurveNodes.push_back(curveNode);
					}
				}
			} else {
				srcCurveList.push_back(maxSourceID);
				srcCurveNodes.push_back(maxSourceNode);
			}

			if (tgtCliqueID >= 0) {
				for (int curveID : mTargetCliques[tgtCliqueID]) {
					int curveNode = -1;
					for (int nodeID : mTargetAllCurvesNodes[curveID]) {
						if (targetNodesValidFlag[nodeID]) {
							curveNode = nodeID;
							break;
						}
					}
					if (curveNode >= 0) {
						tgtCurveList.push_back(curveID);
						tgtCurveNodes.push_back(curveNode);
					}
				}
			} else {
				tgtCurveList.push_back(maxTargetID);
				tgtCurveNodes.push_back(maxTargetNode);
			}

			// slice score matrix

			Eigen::MatrixXd sliceScore(srcCurveList.size(), tgtCurveList.size());
			for (int col = 0; col < (int)tgtCurveList.size(); col++) {
				int mapCol = tgtCurveList[col];
				for (int row = 0; row < (int)srcCurveList.size(); row++) {
					int mapRow = srcCurveList[row];
					sliceScore(row, col) = matScore(mapRow, mapCol);
				}
			}

			// skip straight curve matchings (do this after slicing out score matrix)

			if (mSourceAllCurvesStraightFlags[maxSourceID] && mTargetAllCurvesStraightFlags[maxTargetID]) continue;

			// extract best pairs

			vector<vector<vec3>> groupSourceCurves;
			vector<vector<vec3>> groupTargetCurves;
			vector<int> groupSourceNodes;
			vector<int> groupTargetNodes;
			//cout << "Group:" << endl;
			double bestGroupScore = sliceScore.maxCoeff();
			while (true) {
				int row, col;
				double bestScore = sliceScore.maxCoeff(&row, &col);
				if (bestScore == 0) break;
				double thr = StyleSynthesisConfig::mData_CustomNumber2; // HACK: manual number
				//cout << "pair: " << row << " -- " << col << " : " << bestScore << endl;
				if (bestScore <= bestGroupScore * thr) { // UNDONE: param curve matching score threshold within group
					sliceScore.col(col).setZero();
					continue;
				}
				int sourceID = srcCurveList[row];
				int targetID = tgtCurveList[col];
				int sourceNode = srcCurveNodes[row];
				int targetNode = tgtCurveNodes[col];

				// slice score matrix
				sliceScore.col(col).setZero();
				//sliceScore.row(row).setZero();

				usedSourceCurves.insert(srcCurveList[row]);
				usedTargetCurves.insert(tgtCurveList[col]);

				vector<vec3> clampedSourceCurve;
				vector<vec3> clampedTargetCurve;
				if (!clampCurveToNode(
					mSourceAllCurves[sourceID],
					mpSourceGraph->mAllNodes[sourceNode],
					clampedSourceCurve)) return false;
				if (!clampCurveToNode(
					mTargetAllCurves[targetID],
					mpTargetGraph->mAllNodes[targetNode],
					clampedTargetCurve)) return false;
				if (true) {
					// check curve length again
					double minCurveLength = StyleSynthesisConfig::mDeform_MinimumMatchedCurveLength;
					double srcLen, tgtLen;
					if (!CurveUtil::computeCurveLength(clampedSourceCurve, srcLen)) return false;
					if (!CurveUtil::computeCurveLength(clampedTargetCurve, tgtLen)) return false;
					if (srcLen < minCurveLength || tgtLen < minCurveLength) continue;
				}
				
				groupSourceCurves.push_back(clampedSourceCurve);
				groupTargetCurves.push_back(clampedTargetCurve);
				groupSourceNodes.push_back(sourceNode);
				groupTargetNodes.push_back(targetNode);
			}
			if (groupSourceCurves.empty()) continue; // all pairs are pruned

			mMatchingsSourceCurves.push_back(groupSourceCurves);
			mMatchingsTargetCurves.push_back(groupTargetCurves);
			mMatchingsSourceNodes.push_back(groupSourceNodes);
			mMatchingsTargetNodes.push_back(groupTargetNodes);

			//cout << "Target nodes:";
			//for (int id : groupTargetNodes) cout << " " << id;
			//cout << endl;

			cout << "Extracted " << groupSourceCurves.size() << " pairs of curves" << endl;
		}

		// reset original score matrix (target curves only)
		for (int targetID : usedTargetCurves) {
			matScore.col(targetID).setZero();
		}
	}

	return true;
}

bool ContextPartGraphMatchCurve::findCurveNodes(
	vector<vector<vec3>> &inCurves,
	ContextPartGraph &inGraph,
	vector<vector<int>> &outNodes)
{

	cout << "preparing for node finding" << endl;

	TTriangleMesh &graphMesh = *inGraph.mRootNode->mpGraphMesh;
	vector<vector<vector<int>>> &graphSegment = *inGraph.mRootNode->mpGraphSegments;

	// build KD tree

	SKDTree vertexTree;
	SKDTreeData vertexTreeData;
	if (!SampleUtil::buildKdTree(graphMesh.positions, vertexTree, vertexTreeData)) return false;

	// mark node tag for vertex

	int numMeshVertices = (int)graphMesh.positions.size();
	int numNodes = (int)inGraph.mAllNodes.size();
	vector<set<int>> vertexNodeTag(numMeshVertices, set<int>());
	for (int nodeID = 0; nodeID < numNodes; nodeID++) {
		ContextPartGraphNode *node = inGraph.mAllNodes[nodeID];
		for (int faceID : graphSegment[node->mPartLevelID][node->mPartSegmentID]) {
			vec3i faceIdx = graphMesh.indices[faceID];
			for (int k = 0; k < 3; k++) vertexNodeTag[faceIdx[k]].insert(nodeID);
		}
	}

	// process each curve

	int numCurves = (int)inCurves.size();
	outNodes.resize(numCurves);

	int count = 0;
#pragma omp parallel for
	for (int curveID = 0; curveID < numCurves; curveID++) {

#pragma omp atomic
		count++;
#pragma omp critical
		if(count%100==0) cout << "\rFinding nodes for curve " << count << " / " << numCurves << "      ";

		// find associated vertex

		double eps = 1e-3; // UNDONE: param distance threshold to determine whether a vertex is touched by curves
		set<int> touchedVertexSet;
		for (int pointID = 0; pointID < (int)inCurves[curveID].size(); pointID++) {
			vec3 p = inCurves[curveID][pointID];
			SKDT::NamedPoint queryPoint(p[0], p[1], p[2]);
			Thea::BoundedSortedArray<SKDTree::Neighbor> queryResult(10);
			vertexTree.kClosestElements<Thea::MetricL2>(queryPoint, queryResult, eps);
			for (int qID = 0; qID < queryResult.size(); qID++) {
				int vertexID = (int)vertexTree.getElements()[queryResult[qID].getIndex()].id;
				touchedVertexSet.insert(vertexID);
			}
		}

		// vote for nodes

		map<int, int> nodeVotes; // node ID => vote count
		for (int vertexID : touchedVertexSet) {
			for (int nodeID : vertexNodeTag[vertexID]) {
				auto &it = nodeVotes.find(nodeID);
				if (it == nodeVotes.end()) nodeVotes[nodeID] = 1;
				else it->second++;
			}
		}

		// count the votes

		int minVotes = max(2, (int)(inCurves[curveID].size() * 0.5)); // UNDONE: param minimum votes for curve node
		outNodes[curveID].clear();
		for (auto &it : nodeVotes) {
			int nodeID = it.first;
			int numVotes = it.second;
			if (numVotes >= minVotes) {
				outNodes[curveID].push_back(nodeID);
			}
		}
	}
	cout << "\rFinding nodes for curve " << count << " / " << numCurves << "      " << endl;

	if (true) {
		int nodeCount = 0;
		int minCount = INT_MAX;
		int maxCount = 0;
		for (int curveID = 0; curveID < numCurves; curveID++) {
			int count = (int)outNodes[curveID].size();
			nodeCount += count;
			minCount = min(minCount, count);
			maxCount = max(maxCount, count);
		}
		nodeCount /= numCurves;
		cout << "Average node count = " << nodeCount << endl;
		cout << "Min node count = " << minCount << endl;
		cout << "Max node count = " << maxCount << endl;
	}

	return true;
}

bool ContextPartGraphMatchCurve::detectSymmetry(
	Eigen::Matrix3Xd &curveMatSP,
	Eigen::Matrix3Xd &curveMatSN,
	Eigen::Matrix3Xd &curveMatTP,
	Eigen::Matrix3Xd &curveMatTN,
	bool &isSymmetric)
{

	double maxError = cml::sqr(StyleSynthesisConfig::mCurve_SamplingRadius * 1.0); // UNDONE: param max match error for symmetry detection

	isSymmetric = false;

	Eigen::Affine3d xform;
	if (!MatchRigidICP::runRegularShape(10, curveMatSP, curveMatSN, curveMatTP, curveMatTN, xform)) return false;

	double errorST;
	Eigen::Matrix3Xd matXSP = xform * curveMatSP;
	if (!MatchRigidICP::error(matXSP, curveMatTP, errorST)) return false;

	double errorTS;
	Eigen::Matrix3Xd matXTP = xform.inverse() * curveMatTP;
	if (!MatchRigidICP::error(matXTP, curveMatSP, errorTS)) return false;

	if (errorST < maxError && errorTS < maxError) {
		isSymmetric = true;
		return true;
	}
	
	return true;
}

bool ContextPartGraphMatchCurve::computeCurveDistance(
	int sourceCurveID,
	int targetCurveID,
	vector<double> &outDistances,
	int &outSourceNodeID,
	int &outTargetNodeID,
	double &outNodeScore,
	bool doPruning)
{
	
	outDistances.clear();
	outSourceNodeID = -1;
	outTargetNodeID = -1;
	outNodeScore = 0;

	// quickly prune curve by angle
	if (doPruning) {
		double maxCurveAngle = 60.0; // UNDONE: param maximum angle between matched curve direction

		if ((int)mSourceAllCurves[sourceCurveID].size() <= 1) return true;
		if ((int)mTargetAllCurves[targetCurveID].size() <= 1) return true;
		vec3d srcDir = mSourceAllCurves[sourceCurveID].front() - mSourceAllCurves[sourceCurveID].back();
		vec3d tgtDir = mTargetAllCurves[targetCurveID].front() - mTargetAllCurves[targetCurveID].back();
		double srcDirLen = srcDir.length();
		double tgtDirLen = tgtDir.length();
		if (srcDirLen && tgtDirLen) {
			srcDir /= srcDirLen;
			tgtDir /= tgtDirLen;
			double angle = cml::deg(acos(cml::clamp(fabs(cml::dot(srcDir, tgtDir)), 0.0, 1.0)));
			if (angle > maxCurveAngle) return true; // pruned
		} else if (srcDirLen || tgtDirLen) {
			return true; // one curve is a loop while the othe one is not
		} else {
			// loop
		}
	}

	// compute context score

	for (int srcNodeID : mSourceAllCurvesNodes[sourceCurveID]) {
		for (int tgtNodeID : mTargetAllCurvesNodes[targetCurveID]) {
			double nodeScore = mNodeSimilarity(srcNodeID, tgtNodeID);
			if (nodeScore > outNodeScore) {
				outNodeScore = nodeScore;
				outSourceNodeID = srcNodeID;
				outTargetNodeID = tgtNodeID;
			}
		}
	}

	// compute curve matching distances

	if (!computeMatchingDistances(sourceCurveID, targetCurveID, outDistances)) return false;

	return true;
}

bool ContextPartGraphMatchCurve::computeMatchingDistances(
	int sourceCurveID,
	int targetCurveID,
	vector<double> &outDistances)
{

	double maxDist = 1.0; // upper bound for distance value

	Eigen::Matrix3Xd &matSP = mCurveMatSP[sourceCurveID];
	Eigen::Matrix3Xd &matSN = mCurveMatSN[sourceCurveID];
	Eigen::Matrix3Xd &matTP = mCurveMatTP[targetCurveID];
	Eigen::Matrix3Xd &matTN = mCurveMatTN[targetCurveID];

	// run ICP
	bool success = false;
	Eigen::Affine3d transform;
	if (!MatchCurveICP::prealign(matSP, matTP, transform)) return false;
	double prescale = (transform.linear())(0, 0);
	//if (!MatchCurveICP::visualize("ICP-before.ply", matSP, matTP, transform)) return false;
	if (!MatchCurveICP::run(10, matSP, matSN, matTP, matTN, transform, success)) return false;
	//if (!MatchCurveICP::visualize("ICP-after.ply", matSP, matTP, transform)) return false;

	// compute ICP distance
	double icpDist = maxDist;
	if (success) {
		Eigen::Matrix3Xd matXSP = transform * matSP;
		Eigen::Matrix3Xd matXTP = transform.inverse() * matTP;
		double distSqST, distSqTS;
		if (!MatchCurveICP::error(matXSP, matTP, distSqST)) return false;
		if (!MatchCurveICP::error(matXTP, matSP, distSqTS)) return false;
		//cout << "Curve dist: " << distSqST << ", " << distSqTS << endl;
		double distSqICP = (distSqST + distSqTS) / 2;
		icpDist = cml::sqrt_safe(distSqICP);
	}

	// compute arc-length distance
	double lenDist = 0;
	if (true) {
		double srcLen, tgtLen;
		if (!CurveUtil::computeCurveLength(mSourceAllCurves[sourceCurveID], srcLen)) return false;
		if (!CurveUtil::computeCurveLength(mTargetAllCurves[targetCurveID], tgtLen)) return false;
		lenDist = fabs(srcLen - tgtLen);
	}

	// compute position distance
	double posDist = maxDist;
	if (true) {
		Eigen::Vector3d srcRelPos = mSourceAllCurveRelPos[sourceCurveID];
		Eigen::Vector3d tgtRelPos = mTargetAllCurveRelPos[targetCurveID];
		posDist = (srcRelPos - tgtRelPos).norm();
	}

	outDistances.clear();
	outDistances.push_back(icpDist);
	outDistances.push_back(lenDist);
	outDistances.push_back(posDist);

	return true;
}

bool ContextPartGraphMatchCurve::clampCurveToNode(
	vector<vec3> &inCurve,
	ContextPartGraphNode *inNode,
	vector<vec3> &outCurve)
{

	// HACK: skip clamping so that shared curves are consistent across nodes
	outCurve = inCurve;
	return true;

	Eigen::AlignedBox3d bb = inNode->mNodeDescriptors.mBoundingBox;

	// expand bb a little bit
	Eigen::Vector3d padding;
	padding.setConstant(bb.diagonal().norm() * 0.01); // UNDONE: param bb padding size
	bb.extend(bb.min() - padding);
	bb.extend(bb.max() + padding);

	outCurve.clear();
	for (vec3 &point : inCurve) {
		Eigen::Vector3d v(vec3d(point).data());
		if (bb.contains(v)) outCurve.push_back(point);
	}

	return true;
}

bool ContextPartGraphMatchCurve::loadCurve(string srcCurveFolder, string tgtCurveFolder, vector<vec3> &viewPoints) {

	mCurveViewPoints = viewPoints;

	cout << "Loading curves..." << endl;

	string srcRVName = srcCurveFolder + "data-snap-rv.txt";
	string srcBName = srcCurveFolder + "data-snap-b.txt";
	string srcCName = srcCurveFolder + "data-snap-c.txt";

	string tgtRVName = tgtCurveFolder + "data-snap-rv.txt";
	string tgtBName = tgtCurveFolder + "data-snap-b.txt";
	string tgtCName = tgtCurveFolder + "data-snap-c.txt";

	mSourceCurveIndex.clear();
	mTargetCurveIndex.clear();

	vector<vector<vec3>> srcRV, tgtRV;
	if (!CurveUtil::loadCurves(srcRVName, srcRV)) return false;
	if (!CurveUtil::loadCurves(tgtRVName, tgtRV)) return false;
	mSourceCurveIndex.push_back(vec2i((int)mSourceAllCurves.size(), (int)(mSourceAllCurves.size() + srcRV.size())));
	mTargetCurveIndex.push_back(vec2i((int)mTargetAllCurves.size(), (int)(mTargetAllCurves.size() + tgtRV.size())));
	mSourceAllCurves.insert(mSourceAllCurves.end(), srcRV.begin(), srcRV.end());
	mTargetAllCurves.insert(mTargetAllCurves.end(), tgtRV.begin(), tgtRV.end());

	vector<vector<vec3>> srcB, tgtB;
	if (!CurveUtil::loadCurves(srcBName, srcB)) return false;
	if (!CurveUtil::loadCurves(tgtBName, tgtB)) return false;
	mSourceCurveIndex.push_back(vec2i((int)mSourceAllCurves.size(), (int)(mSourceAllCurves.size() + srcB.size())));
	mTargetCurveIndex.push_back(vec2i((int)mTargetAllCurves.size(), (int)(mTargetAllCurves.size() + tgtB.size())));
	mSourceAllCurves.insert(mSourceAllCurves.end(), srcB.begin(), srcB.end());
	mTargetAllCurves.insert(mTargetAllCurves.end(), tgtB.begin(), tgtB.end());

	int numViewPoints;

	ifstream srcCFile(srcCName, ios::binary);
	if (!srcCFile.is_open()) return error("cannot open file " + srcCName);
	srcCFile.read((char*)&numViewPoints, sizeof(numViewPoints));
	if (numViewPoints != (int)mCurveViewPoints.size()) return error("incompatible contour file " + srcCName);
	for (int viewID = 0; viewID < numViewPoints; viewID++) {
		vector<vector<vec3>> srcC;
		if (!CurveUtil::loadCurves(srcCFile, srcC)) return false;
		mSourceCurveIndex.push_back(vec2i((int)mSourceAllCurves.size(), (int)(mSourceAllCurves.size() + srcC.size())));
		mSourceAllCurves.insert(mSourceAllCurves.end(), srcC.begin(), srcC.end());
	}
	srcCFile.close();

	ifstream tgtCFile(tgtCName, ios::binary);
	if (!tgtCFile.is_open()) return error("cannot open file " + tgtCName);
	tgtCFile.read((char*)&numViewPoints, sizeof(numViewPoints));
	if (numViewPoints != (int)mCurveViewPoints.size()) return error("incompatible contour file " + tgtCName);
	for (int viewID = 0; viewID < numViewPoints; viewID++) {
		vector<vector<vec3>> tgtC;
		if (!CurveUtil::loadCurves(tgtCFile, tgtC)) return false;
		mTargetCurveIndex.push_back(vec2i((int)mTargetAllCurves.size(), (int)(mTargetAllCurves.size() + tgtC.size())));
		mTargetAllCurves.insert(mTargetAllCurves.end(), tgtC.begin(), tgtC.end());
	}
	tgtCFile.close();

	return true;
}

bool ContextPartGraphMatchCurve::loadCandidates(vector<vec2i> &curvePairs) {

	// copy candidate curve pairs

	mAllPairIndices = curvePairs;

	// compute view distance

	int numPairs = (int)curvePairs.size();
	int numSourceTypes = (int)mSourceCurveIndex.size();
	int numTargetTypes = (int)mTargetCurveIndex.size();

	mAllPairViewDistances.assign(numPairs, 0.0);
	for (int pairID = 0; pairID < numPairs; pairID++) {
		int sourceID = curvePairs[pairID][0];
		int targetID = curvePairs[pairID][1];
		int srcTypeID, tgtTypeID;
		for (srcTypeID = 0; srcTypeID < numSourceTypes; srcTypeID++) {
			if (sourceID < mSourceCurveIndex[srcTypeID][1]) break;
		}
		for (tgtTypeID = 0; tgtTypeID < numTargetTypes; tgtTypeID++) {
			if (targetID < mTargetCurveIndex[tgtTypeID][1]) break;
		}

		if (srcTypeID >= 2 && tgtTypeID >= 2) {
			mAllPairViewDistances[pairID] = (vec3d(mCurveViewPoints[srcTypeID - 2]) - vec3d(mCurveViewPoints[tgtTypeID - 2])).length();
		}
	}

	return true;
}

bool ContextPartGraphMatchCurve::visualize(string fileName) {

	cout << "Visualizing..." << endl;

	TTriangleMesh &srcMesh = *mpSourceGraph->mRootNode->mpGraphMesh;
	TTriangleMesh &tgtMesh = *mpTargetGraph->mRootNode->mpGraphMesh;

	vec3 srcBBMin, srcBBMax;
	vec3 tgtBBMin, tgtBBMax;
	if (!MeshUtil::computeAABB(srcMesh, srcBBMin, srcBBMax)) return false;
	if (!MeshUtil::computeAABB(tgtMesh, tgtBBMin, tgtBBMax)) return false;
	float hspace = (srcBBMax[0] - srcBBMin[0] + tgtBBMax[0] - tgtBBMin[0]) * 1.2f;
	float vspace = max(srcBBMax[1] - srcBBMin[1], tgtBBMax[1] - tgtBBMin[1]) * 1.2f;
	float stspace = (srcBBMax[0] - tgtBBMin[0]) * 1.2f;

	PlyExporter pe;
	vec3i baseColor(127, 127, 127);
	vec3i nodeColor(255, 255, 255);
	float tubeRadius = 0.01f;

	for (int groupID = 0; groupID < (int)mMatchingsSourceCurves.size(); groupID++) {

		vec3 srcOffset((groupID % 5) * hspace, -(groupID / 5)*vspace, 0.0f);
		vec3 tgtOffset = srcOffset + vec3(stspace, 0.0f, 0.0f);

		// add curve mesh

		int numPairs = (int)mMatchingsSourceCurves[groupID].size();
		for (int pairID = 0; pairID < numPairs; pairID++) {

			TTriangleMesh srcTube;
			TTriangleMesh tgtTube;
			if (!CurveUtil::makeTube(mMatchingsSourceCurves[groupID][pairID], srcTube, tubeRadius)) return false;
			if (!CurveUtil::makeTube(mMatchingsTargetCurves[groupID][pairID], tgtTube, tubeRadius)) return false;

			//vec3i srcColor = SegmentUtil::colorMapping(pairID);
			//vec3i tgtColor = SegmentUtil::colorMapping(pairID);
			vec3i srcColor(255, 0, 0);
			vec3i tgtColor(0, 255, 0);
			if (!pe.addMesh(&srcTube.positions, &srcTube.normals, &srcTube.indices, srcColor, srcOffset)) return false;
			if (!pe.addMesh(&tgtTube.positions, &tgtTube.normals, &tgtTube.indices, tgtColor, tgtOffset)) return false;
		}

		// compute face color

		vector<vec3i> srcFaceColor(srcMesh.indices.size(), baseColor);
		vector<vec3i> tgtFaceColor(tgtMesh.indices.size(), baseColor);

		for (int nodeID : mMatchingsSourceNodes[groupID]) {
			auto &node = mpSourceGraph->mAllNodes[nodeID];
			auto &segment = (*node->mpGraphSegments)[node->mPartLevelID][node->mPartSegmentID];
			for (int faceID : segment) srcFaceColor[faceID] = nodeColor;
		}
		for (int nodeID : mMatchingsTargetNodes[groupID]) {
			auto &node = mpTargetGraph->mAllNodes[nodeID];
			auto &segment = (*node->mpGraphSegments)[node->mPartLevelID][node->mPartSegmentID];
			for (int faceID : segment) tgtFaceColor[faceID] = nodeColor;
		}

		// add base mesh

		if (!pe.addMesh(&srcMesh.positions, &srcMesh.normals, &srcMesh.indices, &srcFaceColor, srcOffset)) return false;
		if (!pe.addMesh(&tgtMesh.positions, &tgtMesh.normals, &tgtMesh.indices, &tgtFaceColor, tgtOffset)) return false;
	}

	if (!pe.output(fileName)) return false;
	
	return true;
}

bool ContextPartGraphMatchCurve::output(
	vector<vector<vector<vec3>>> &srcCurves,
	vector<vector<vector<vec3>>> &tgtCurves,
	vector<vector<int>> &srcNodes,
	vector<vector<int>> &tgtNodes)
{
	srcCurves = mMatchingsSourceCurves;
	tgtCurves = mMatchingsTargetCurves;
	srcNodes = mMatchingsSourceNodes;
	tgtNodes = mMatchingsTargetNodes;

	return true;
}

bool ContextPartGraphMatchCurve::loadCurveDistances(string fileName) {

	ifstream file(fileName);
	if (!file.is_open()) return error("cannot open file " + fileName);

	int numPairs;
	file >> numPairs;
	mAllPairDistances.resize(numPairs);
	mAllPairNodes.resize(numPairs);
	mAllPairNodeScores.resize(numPairs);
	for (int pairID = 0; pairID < numPairs; pairID++) {
		auto &pairDist = mAllPairDistances[pairID];
		int numEntries;
		file >> numEntries;
		pairDist.resize(numEntries);
		for (auto &entry : pairDist) {
			file >> entry;
		}
		vec2i nodePair;
		double nodePairScore;
		file >> nodePair[0] >> nodePair[1] >> nodePairScore;
		mAllPairNodes[pairID] = nodePair;
		mAllPairNodeScores[pairID] = nodePairScore;
	}

	file.close();

	return true;
}

bool ContextPartGraphMatchCurve::saveCurveDistances(string fileName) {

	ofstream file(fileName);
	if (!file.is_open()) return error("cannot write to file " + fileName);

	int numPairs = (int)mAllPairDistances.size();
	file << numPairs << endl;
	for (int pairID = 0; pairID < numPairs; pairID++) {
		auto &pairDist = mAllPairDistances[pairID];
		int numEntries = (int)pairDist.size();
		file << numEntries;
		for (auto &dist : pairDist) {
			file << " " << dist;
		}
		file << endl;
		file << mAllPairNodes[pairID][0] << " " << mAllPairNodes[pairID][1];
		file << " " << mAllPairNodeScores[pairID] << endl;
	}

	file.close();

	return true;
}