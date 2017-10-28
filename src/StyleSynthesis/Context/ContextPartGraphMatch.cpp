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

#include "ContextPartGraphMatch.h"

#include <iostream>
#include <fstream>

#include "Match/MatchRigidICP.h"

#include "Mesh/MeshUtil.h"
#include "Sample/SampleUtil.h"
#include "Segment/SegmentUtil.h"

#include "Utility/PlyExporter.h"
#include "Utility/Timer.h"

#include "Data/DataUtil.h"
#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

//#define OUTPUT_PROGRESS

ContextPartGraphMatch::ContextPartGraphMatch() {

	mNumIterations = StyleSynthesisConfig::mContext_GraphPropagationIteration;
	mPathWeights.resize(mNumIterations + 1, 1.0 / double(mNumIterations+1) );
	mNodeWeights.clear();
	mEdgeWeights.clear();
	mNodeSigma.clear();
	mEdgeSigma.clear();
	mNodeSigmaMultipliers.clear();
	mEdgeSigmaMultipliers.clear();

	mMatchingMode = 0;
}

ContextPartGraphMatch::~ContextPartGraphMatch() {
}

bool ContextPartGraphMatch::loadWeights(string weightsFolder) {

	string graphNodeWeightsName = weightsFolder + "weights-graph-node.txt";
	string graphEdgeWeightsName = weightsFolder + "weights-graph-edge.txt";
	string graphNodeSigmaName = weightsFolder + "sigma-graph-node.txt";
	string graphEdgeSigmaName = weightsFolder + "sigma-graph-edge.txt";
	string graphNodeSigmaMultipliersName = weightsFolder + "sigma-multipliers-graph-node.txt";
	string graphEdgeSigmaMultipliersName = weightsFolder + "sigma-multipliers-graph-edge.txt";
	string graphPathWeightsName = weightsFolder + "weights-graph-path.txt";

	if (!DataUtil::loadValueListASCII(graphNodeWeightsName, mNodeWeights)) return false;
	if (!DataUtil::loadValueListASCII(graphEdgeWeightsName, mEdgeWeights)) return false;
	if (!DataUtil::loadValueListASCII(graphNodeSigmaName, mNodeSigma)) return false;
	if (!DataUtil::loadValueListASCII(graphEdgeSigmaName, mEdgeSigma)) return false;
	if (!DataUtil::loadValueListASCII(graphNodeSigmaMultipliersName, mNodeSigmaMultipliers)) return false;
	if (!DataUtil::loadValueListASCII(graphEdgeSigmaMultipliersName, mEdgeSigmaMultipliers)) return false;
	if (!DataUtil::loadValueListASCII(graphPathWeightsName, mPathWeights)) return false;

	mNumIterations = (int)mPathWeights.size() - 1;

	return true;
}

bool ContextPartGraphMatch::loadGraph(TGraph &source, TGraph &target) {

	mpSourceGraph = &source;
	mpTargetGraph = &target;

	mNumSourceNodes = (int)mpSourceGraph->mAllNodes.size();
	mNumTargetNodes = (int)mpTargetGraph->mAllNodes.size();

	/*
	if (true) {
		cout << "Source nodes:" << endl;
		for (int nodeID = 0; nodeID < mNumSourceNodes; nodeID++) {
			TNode *node = mpSourceGraph->mAllNodes[nodeID];
			int parentID = node->mParent->mUID;
			cout << "\t" << nodeID << ": " << parentID << " | ";
			for (TNode *sym : node->mSymmetry) {
				cout << " " << sym->mUID;
			}
			cout << endl;
		}
		cout << "Target nodes:" << endl;
		for (int nodeID = 0; nodeID < mNumTargetNodes; nodeID++) {
			TNode *node = mpTargetGraph->mAllNodes[nodeID];
			int parentID = node->mParent->mUID;
			cout << "\t" << nodeID << ": " << parentID << " | ";
			for (TNode *sym : node->mSymmetry) {
				cout << " " << sym->mUID;
			}
			cout << endl;
		}
	}
	*/

	return true;
}

bool ContextPartGraphMatch::loadMatchingMode(int mode) {

	mMatchingMode = mode;

	return true;
}

bool ContextPartGraphMatch::process() {

	if (!computeNodeDistance()) return false;
	if (!evaluateNodeSimilarity()) return false;
	if (!computeEdgeDistance()) return false;
	if (!evaluateEdgeSimilarity()) return false;
	if (!computeGraphSimilarity()) return false;

	return true;
}

bool ContextPartGraphMatch::computeNodeDistance() {

	// compute raw distance

#ifdef OUTPUT_PROGRESS
	cout << "Computing node distance" << endl;
#endif

	int numAllNodePairs = mNumSourceNodes * mNumTargetNodes;
	mRawNodeDistance.resize(numAllNodePairs);

#pragma omp parallel for
	for(int pairID=0; pairID < numAllNodePairs; pairID++) {
		int sourceID = pairID / mNumTargetNodes;
		int targetID = pairID % mNumTargetNodes;
		TNode *sourceNode = mpSourceGraph->mAllNodes[sourceID];
		TNode *targetNode = mpTargetGraph->mAllNodes[targetID];

		if (!TNode::computeNodeDistance(sourceNode, targetNode, mRawNodeDistance[pairID])) error("compute node distance");
	}

	if (!mRawNodeDistance.empty() && mNodeWeights.empty()) {
		int dimNodeDistance = (int)mRawNodeDistance[0].size();
		double weight = 1.0 / dimNodeDistance;
		mNodeWeights.assign(dimNodeDistance, weight);
	}
	if (mNodeSigma.empty()) {
		if (!computeNodeSigma()) return false;
	}
	if (mNodeSigmaMultipliers.empty())
	{
		int dimNodeDistance = (int)mRawNodeDistance[0].size();
		mNodeSigmaMultipliers.assign(dimNodeDistance, 1.0);
	}

	return true;
}

bool ContextPartGraphMatch::computeNodeSigma() {

	// compute node distance sigma

#ifdef OUTPUT_PROGRESS
	cout << "Computing node sigma" << endl;
#endif

	int dimNodeDist = (int)mRawNodeDistance[0].size();
	int numAllNodePairs = (int)mRawNodeDistance.size();

	mNodeSigma.resize(dimNodeDist);
	for (int dim = 0; dim < dimNodeDist; dim++) {
		vector<double> allDistances(numAllNodePairs);
		for (int pairID = 0; pairID < numAllNodePairs; pairID++) {
			allDistances[pairID] = mRawNodeDistance[pairID][dim];
		}
		int n = (int)(allDistances.size() * 0.5);
		nth_element(allDistances.begin(), allDistances.begin() + n, allDistances.end());
		mNodeSigma[dim] = allDistances[n];
		if (mNodeSigma[dim] == 0) mNodeSigma[dim] = 1.0; // NOTE: be careful!!
	}

	mEdgeSigma = mNodeSigma; // use node sigma for edge distance if not computed later
	mEdgeSigmaMultipliers = mNodeSigmaMultipliers;  // use node sigma for edge distance if not computed later

	return true;
}

bool ContextPartGraphMatch::evaluateNodeSimilarity() {

	// evaluate node similarity

#ifdef OUTPUT_PROGRESS
	cout << "Evaluating node similarity" << endl;
#endif

	int numAllNodePairs = mNumSourceNodes * mNumTargetNodes;

	mNodeSimilarityMatrix.resize(mNumSourceNodes, mNumTargetNodes);
	for (int pairID = 0; pairID < numAllNodePairs; pairID++) {
#ifdef OUTPUT_PROGRESS
		if ((pairID + 1) % max(1, numAllNodePairs / 10) == 0) cout << "\rProcessed " << (int)((pairID + 1) * 100.0 / numAllNodePairs) << "%         ";
#endif
		int sourceID = pairID / mNumTargetNodes;
		int targetID = pairID % mNumTargetNodes;
		double similarity;
		if (!TNode::computeSimilarity(mRawNodeDistance[pairID], mNodeSigma, mNodeSigmaMultipliers, mNodeWeights, similarity)) return false;
		mNodeSimilarityMatrix(sourceID, targetID) = similarity;
	}
#ifdef OUTPUT_PROGRESS
	cout << endl;
#endif

	return true;
}

bool ContextPartGraphMatch::computeEdgeDistance() {

	int numAllNodePairs = mNumSourceNodes * mNumTargetNodes;

	// extract edge quads

#ifdef OUTPUT_PROGRESS
	cout << "Extracting edge quads" << endl;
#endif

	for (; mMatchingMode < 5; mMatchingMode++) {

		int alterID = mMatchingMode;
		// alternative:
		//		0: all edges
		//		1: remove symmetry edges
		//		2: remove co-centric edges
		//		3: remove adjacent edges
		//		4: remove symmetry+cocentric+adjacent edges
		//		5: fail...

		if (alterID == 5) return error("Graph is too huge...");

		mQuadsForPairs.assign(numAllNodePairs, vector<int>(0));
		mAllEdgeQuads.clear();
		mAllEdgeQuads.reserve(numAllNodePairs * 10);

		bool success = true;
		for (int pairID = 0; pairID < numAllNodePairs; pairID++) { // NOTE: don't parallelize it

			int sourceID = pairID / mNumTargetNodes;
			int targetID = pairID % mNumTargetNodes;
			TNode *sourceNode = mpSourceGraph->mAllNodes[sourceID];
			TNode *targetNode = mpTargetGraph->mAllNodes[targetID];

			int numQuads = (int)mAllEdgeQuads.size();

			// check parent node
			if (true) {
				TNode *srcNeighbor = (sourceNode->mParent != mpSourceGraph->mRootNode) ? sourceNode->mParent : sourceNode;
				TNode *tgtNeighbor = (targetNode->mParent != mpTargetGraph->mRootNode) ? targetNode->mParent : targetNode;
				int srcNbID = srcNeighbor->mUID;
				int tgtNbID = tgtNeighbor->mUID;
				//if (srcNbID != sourceID || tgtNbID != targetID) {
				if (srcNbID != sourceID && tgtNbID != targetID) {
					mAllEdgeQuads.push_back(vec4i(sourceID, targetID, srcNbID, tgtNbID));
				}
			}

			// check children node
			if (true) {
				for (TNode *srcNeighbor : sourceNode->mChildren) {
					for (TNode *tgtNeighbor : targetNode->mChildren) {
						mAllEdgeQuads.push_back(vec4i(sourceID, targetID, srcNeighbor->mUID, tgtNeighbor->mUID));
					}
				}
				/*
				if (sourceNode->mChildren.empty()) {
				for (TNode *tgtNeighbor : targetNode->mChildren) {
				mAllEdgeQuads.push_back(vec4i(sourceID, targetID, sourceID, tgtNeighbor->mUID));
				}
				}
				if (targetNode->mChildren.empty()) {
				for (TNode *srcNeighbor : sourceNode->mChildren) {
				mAllEdgeQuads.push_back(vec4i(sourceID, targetID, srcNeighbor->mUID, targetID));
				}
				}
				*/
			}

			// check symmetry node
			if (alterID != 1 && alterID != 4) {
				for (TNode *srcNeighbor : sourceNode->mSymmetry) {
					for (TNode *tgtNeighbor : targetNode->mSymmetry) {
						mAllEdgeQuads.push_back(vec4i(sourceID, targetID, srcNeighbor->mUID, tgtNeighbor->mUID));
					}
				}
			}

			// check co-centric node
			if (alterID != 2 && alterID != 4) {
				for (TNode *srcNeighbor : sourceNode->mCocentric) {
					for (TNode *tgtNeighbor : targetNode->mCocentric) {
						mAllEdgeQuads.push_back(vec4i(sourceID, targetID, srcNeighbor->mUID, tgtNeighbor->mUID));
					}
				}
			}

			// check adjacent node
			if (alterID != 3 && alterID != 4) {
				for (TNode *srcNeighbor : sourceNode->mAdjacent) {
					for (TNode *tgtNeighbor : targetNode->mAdjacent) {
						mAllEdgeQuads.push_back(vec4i(sourceID, targetID, srcNeighbor->mUID, tgtNeighbor->mUID));
					}
				}
			}

			// check contact node
			if (false) { // HACK: skip this relationship
				for (TNode *srcNeighbor : sourceNode->mContact) {
					for (TNode *tgtNeighbor : targetNode->mContact) {
						mAllEdgeQuads.push_back(vec4i(sourceID, targetID, srcNeighbor->mUID, tgtNeighbor->mUID));
					}
				}
			}

			// check support node
			if (false) { // HACK: skip this relationship
				for (TNode *srcNeighbor : sourceNode->mSupport) {
					for (TNode *tgtNeighbor : targetNode->mSupport) {
						mAllEdgeQuads.push_back(vec4i(sourceID, targetID, srcNeighbor->mUID, tgtNeighbor->mUID));
					}
				}
			}

			// update quads list for node pairs
			for (int quadID = numQuads; quadID < (int)mAllEdgeQuads.size(); quadID++) {
				mQuadsForPairs[pairID].push_back(quadID);
			}

			if (mAllEdgeQuads.size() > 100000000) {
				// graph is too huge... redo it with less edges
				cout << "Graph is too huge... " << mAllEdgeQuads.size() << endl;
				cout << "Pruning graph..." << endl;
				success = false;
				break;
			}
		}

		if (success) break;
	}

	// compute path count

	mPathCountMatrix.resize(mNumIterations + 1);
	mPathCountMatrix[0].resize(mNumSourceNodes, mNumTargetNodes);
	mPathCountMatrix[0].setOnes();

	for (int iterID = 0; iterID < mNumIterations; iterID++) {

#ifdef OUTPUT_PROGRESS
		cout << "\rComputing path count " << (iterID + 1) << " / " << mNumIterations << "       ";
#endif

		EigenMatrixXll newPathCount = mPathCountMatrix[iterID];

#pragma omp parallel for
		for (int pairID = 0; pairID < numAllNodePairs; pairID++) {
			int sourceID = pairID / mNumTargetNodes;
			int targetID = pairID % mNumTargetNodes;
			long long numPaths = 0;
			for (int quadID : mQuadsForPairs[pairID]) {
				vec4i quad = mAllEdgeQuads[quadID];
				int srcNbID = quad[2];
				int tgtNbID = quad[3];
				numPaths += mPathCountMatrix[iterID](srcNbID, tgtNbID);
			}
			newPathCount(sourceID, targetID) = max(1ll, numPaths);
		}

		mPathCountMatrix[iterID + 1] = newPathCount;
		}
#ifdef OUTPUT_PROGRESS
	cout << endl;
#endif

	// compute edge distance

#ifdef OUTPUT_PROGRESS
	cout << "Computing edge distance" << endl;
#endif

	int numAllEdgeQuads = (int)mAllEdgeQuads.size();
	mRawEdgeDistance.resize(numAllEdgeQuads);

#pragma omp parallel for
	for (int quadID = 0; quadID < numAllEdgeQuads; quadID++) {
		vec4i quad = mAllEdgeQuads[quadID];
		TNode* sourceNode = mpSourceGraph->mAllNodes[quad[0]];
		TNode* targetNode = mpTargetGraph->mAllNodes[quad[1]];
		TNode* sourceNeighbor = mpSourceGraph->mAllNodes[quad[2]];
		TNode* targetNeighbor = mpTargetGraph->mAllNodes[quad[3]];

		if (!TNode::computeEdgeDistance(
			sourceNode, sourceNeighbor, targetNode, targetNeighbor,
			mRawEdgeDistance[quadID])) error("compute edge distance");
	}

	if (!mRawEdgeDistance.empty() && mEdgeWeights.empty()) {
		int dimEdgeDistance = (int)mRawEdgeDistance[0].size();
		double weight = 1.0 / dimEdgeDistance;
		mEdgeWeights.assign(dimEdgeDistance, weight);
	}
	if (!mRawEdgeDistance.empty() && mEdgeSigma.empty()) {
		if (!computeEdgeSigma()) return false;
	}
	if (mEdgeSigmaMultipliers.empty()) {
		int dimEdgeDistance = (int)mRawEdgeDistance[0].size();
		mEdgeSigmaMultipliers.assign(dimEdgeDistance, 1.0);
	}

	return true;
}

bool ContextPartGraphMatch::computeEdgeSigma() {

	// compute edge distance sigma

#ifdef OUTPUT_PROGRESS
	cout << "Computing edge sigma" << endl;
#endif

	if (mRawEdgeDistance.empty()) return true; // no edge

	int dimEdgeDist = (int)mRawEdgeDistance[0].size();
	int numAllEdgeQuads = (int)mRawEdgeDistance.size();

	mEdgeSigma.resize(dimEdgeDist);
#pragma omp parallel for
	for (int dim = 0; dim < dimEdgeDist; dim++) {
		vector<double> allDistances(numAllEdgeQuads);
		for (int quadID = 0; quadID < numAllEdgeQuads; quadID++) {
			allDistances[quadID] = mRawEdgeDistance[quadID][dim];
		}
		int n = (int)(allDistances.size() * 0.5);
		nth_element(allDistances.begin(), allDistances.begin() + n, allDistances.end());
		mEdgeSigma[dim] = allDistances[n];
		if (mEdgeSigma[dim] == 0) mEdgeSigma[dim] = 1.0; // NOTE: be careful!!
	}

	return true;
}

bool ContextPartGraphMatch::evaluateEdgeSimilarity() {
	// evaluate edge similarity

#ifdef OUTPUT_PROGRESS
	cout << "Evaluating edge similarity" << endl;
#endif

	int numAllEdgeQuads = (int)mAllEdgeQuads.size();

	mEdgeSimilarityVector.resize(numAllEdgeQuads);
#pragma omp parallel for
	for (int quadID = 0; quadID < numAllEdgeQuads; quadID++) {
		double similarity;
		if (!TNode::computeSimilarity(mRawEdgeDistance[quadID], mEdgeSigma, mEdgeSigmaMultipliers, mEdgeWeights, similarity)) error("edge similarity");
		mEdgeSimilarityVector[quadID] = similarity;
	}

	return true;
}


bool ContextPartGraphMatch::computeGraphSimilarity() {

	int numAllNodePairs = mNumSourceNodes * mNumTargetNodes;

	mGraphSimilarityMatrixPerPathLength.resize(mNumIterations + 1);
	mGraphSimilarityMatrixPerPathLength[0] = mNodeSimilarityMatrix;

	for (int iterID = 0; iterID < mNumIterations; iterID++) {

#ifdef OUTPUT_PROGRESS
		cout << "\rComputing graph similarity " << (iterID + 1) << " / " << mNumIterations << "       ";
#endif
		Eigen::MatrixXd newSimilarityMatrix = mGraphSimilarityMatrixPerPathLength[iterID];

#pragma omp parallel for
		for (int pairID = 0; pairID < numAllNodePairs; pairID++) {
			int sourceID = pairID / mNumTargetNodes;
			int targetID = pairID % mNumTargetNodes;

			double geometrySimilarity = mNodeSimilarityMatrix(sourceID, targetID);
			double pathSimilarity = 0;

			for (int quadID : mQuadsForPairs[pairID]) {
				vec4i quad = mAllEdgeQuads[quadID];
				int srcNbID = quad[2];
				int tgtNbID = quad[3];
				double nodeSimilarity = mGraphSimilarityMatrixPerPathLength[iterID](srcNbID, tgtNbID);
				double edgeSimilarity = mEdgeSimilarityVector(quadID);

				pathSimilarity += nodeSimilarity * edgeSimilarity;
			}

			newSimilarityMatrix(sourceID, targetID) = geometrySimilarity * pathSimilarity;
		}

		mGraphSimilarityMatrixPerPathLength[iterID+1] = newSimilarityMatrix;
	}
#ifdef OUTPUT_PROGRESS
	cout << endl;
#endif

	for (int iterID = 0; iterID <= mNumIterations; iterID++)
		mGraphSimilarityMatrixPerPathLength[iterID].array() /= mPathCountMatrix[iterID].cast<double>().array();

	mGraphSimilarityMatrix = mPathWeights[0] * mGraphSimilarityMatrixPerPathLength[0];
	for (int iterID = 1; iterID <= mNumIterations; iterID++)
		mGraphSimilarityMatrix += mPathWeights[iterID] * mGraphSimilarityMatrixPerPathLength[iterID];

	/*
	if (true) {
		int sourceID1 = (int)StyleSynthesisConfig::mData_CustomNumberList1.values[0];
		int targetID1 = (int)StyleSynthesisConfig::mData_CustomNumberList1.values[1];
		int sourceID2 = (int)StyleSynthesisConfig::mData_CustomNumberList2.values[0];
		int targetID2 = (int)StyleSynthesisConfig::mData_CustomNumberList2.values[1];
		for (int iterID = 0; iterID <= mNumIterations; iterID++) {
			cout << "Iteration " << iterID << ": ";
			cout << mGraphSimilarityMatrixPerPathLength[iterID](sourceID1, targetID1) << "\t";
			cout << mGraphSimilarityMatrixPerPathLength[iterID](sourceID2, targetID2) << endl;
		}
		cout << "Total: ";
		cout << mGraphSimilarityMatrix(sourceID1, targetID1) << "\t";
		cout << mGraphSimilarityMatrix(sourceID2, targetID2) << endl;
	}
	*/

	return true;
}

bool ContextPartGraphMatch::normalizeGraphSimilarity(Eigen::MatrixXd &sourceSimilarity, Eigen::MatrixXd &targetSimilarity) {

	if ((int)sourceSimilarity.rows() != mNumSourceNodes || (int)sourceSimilarity.cols() != mNumSourceNodes) {
		return error("incorrect size of source similarity matrix");
	}

	if ((int)targetSimilarity.rows() != mNumTargetNodes || (int)targetSimilarity.cols() != mNumTargetNodes) {
		return error("incorrect size of target similarity matrix");
	}

	Eigen::VectorXd sourceSelfSimilarity = sourceSimilarity.diagonal();
	Eigen::VectorXd targetSelfSimilarity = targetSimilarity.diagonal();
	if (sourceSelfSimilarity.minCoeff() <= 0 || targetSelfSimilarity.minCoeff() <= 0) {
		return error("non-positive self-similarity");
	}

	Eigen::VectorXd sourceNorm = sourceSelfSimilarity.cwiseInverse().cwiseSqrt();
	Eigen::VectorXd targetNorm = targetSelfSimilarity.cwiseInverse().cwiseSqrt();

	mGraphSimilarityMatrix = sourceNorm.asDiagonal() * mGraphSimilarityMatrix * targetNorm.asDiagonal();

	return true;
}

bool ContextPartGraphMatch::matchGraphNodes() {

#ifdef OUTPUT_PROGRESS
	cout << "Matching graph nodes" << endl;
#endif

	Eigen::MatrixXd matSim = mGraphSimilarityMatrix;
	double similarityThreshold = StyleSynthesisConfig::mContext_MatchNodeSimilarityThreshold;

	// HACK: match nodes only within specific levels
	if (true) {
		matSim.setZero();

		vector<int> sourceNodeLevel(mNumSourceNodes, -1);
		vector<int> targetNodeLevel(mNumTargetNodes, -1);
		for (int nodeID = 0; nodeID < mNumSourceNodes; nodeID++) {
			TNode *node = mpSourceGraph->mAllNodes[nodeID];
			if (node->mParent == mpSourceGraph->mRootNode) sourceNodeLevel[nodeID] = 0;
			else sourceNodeLevel[nodeID] = sourceNodeLevel[node->mParent->mUID] + 1;
		}
		for (int nodeID = 0; nodeID < mNumTargetNodes; nodeID++) {
			TNode *node = mpTargetGraph->mAllNodes[nodeID];
			if (node->mParent == mpTargetGraph->mRootNode) targetNodeLevel[nodeID] = 0;
			else targetNodeLevel[nodeID] = targetNodeLevel[node->mParent->mUID] + 1;
		}

		int matchLevels = StyleSynthesisConfig::mContext_MatchNodeLevels;
		for (int sourceID = 0; sourceID < mNumSourceNodes; sourceID++) {
			if (sourceNodeLevel[sourceID] >= matchLevels) continue;
			for (int targetID = 0; targetID < mNumTargetNodes; targetID++) {
				if (targetNodeLevel[targetID] >= matchLevels) continue;
				matSim(sourceID, targetID) = mGraphSimilarityMatrix(sourceID, targetID);
			}
		}
	}

	vector<vector<vec2i>> allMatchingGroups(0);
	vector<double> allMatchingScores(0);
	while (matSim.sum() > 0) {

		// find top matched node pair

		int maxSourceID;
		int maxTargetID;
		double maxSimilarity = matSim.maxCoeff(&maxSourceID, &maxTargetID);

		vector<vec2i> matchingGroup(1, vec2i(maxSourceID, maxTargetID));
		double matchingScore = mNodeSimilarityMatrix(maxSourceID, maxTargetID);
		if (!sliceOutSimilarityMatrix(maxSourceID, maxTargetID, matSim)) return false;

		// enforce symmetry for matching

		if (true) {
			vector<int> sourceNodeCandidates(0);
			for (TNode* srcNode : mpSourceGraph->mAllNodes[maxSourceID]->mSymmetry) {
				sourceNodeCandidates.push_back(srcNode->mUID);
			}
			sourceNodeCandidates.push_back(maxSourceID);

			for (TNode* tgtNode : mpTargetGraph->mAllNodes[maxTargetID]->mSymmetry) {
				int tgtID = tgtNode->mUID;

				double maxCandidateScore = 0;
				int maxCandidateID = -1;
				for (int srcID : sourceNodeCandidates) {
					double score = matSim(srcID, tgtID);
					if (score > maxCandidateScore) {
						maxCandidateScore = score;
						maxCandidateID = srcID;
					}
				}
				maxSimilarity = max(maxSimilarity, maxCandidateScore);

				if (maxCandidateID >= 0) {
					matchingGroup.push_back(vec2i(maxCandidateID, tgtID));
					matchingScore = max(matchingScore, mNodeSimilarityMatrix(maxCandidateID, tgtID));
					if (!sliceOutSimilarityMatrix(maxCandidateID, tgtID, matSim)) return false;
				}
			}
		}

		if (maxSimilarity < similarityThreshold) continue;

		cout << "Node pair similarity = " << maxSimilarity << endl;

		allMatchingGroups.push_back(matchingGroup);
		allMatchingScores.push_back(maxSimilarity); // graph similarity
		//allMatchingScores.push_back(matchingScore); // node similarity
	}
	int numMatchingGroups = (int)allMatchingGroups.size();

	vector<int> groupOrder(numMatchingGroups);
	for (int k = 0; k < numMatchingGroups; k++) groupOrder[k] = k;
	sort(groupOrder.begin(), groupOrder.end(),
		[&allMatchingScores](int lhs, int rhs){return allMatchingScores[lhs] > allMatchingScores[rhs]; });

	mMatchings.clear();
	for (int k = 0; k < numMatchingGroups; k++) {
		mMatchings.push_back(allMatchingGroups[groupOrder[k]]);
		mMatchingsSimilarity.push_back(allMatchingScores[groupOrder[k]]);
	}

	return true;
}

bool ContextPartGraphMatch::sliceOutSimilarityMatrix(
	int sourceID, int targetID,
	Eigen::MatrixXd &similarityMatrix)
{
	
	// find out source sliced rows
	vector<int> slicedSourceRows(1, sourceID);
	if (true) {

		vector<int> childQueue(1, sourceID);
		int head = 0;
		while (head < (int)childQueue.size()) {
			int nodeID = childQueue[head];
			TNode *node = mpSourceGraph->mAllNodes[nodeID];
			for (TNode *child : node->mChildren) {
				childQueue.push_back(child->mUID);
				slicedSourceRows.push_back(child->mUID);
			}
			head++;
		}

		int parentID = mpSourceGraph->mAllNodes[sourceID]->mParent->mUID;
		while (parentID >= 0) {
			slicedSourceRows.push_back(parentID);
			parentID = mpSourceGraph->mAllNodes[parentID]->mParent->mUID;
		}
	}

	// slice source rows
	for (int row : slicedSourceRows) {
		similarityMatrix.row(row) *= 0.99; // demote a bit after used
	}
	
	// find out target sliced cols
	vector<int> slicedTargetCols(1, targetID);
	if (true) {

		vector<int> childQueue(1, targetID);
		int head = 0;
		while (head < (int)childQueue.size()) {
			int nodeID = childQueue[head];
			TNode *node = mpTargetGraph->mAllNodes[nodeID];
			for (TNode *child : node->mChildren) {
				childQueue.push_back(child->mUID);
				slicedTargetCols.push_back(child->mUID);
			}
			head++;
		}

		int parentID = mpTargetGraph->mAllNodes[targetID]->mParent->mUID;
		while (parentID >= 0) {
			slicedTargetCols.push_back(parentID);
			parentID = mpTargetGraph->mAllNodes[parentID]->mParent->mUID;
		}
	}

	// slice target cols
	for (int col : slicedTargetCols) {
		similarityMatrix.col(col).setZero(); // clear entries
	}

	return true;
}

bool ContextPartGraphMatch::visualize(string fileName) {

	int numTopPairs = 100;

	// visualize matched node pairs

	auto &srcMesh = *(mpSourceGraph->mRootNode->mpGraphMesh);
	auto &tgtMesh = *(mpTargetGraph->mRootNode->mpGraphMesh);
	auto &srcSegment = *(mpSourceGraph->mRootNode->mpGraphSegments);
	auto &tgtSegment = *(mpTargetGraph->mRootNode->mpGraphSegments);

	vec3 srcBBMin, srcBBMax;
	vec3 tgtBBMin, tgtBBMax;
	if (!MeshUtil::computeAABB(srcMesh, srcBBMin, srcBBMax)) return false;
	if (!MeshUtil::computeAABB(tgtMesh, tgtBBMin, tgtBBMax)) return false;
	vec3 bbMin = srcBBMin;
	vec3 bbMax = srcBBMax;
	bbMin.minimize(tgtBBMin);
	bbMax.maximize(tgtBBMax);
	vec3 bbLen = bbMax - bbMin;
	float vspace = bbLen[1] * 2.5f;
	float hspace = bbLen[0] * 1.5f;
	float stspace = bbLen[1] * 1.1f;

	vector<vec2i> matchedPairs(0);
	for (auto &matchingGroup : mMatchings) {
		matchedPairs.insert(matchedPairs.end(), matchingGroup.begin(), matchingGroup.end());
	}

	PlyExporter pe;
	numTopPairs = min(numTopPairs, (int)matchedPairs.size());
	for(int orderID=0; orderID < numTopPairs; orderID++) {
		vec2i pairID = matchedPairs[orderID];
		int sourceID = pairID[0];
		int targetID = pairID[1];
		TNode *srcNode = mpSourceGraph->mAllNodes[sourceID];
		TNode *tgtNode = mpTargetGraph->mAllNodes[targetID];
		int srcLevel = 0;
		if (true) {
			TNode *pNode = srcNode;
			while (pNode->mParent != mpSourceGraph->mRootNode) {
				pNode = pNode->mParent;
				srcLevel++;
			}
		}
		int tgtLevel = 0;
		if (true) {
			TNode *pNode = tgtNode;
			while (pNode->mParent != mpTargetGraph->mRootNode) {
				pNode = pNode->mParent;
				tgtLevel++;
			}
		}

		int rowID = orderID / 10;
		int colID = orderID % 10;
		vec3 srcOffset(colID * hspace, -rowID * vspace, 0.0f);
		vec3 tgtOffset = srcOffset + vec3(0.0f, stspace, 0.0f);
		vec3i srcColor = SegmentUtil::colorMapping(srcLevel);
		vec3i tgtColor = SegmentUtil::colorMapping(tgtLevel);
		//vec3i srcColor = SegmentUtil::colorMapping(0);
		//vec3i tgtColor = SegmentUtil::colorMapping(0);

		vector<vec3i> srcColors(srcMesh.indices.size(), vec3i(127, 127, 127));
		vector<vec3i> tgtColors(tgtMesh.indices.size(), vec3i(127, 127, 127));
		for (int faceID : srcSegment[srcNode->mPartLevelID][srcNode->mPartSegmentID]) srcColors[faceID] = srcColor;
		for (int faceID : tgtSegment[tgtNode->mPartLevelID][tgtNode->mPartSegmentID]) tgtColors[faceID] = tgtColor;

		if(!pe.addMesh(&srcMesh.positions, &srcMesh.normals, &srcMesh.indices, &srcColors, srcOffset)) return false;		
		if(!pe.addMesh(&tgtMesh.positions, &tgtMesh.normals, &tgtMesh.indices, &tgtColors, tgtOffset)) return false;
	}
	if(!pe.output(fileName)) return false;

	return true;
}

bool ContextPartGraphMatch::exportSimilarityMatrix(Eigen::MatrixXd &outMatrix) {

	outMatrix = mGraphSimilarityMatrix;

	return true;
}

bool ContextPartGraphMatch::exportGraphSimilarity(double &similarity) {

	Eigen::MatrixXd matSim = mGraphSimilarityMatrix;

	// HACK: match nodes only within specific levels
	if (true) {
		// compute node levels
		vector<int> sourceNodeLevel(mNumSourceNodes, -1);
		vector<int> targetNodeLevel(mNumTargetNodes, -1);
		for (int nodeID = 0; nodeID < mNumSourceNodes; nodeID++) {
			TNode *node = mpSourceGraph->mAllNodes[nodeID];
			if (node->mParent == mpSourceGraph->mRootNode) sourceNodeLevel[nodeID] = 0;
			else sourceNodeLevel[nodeID] = sourceNodeLevel[node->mParent->mUID] + 1;
		}
		for (int nodeID = 0; nodeID < mNumTargetNodes; nodeID++) {
			TNode *node = mpTargetGraph->mAllNodes[nodeID];
			if (node->mParent == mpTargetGraph->mRootNode) targetNodeLevel[nodeID] = 0;
			else targetNodeLevel[nodeID] = targetNodeLevel[node->mParent->mUID] + 1;
		}

		// get slicing mapping
		int matchLevels = StyleSynthesisConfig::mContext_MatchNodeLevels;
		vector<int> rowMap(0);
		vector<int> colMap(0);
		for (int sourceID = 0; sourceID < mNumSourceNodes; sourceID++) {
			if (sourceNodeLevel[sourceID] < matchLevels) rowMap.push_back(sourceID);
		}
		for (int targetID = 0; targetID < mNumTargetNodes; targetID++) {
			if (targetNodeLevel[targetID] < matchLevels) colMap.push_back(targetID);
		}

		// slice matrix
		matSim.resize((int)rowMap.size(), (int)colMap.size());
		for (int row = 0; row < (int)rowMap.size(); row++) {
			for (int col = 0; col < (int)colMap.size(); col++) {
				matSim(row, col) = mGraphSimilarityMatrix(rowMap[row], colMap[col]);
			}
		}
	}

	//double simST = matSim.rowwise().maxCoeff().mean();
	//double simTS = matSim.colwise().maxCoeff().mean();
	//similarity = (simST + simTS) / 2;

	//similarity = matSim.colwise().maxCoeff().mean();

	similarity = matSim.colwise().maxCoeff().sum();

	return true;
}

bool ContextPartGraphMatch::exportMatchings(vector<vector<vec2i>> &matchings) {

	matchings = mMatchings;

	return true;
}

bool ContextPartGraphMatch::exportMatchingsSimilarity(vector<double> &matchingsSimilarity) {

	matchingsSimilarity = mMatchingsSimilarity;

	return true;
}

bool ContextPartGraphMatch::exportMatchingMode(int &mode) {

	mode = mMatchingMode;

	return true;
}

bool ContextPartGraphMatch::saveMatchings(string fileName, vector<vector<vec2i>> &matchings) {

	ofstream file(fileName);
	if (!file.is_open()) return error("cannot write to file " + fileName);

	int numGroups = (int)matchings.size();
	file << numGroups << endl;
	for (auto &group : matchings) {
		int numPairs = (int)group.size();
		file << numPairs;
		for (vec2i matchedPair : group) {
			file << " " << matchedPair;
		}
		file << endl;
	}

	file.close();

	return true;
}

bool ContextPartGraphMatch::loadMatchings(string fileName, vector<vector<vec2i>> &matchings) {

	ifstream file(fileName);
	if (!file.is_open()) return error("cannot open file " + fileName);

	int numGroups;
	file >> numGroups;
	matchings.resize(numGroups);
	for (auto &group : matchings) {
		int numPairs;
		file >> numPairs;
		group.resize(numPairs);
		for (vec2i &matchedPair : group) {
			file >> matchedPair[0] >> matchedPair[1];
		}
	}

	file.close();

	return true;
}




////////////////////////////////////////
// V addition - export derivatives
////////////////////////////////////////

bool ContextPartGraphMatch::exportDerivatives(vector<Eigen::MatrixXd>& derivative_similarity_matrix_wrt_node_weights, 
											  vector<Eigen::MatrixXd>& derivative_similarity_matrix_wrt_node_sigma_multipliers,
											  vector<Eigen::MatrixXd>& derivative_similarity_matrix_wrt_edge_weights, 
											  vector<Eigen::MatrixXd>& derivative_similarity_matrix_wrt_edge_sigma_multipliers,
											  vector<Eigen::MatrixXd>& derivative_similarity_matrix_wrt_path_weights)
{
	if (mGraphSimilarityMatrixPerPathLength.empty())
	{
		return error("computeGraphSimilarity() should be called before exportDerivatives()");
	}
	// if exportDerivatives was not called before, initialize the derivative vectors & matrices
	if (derivative_similarity_matrix_wrt_node_weights.empty())
	{
		derivative_similarity_matrix_wrt_node_weights.resize(mNodeWeights.size());
		for (int d = 0; d < (int)mNodeWeights.size(); d++)
			derivative_similarity_matrix_wrt_node_weights[d] = Eigen::MatrixXd(mGraphSimilarityMatrix.rows(), mGraphSimilarityMatrix.cols());
	}
	if (derivative_similarity_matrix_wrt_edge_weights.empty())
	{
		derivative_similarity_matrix_wrt_edge_weights.resize(mEdgeWeights.size());
		for (int d = 0; d < (int)mEdgeWeights.size(); d++)
			derivative_similarity_matrix_wrt_edge_weights[d] = Eigen::MatrixXd(mGraphSimilarityMatrix.rows(), mGraphSimilarityMatrix.cols());
	}	
	if (derivative_similarity_matrix_wrt_node_sigma_multipliers.empty())
	{
		derivative_similarity_matrix_wrt_node_sigma_multipliers.resize(mNodeWeights.size());
		for (int d = 0; d < (int)mNodeWeights.size(); d++)
			derivative_similarity_matrix_wrt_node_sigma_multipliers[d] = Eigen::MatrixXd(mGraphSimilarityMatrix.rows(), mGraphSimilarityMatrix.cols());
	}
	if (derivative_similarity_matrix_wrt_edge_sigma_multipliers.empty())
	{
		derivative_similarity_matrix_wrt_edge_sigma_multipliers.resize(mEdgeWeights.size());
		for (int d = 0; d < (int)mEdgeWeights.size(); d++)
			derivative_similarity_matrix_wrt_edge_sigma_multipliers[d] = Eigen::MatrixXd(mGraphSimilarityMatrix.rows(), mGraphSimilarityMatrix.cols());
	}
	if (derivative_similarity_matrix_wrt_path_weights.empty())
	{
		derivative_similarity_matrix_wrt_path_weights.resize(mPathWeights.size());
		for (int d = 0; d < (int)mPathWeights.size(); d++)
			derivative_similarity_matrix_wrt_path_weights[d] = Eigen::MatrixXd(mGraphSimilarityMatrix.rows(), mGraphSimilarityMatrix.cols());
	}

	// initialize all derivatives to 0
	for (int d = 0; d < (int)mNodeWeights.size(); d++)
		derivative_similarity_matrix_wrt_node_weights[d].fill(0.0);
	for (int d = 0; d < (int)mEdgeWeights.size(); d++)
		derivative_similarity_matrix_wrt_edge_weights[d].fill(0.0);
	for (int d = 0; d < (int)mNodeWeights.size(); d++)
		derivative_similarity_matrix_wrt_node_sigma_multipliers[d].fill(0.0);
	for (int d = 0; d < (int)mEdgeWeights.size(); d++)
		derivative_similarity_matrix_wrt_edge_sigma_multipliers[d].fill(0.0);
	for (int d = 0; d < (int)mPathWeights.size(); d++)
		derivative_similarity_matrix_wrt_path_weights[d].fill(0.0);

	int numAllNodePairs = mNumSourceNodes * mNumTargetNodes;

	vector< Eigen::MatrixXd > unnormalized_graph_similarity_per_path_length(mNumIterations + 1);
	for (int iterID = 0; iterID <= mNumIterations; iterID++)
		unnormalized_graph_similarity_per_path_length[iterID] = mGraphSimilarityMatrixPerPathLength[iterID].array() * mPathCountMatrix[iterID].cast<double>().array();


	// ************************
	// *** node derivatives ***
	// ************************
	for (int d = 0; d < (int)mNodeWeights.size(); d++)
	{
		Eigen::MatrixXd derivative_node_similarity_matrix = Eigen::MatrixXd::Zero(mNodeSimilarityMatrix.rows(), mNodeSimilarityMatrix.cols());
		Eigen::MatrixXd derivative_node_similarity_matrix2 = Eigen::MatrixXd::Zero(mNodeSimilarityMatrix.rows(), mNodeSimilarityMatrix.cols());
		for (int nodeID = 0; nodeID < numAllNodePairs; nodeID++)
		{
			int sourceID = nodeID / mNumTargetNodes;
			int targetID = nodeID % mNumTargetNodes;
			derivative_node_similarity_matrix(sourceID, targetID) = exp( -cml::sqr(mRawNodeDistance[nodeID][d] / mNodeSigma[d]) * mNodeSigmaMultipliers[d] );
			derivative_node_similarity_matrix2(sourceID, targetID) = mNodeWeights[d] *
														 		     exp( -cml::sqr(mRawNodeDistance[nodeID][d] / mNodeSigma[d]) * mNodeSigmaMultipliers[d] ) *
  																	    (-cml::sqr(mRawNodeDistance[nodeID][d] / mNodeSigma[d]));
		}

		vector< Eigen::MatrixXd > derivative_similarity_matrix(mNumIterations + 1);
		vector< Eigen::MatrixXd > derivative_similarity_matrix2(mNumIterations + 1);
		derivative_similarity_matrix[0] = derivative_node_similarity_matrix;
		derivative_similarity_matrix2[0] = derivative_node_similarity_matrix2;

		for (int iterID = 0; iterID < mNumIterations; iterID++)
		{
				Eigen::MatrixXd temp_derivative_matrix = derivative_node_similarity_matrix.array()
													   * (unnormalized_graph_similarity_per_path_length[iterID + 1].array() / mNodeSimilarityMatrix.array());
				Eigen::MatrixXd temp_derivative_matrix2 = derivative_node_similarity_matrix2.array()
												       * (unnormalized_graph_similarity_per_path_length[iterID + 1].array() / mNodeSimilarityMatrix.array());

#pragma omp parallel for
			for (int pairID = 0; pairID < numAllNodePairs; pairID++) {
				int sourceID = pairID / mNumTargetNodes;
				int targetID = pairID % mNumTargetNodes;

				double geometrySimilarity = mNodeSimilarityMatrix(sourceID, targetID);
				double pathSimilarity = 0.0;
				double pathSimilarity2 = 0.0;

				for (int quadID : mQuadsForPairs[pairID]) {
					vec4i quad = mAllEdgeQuads[quadID];
					int srcNbID = quad[2];
					int tgtNbID = quad[3];
					double edgeSimilarity = mEdgeSimilarityVector(quadID);

					pathSimilarity += edgeSimilarity * derivative_similarity_matrix[iterID](srcNbID, tgtNbID);
					pathSimilarity2 += edgeSimilarity * derivative_similarity_matrix2[iterID](srcNbID, tgtNbID);
				}

				temp_derivative_matrix(sourceID, targetID) += geometrySimilarity * pathSimilarity;
				temp_derivative_matrix2(sourceID, targetID) += geometrySimilarity * pathSimilarity2;
			}

			derivative_similarity_matrix[iterID + 1] = temp_derivative_matrix;
			derivative_similarity_matrix2[iterID + 1] = temp_derivative_matrix2;
		}


		derivative_similarity_matrix_wrt_node_weights[d] = mPathWeights[0] * (derivative_similarity_matrix[0].array() / mPathCountMatrix[0].cast<double>().array());
		derivative_similarity_matrix_wrt_node_sigma_multipliers[d] = mPathWeights[0] * (derivative_similarity_matrix2[0].array() / mPathCountMatrix[0].cast<double>().array());
		for (int iterID = 1; iterID <= mNumIterations; iterID++)
		{
			derivative_similarity_matrix_wrt_node_weights[d].array() += mPathWeights[iterID] * (derivative_similarity_matrix[iterID].array() / mPathCountMatrix[iterID].cast<double>().array());
			derivative_similarity_matrix_wrt_node_sigma_multipliers[d].array() += mPathWeights[iterID] * (derivative_similarity_matrix2[iterID].array() / mPathCountMatrix[iterID].cast<double>().array());
		}
	}




	 // ************************
	 // *** edge derivatives ***
	 // ************************
	for (int d = 0; d < (int)mEdgeWeights.size(); d++)
	{
		vector< Eigen::MatrixXd > derivative_similarity_matrix(mNumIterations + 1);
		vector< Eigen::MatrixXd > derivative_similarity_matrix2(mNumIterations + 1);
		derivative_similarity_matrix[0] = Eigen::MatrixXd::Zero(mNodeSimilarityMatrix.rows(), mNodeSimilarityMatrix.cols());
		derivative_similarity_matrix2[0] = Eigen::MatrixXd::Zero(mNodeSimilarityMatrix.rows(), mNodeSimilarityMatrix.cols());

		for (int iterID = 0; iterID < mNumIterations; iterID++)
		{
			Eigen::MatrixXd temp_derivative_matrix = Eigen::MatrixXd::Zero(mNodeSimilarityMatrix.rows(), mNodeSimilarityMatrix.cols());
			Eigen::MatrixXd temp_derivative_matrix2 = Eigen::MatrixXd::Zero(mNodeSimilarityMatrix.rows(), mNodeSimilarityMatrix.cols());

#pragma omp parallel for
			for (int pairID = 0; pairID < numAllNodePairs; pairID++) {
				int sourceID = pairID / mNumTargetNodes;
				int targetID = pairID % mNumTargetNodes;

				double geometrySimilarity = mNodeSimilarityMatrix(sourceID, targetID);
				double pathSimilarity = 0.0;
				double pathSimilarity2 = 0.0;

				for (int quadID : mQuadsForPairs[pairID]) {
					vec4i quad = mAllEdgeQuads[quadID];
					int srcNbID = quad[2];
					int tgtNbID = quad[3];
					double edgeSimilarity = mEdgeSimilarityVector(quadID);
					double edgeSimilarityDerivative = exp(-cml::sqr(mRawEdgeDistance[quadID][d] / mEdgeSigma[d]) * mEdgeSigmaMultipliers[d] ); 
					double edgeSimilarityDerivative2 = mEdgeWeights[d] *
													   exp(-cml::sqr(mRawEdgeDistance[quadID][d] / mEdgeSigma[d]) * mEdgeSigmaMultipliers[d]) *
							  						      (-cml::sqr(mRawEdgeDistance[quadID][d] / mEdgeSigma[d]));
					
					pathSimilarity += edgeSimilarityDerivative * unnormalized_graph_similarity_per_path_length[iterID](srcNbID, tgtNbID)
									  + edgeSimilarity * derivative_similarity_matrix[iterID](srcNbID, tgtNbID);
					pathSimilarity2 += edgeSimilarityDerivative2 * unnormalized_graph_similarity_per_path_length[iterID](srcNbID, tgtNbID)
								      + edgeSimilarity * derivative_similarity_matrix2[iterID](srcNbID, tgtNbID);
				}

				temp_derivative_matrix(sourceID, targetID) += geometrySimilarity * pathSimilarity;
				temp_derivative_matrix2(sourceID, targetID) += geometrySimilarity * pathSimilarity2;
			}

			derivative_similarity_matrix[iterID + 1] = temp_derivative_matrix;
			derivative_similarity_matrix2[iterID + 1] = temp_derivative_matrix2;
		}


		derivative_similarity_matrix_wrt_edge_weights[d] = mPathWeights[0] * (derivative_similarity_matrix[0].array() / mPathCountMatrix[0].cast<double>().array());
		derivative_similarity_matrix_wrt_edge_sigma_multipliers[d] = mPathWeights[0] * (derivative_similarity_matrix2[0].array() / mPathCountMatrix[0].cast<double>().array());
		for (int iterID = 1; iterID <= mNumIterations; iterID++)
		{
			derivative_similarity_matrix_wrt_edge_weights[d].array() += mPathWeights[iterID] * (derivative_similarity_matrix[iterID].array() / mPathCountMatrix[iterID].cast<double>().array());
			derivative_similarity_matrix_wrt_edge_sigma_multipliers[d].array() += mPathWeights[iterID] * (derivative_similarity_matrix2[iterID].array() / mPathCountMatrix[iterID].cast<double>().array());
		}
	}


	// ************************
	// *** path derivatives ***
	// ************************
	for (int iterID = 0; iterID <= mNumIterations; iterID++)
		derivative_similarity_matrix_wrt_path_weights[iterID] = mGraphSimilarityMatrixPerPathLength[iterID];

	return true;
}