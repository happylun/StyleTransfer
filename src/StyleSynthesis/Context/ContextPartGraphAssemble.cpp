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

#include "ContextPartGraphAssemble.h"

#include <iostream>
#include <fstream>
#include <set>

#include "ContextPartGraph.h"
#include "ContextPartGraphAssembleResult.h"

#include "Mesh/MeshUtil.h"
#include "Sample/SampleUtil.h"
#include "Segment/SegmentUtil.h"

#include "Match/MatchSimpleICP.h"
#include "Match/MatchPrimitiveICP.h"

#include "Utility/PlyExporter.h"

#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

#define OUTPUT_PROGRESS

ContextPartGraphAssemble::ContextPartGraphAssemble() {
	TUtil::gDebugVisualization = StyleSynthesisConfig::mAssemble_DebugVisualization;
}

ContextPartGraphAssemble::~ContextPartGraphAssemble() {
}

bool ContextPartGraphAssemble::loadGraph(TNodeGen *nodeGen, TGraph *source, TGraph *target) {
	mpNodeGenerator = nodeGen;
	mpSourceGraph = source;
	mpTargetGraph = target;
	return true;
}

bool ContextPartGraphAssemble::loadOldSolution(TSolution *initSolution, TSolution *prevSolution) {
	mpInitialSolution = initSolution;
	mpPreviousSolution = prevSolution;
	return true;
}

bool ContextPartGraphAssemble::loadMatching(vector<vec2i> &matching) {
	mMatchPairs = matching;
	return true;
}

bool ContextPartGraphAssemble::loadRemoving(vector<int> &removing) {
	mRemovingNodes = removing;
	return true;
}

bool ContextPartGraphAssemble::loadAdding(vector<int> &adding) {
	mAddingNodes = adding;
	return true;
}

bool ContextPartGraphAssemble::process(bool &outSuccessFlag) {

	outSuccessFlag = false;
	mAlignmentMode = 0;

	if (!initData()) return false;
	if (!findSlots()) return false;

	double maxAlignError = StyleSynthesisConfig::mAssemble_SlotsAlignmentErrorThreshold;
	double maxAddingError = StyleSynthesisConfig::mAssemble_SlotsAddingErrorThreshold;
	double bestAlignError = DBL_MAX;
	double bestAddingError = DBL_MAX;
	vector<Eigen::Affine3d> bestMatchTransformations;
	vector<Eigen::Affine3d> bestNodeTransformations;

	mNumAlignmentMode = max(mNumAlignmentMode, mAlignmentMode + 1);
	for (; mAlignmentMode < mNumAlignmentMode; mAlignmentMode++) {

#ifdef OUTPUT_PROGRESS
		cout << "============== Alignment Mode: " << mAlignmentMode << " ==============" << endl;
#endif

		if (!transformParts()) return false;
		if (!matchSlots()) return false;
		if (!reuseSlots()) return false;

		double alignError;
		if (!alignSlots(alignError)) return false;

		double addingError;
		if (!alignAddings(addingError)) return false;

		if (alignError < bestAlignError) {
			bestAlignError = alignError;
			bestMatchTransformations = mMatchTransformations;
			bestNodeTransformations = mNodeTransformations;
		}
		bestAddingError = min(bestAddingError, addingError);

		if (bestAlignError < maxAlignError && bestAddingError < maxAddingError) {
			if (mNumAlignmentMode < 7 || (mAlignmentMode != 1)) {
				break;
			} // else: checking one principal direction for plane; should also check another one
		}

		if (mUseBestGuess) {
			mUseBestGuess = false;
			mAlignmentMode = 0;
		}

		//system("pause");
	}

	mMatchTransformations.swap(bestMatchTransformations);
	mNodeTransformations.swap(bestNodeTransformations);

	outSuccessFlag = bestAlignError < maxAlignError && bestAddingError < maxAddingError;
#ifdef OUTPUT_PROGRESS
	cout << "Best align error: " << bestAlignError << " vs. " << maxAlignError << endl;
#endif

	return true;
}

bool ContextPartGraphAssemble::postProcess(bool &outSuccessFlag) {

	//TUtil::gDebugVisualization = true; // HACK: just debug post-process

	// call this after removing and before adding

	outSuccessFlag = false;
	mAlignmentMode = 0;

	double maxAlignError = StyleSynthesisConfig::mAssemble_SlotsAlignmentErrorThreshold;

	if (!initData()) return false;

	for (int nodeID : mTargetWorkingNodes) {
		if (mpPreviousSolution->mReplaceMapping[nodeID][0] < 0) {
			cout << "Not all nodes are replaced" << endl;
			cout << "Node " << nodeID << " is not replaced" << endl;
			return true; // not all nodes are replaced -- skip post-alignment
		}
	}

	if (!findSlots()) return false;

	double alignError;
	if (!postAlignSlots(alignError)) return false;

	outSuccessFlag = alignError < maxAlignError;
#ifdef OUTPUT_PROGRESS
	cout << "Post align error: " << alignError << " vs. " << maxAlignError << endl;
#endif

	return true;
}

bool ContextPartGraphAssemble::initData() {

#ifdef OUTPUT_PROGRESS
	cout << "Initializing data..." << endl;
#endif

	int numMatchPairs = (int)mMatchPairs.size();
	int numSourceNodes = (int)mpSourceGraph->mAllNodes.size();
	int numTargetNodes = (int)mpTargetGraph->mAllNodes.size();

	mMatchTransformations.assign(numMatchPairs, Eigen::Affine3d::Identity());
	mNodeTransformations.assign(numTargetNodes, Eigen::Affine3d::Identity());

	// mark node type

	if (true) {
		// get current replacing list
		set<int> srcReplaceSet;
		set<int> tgtReplaceSet;
		for (vec2i matchPair : mMatchPairs) {
			srcReplaceSet.insert(matchPair[0]);
			tgtReplaceSet.insert(matchPair[1]);
		}
		// get previous replacing list
		vector<int> prevReplaceList;
		for (int nodeID = 0; nodeID < numTargetNodes; nodeID++) {
			int sourceID = mpPreviousSolution->mReplaceMapping[nodeID][0];
			if (sourceID >= 0) {
				prevReplaceList.push_back(nodeID);
				srcReplaceSet.insert(sourceID);
			}
		}
		vector<int> srcReplaceList(srcReplaceSet.begin(), srcReplaceSet.end());
		vector<int> tgtReplaceList(tgtReplaceSet.begin(), tgtReplaceSet.end());
		// mark node type
		if (!TUtil::markNodeTypes(
			mpSourceGraph,
			srcReplaceList,
			vector<int>(0), // no removing
			vector<int>(0), // no previous operations
			mSourceNodeTypes)) return false;
		if (!TUtil::markNodeTypes(
			mpTargetGraph,
			tgtReplaceList,
			mRemovingNodes,
			prevReplaceList,
			mTargetNodeTypes)) return false;

		// push back working node IDs
		// working nodes = replaced nodes + active nodes
		// replaced node always precedes active node (for better global alignment order)

		mSourceWorkingNodes.clear();
		for (int nodeID = 0; nodeID < numSourceNodes; nodeID++) {
			if (mSourceNodeTypes[nodeID] == TUtil::TNodeType_REPLACED) mSourceWorkingNodes.push_back(nodeID);
		}
		for (int nodeID = 0; nodeID < numSourceNodes; nodeID++) {
			if (mSourceNodeTypes[nodeID] == TUtil::TNodeType_ACTIVE) mSourceWorkingNodes.push_back(nodeID);
		}

		mTargetWorkingNodes.clear();
		for (int nodeID = 0; nodeID < numTargetNodes; nodeID++) {
			if (mTargetNodeTypes[nodeID] == TUtil::TNodeType_REPLACED) mTargetWorkingNodes.push_back(nodeID);
		}
		for (int nodeID = 0; nodeID < numTargetNodes; nodeID++) {
			if (mTargetNodeTypes[nodeID] == TUtil::TNodeType_ACTIVE) mTargetWorkingNodes.push_back(nodeID);
		}
	}

	// extract working nodes for part adding

	if (true) {
		set<int> prevReplaceNodeSet;
		for (vec2i replacePair : mpPreviousSolution->mReplaceMapping) {
			if(replacePair[1] >= 0) prevReplaceNodeSet.insert(replacePair[0]);
		}
		mSourceAddingWorkingNodes.clear();
		set<int> addingWorkingNodeSet;
		for (int nodeID : mAddingNodes) {
			mSourceAddingWorkingNodes.push_back(nodeID);
			for (TNode *adjNode : mpSourceGraph->mAllNodes[nodeID]->mAdjacent) {
				int adjID = adjNode->mUID;
				if (prevReplaceNodeSet.find(adjID) != prevReplaceNodeSet.end()) {
					addingWorkingNodeSet.insert(adjID);
				}
			}
		}
		for (int nodeID : addingWorkingNodeSet) {
			mSourceAddingWorkingNodes.push_back(nodeID);
		}
	}

	// extract node meshes/samples

	mSourceNodeSamples.resize(numSourceNodes);
	mTargetNodeSamples.resize(numTargetNodes);

	//cout << "Source working nodes:";
	for (int nodeID : mSourceWorkingNodes) {
		//cout << " " << nodeID;
		TTriangleMesh nodeMesh;
		if (!TUtil::extractNodeMesh(mpSourceGraph->mAllNodes[nodeID], nodeMesh)) return false;
		if (!TUtil::extractMeshSamples(nodeMesh, mSourceNodeSamples[nodeID])) return false;
	}
	//cout << endl;
	//cout << "Target working nodes:";
	for (int nodeID : mTargetWorkingNodes) {
		//cout << " " << nodeID;
		TTriangleMesh nodeMesh;
		if (!TUtil::extractNodeMesh(mpTargetGraph->mAllNodes[nodeID], nodeMesh)) return false;
		if (!TUtil::extractMeshSamples(nodeMesh, mTargetNodeSamples[nodeID])) return false;
	}
	//cout << endl;
	for (int nodeID : mSourceAddingWorkingNodes) {
		TTriangleMesh nodeMesh;
		if (!TUtil::extractNodeMesh(mpSourceGraph->mAllNodes[nodeID], nodeMesh)) return false;
		if (!TUtil::extractMeshSamples(nodeMesh, mSourceNodeSamples[nodeID])) return false;
	}

	// compute number of alignment modes

	if (true) {
		set<int> primitiveSet;
		for (vec2i matchPair : mMatchPairs) {
			int primitive = mpSourceGraph->mAllNodes[matchPair[0]]->mNodeDescriptors.mPrimitive;
			primitiveSet.insert(primitive);
		}
		for (int sourceID : mAddingNodes) {
			int primitive = mpSourceGraph->mAllNodes[sourceID]->mNodeDescriptors.mPrimitive;
			primitiveSet.insert(primitive);
		}

		mAlignmentMode = StyleSynthesisConfig::mAssemble_InitialAlignmentMode;
		mUseBestGuess = mAlignmentMode < 0;
		auto &bestGuessModes = StyleSynthesisConfig::mAssemble_BestGuessAlignmentMode.values;
		mNumAlignmentMode = 0;
		for (int primitive : primitiveSet) {
			int numModes = 0;
			if (primitive == 0) { // stick
				numModes = 6;
			} else if (primitive == 1) { // plane
				numModes = 8;
			} else { // sphere
				numModes = 5;
			}
			if (mUseBestGuess && primitive < (int)bestGuessModes.size()) {
				mAlignmentMode = bestGuessModes[primitive];
			}
			if (numModes > mNumAlignmentMode) {
				mNumAlignmentMode = numModes;
#ifdef OUTPUT_PROGRESS
				cout << "Primitive: " << primitive << endl;
#endif
			}
		}
		if (mNumAlignmentMode == 0) mNumAlignmentMode = 8; // use max by default
	}

	return true;
}

bool ContextPartGraphAssemble::findSlots() {

#ifdef OUTPUT_PROGRESS
	cout << "Finding slots..." << endl;
#endif

	if (!TUtil::extractNodeSlots(
		mpSourceGraph,
		mSourceNodeSamples,
		mSourceWorkingNodes,
		mSourceNodeSlots)) return false;

	if (!TUtil::extractNodeSlots(
		mpTargetGraph,
		mTargetNodeSamples,
		mTargetWorkingNodes,
		mTargetNodeSlots)) return false;

	if (!TUtil::extractNodeSlots(
		mpSourceGraph,
		mSourceNodeSamples,
		mSourceAddingWorkingNodes,
		mAddingNodeSlots)) return false;

	if (TUtil::gDebugVisualization) {
		if (!visualizeSlots(
			"sourceSlots.ply",
			mSourceWorkingNodes,
			mSourceNodeSlots)) return false;
		if (!visualizeSlots(
			"targetSlots.ply",
			mTargetWorkingNodes,
			mTargetNodeSlots)) return false;
		if (!visualizeSlots(
			"addingSlots.ply",
			mSourceAddingWorkingNodes,
			mAddingNodeSlots)) return false;
	}

	return true;
}

bool ContextPartGraphAssemble::transformParts() {

	int numMatchPairs = (int)mMatchPairs.size();

	if (numMatchPairs == 0) return true;

	// compute transformation between matched parts

	mMatchTransformations.resize(numMatchPairs);
	mFittingTransformations.resize(numMatchPairs);

	for (int matchID = 0; matchID < numMatchPairs; matchID++) {
#ifdef OUTPUT_PROGRESS
		cout << "Applying transformation " << matchID << endl;
#endif
		int sourceID = mMatchPairs[matchID][0];
		int targetID = mMatchPairs[matchID][1];

		Eigen::Affine3d matchXform;
		Eigen::Affine3d fittingXform;
		if (!TUtil::alignMatchedParts(
			mAlignmentMode,
			mpSourceGraph->mAllNodes[sourceID],
			mpTargetGraph->mAllNodes[targetID],
			mSourceNodeSamples[sourceID],
			mTargetNodeSamples[targetID],
			matchXform, fittingXform)) return false;

		mMatchTransformations[matchID] = matchXform;
		mFittingTransformations[matchID] = fittingXform;
	}

	// enforce symetry for part transformation

	if (false) {

		double eps = StyleSynthesisConfig::mAssemble_SlotsSamplingRadius * 2.0;

		Eigen::Matrix3Xd matTP, matTN;
		if (true) {
			int pivotID = 0; // use first target part in the matching list as pivot
			vec2i pivotPair = mMatchPairs[pivotID];
			Eigen::Affine3d pivotXform = mMatchTransformations[pivotID];
			if (!SampleUtil::buildMatrices(mSourceNodeSamples[pivotPair[0]], matTP, matTN)) return false;
			matTP = pivotXform * matTP;
			matTN = pivotXform.rotation() * matTN;
		}

		for (int pairID = 1; pairID < numMatchPairs; pairID++) {
#ifdef OUTPUT_PROGRESS
			cout << "Enforcing symmetry " << pairID << endl;
#endif

			Eigen::Matrix3Xd matSP, matSN;
			if (true) {
				vec2i otherPair = mMatchPairs[pairID];
				Eigen::Affine3d otherXform = mMatchTransformations[pairID];
				if (!SampleUtil::buildMatrices(mSourceNodeSamples[otherPair[0]], matSP, matSN)) return false;
				matSP = otherXform * matSP;
				matSN = otherXform.rotation() * matSN;
			}

			Eigen::Affine3d xformICP, xformReg;
			xformICP.setIdentity();
			if (!MatchSimpleICP::run(20, matSP, matSN, matTP, matTN, xformICP)) return false;

			if (!TUtil::regularizeTransformation(xformICP, xformReg, eps)) return false;

			Eigen::Affine3d xformComp = xformReg.inverse() * xformICP;
			mMatchTransformations[pairID] = xformComp * mMatchTransformations[pairID];
		}
	}

	return true;
}

bool ContextPartGraphAssemble::matchSlots() {

	mAlignedNodeSlots = mTargetNodeSlots;

	int numMatchPairs = (int)mMatchPairs.size();
	if (numMatchPairs == 0) return true;

	// align slots between matched parts

	for (int matchID = 0; matchID < numMatchPairs; matchID++) {
#ifdef OUTPUT_PROGRESS
		cout << "Matching slots " << matchID << endl;
#endif
		int sourceID = mMatchPairs[matchID][0];
		int targetID = mMatchPairs[matchID][1];

		// apply transformation
		Eigen::Affine3d &fittingXform = mFittingTransformations[matchID];
		Eigen::Affine3d &matchXform = mMatchTransformations[matchID];
		Eigen::Matrix3Xd matSourcePoints;
		if (!SampleUtil::buildMatrices(mSourceNodeSamples[sourceID].positions, matSourcePoints)) return false;
		Eigen::Matrix3Xd fittingSourcePoints = fittingXform * matSourcePoints;
		Eigen::Matrix3Xd xformedSourcePoints = matchXform * matSourcePoints;
		vector<TSlot> fittingSourceSlots = mSourceNodeSlots[sourceID];
		vector<TSlot> xformedSourceSlots = mSourceNodeSlots[sourceID];
		for (TSlot &slot : fittingSourceSlots) {
			slot.samples = fittingXform * slot.samples;
			slot.center = fittingXform * slot.center;
		}
		for (TSlot &slot : xformedSourceSlots) {
			slot.samples = matchXform * slot.samples;
			slot.center = matchXform * slot.center;
		}

		// match slots
		if (!TUtil::alignMatchedSlots(
			fittingSourcePoints,
			xformedSourcePoints,
			fittingSourceSlots,
			xformedSourceSlots,
			mTargetNodeSlots[targetID],
			mAlignedNodeSlots[targetID])) return false;
	}

	if (TUtil::gDebugVisualization) {
		if (!visualizeSlots(
			"alignedSlots.ply",
			mTargetWorkingNodes,
			mAlignedNodeSlots)) return false;
	}

	return true;
}

bool ContextPartGraphAssemble::reuseSlots() {

	// check if previous slots from source shape can be re-used (on previously replaced nodes)

	double reusingDistanceThreshold = StyleSynthesisConfig::mAssemble_SlotsReusingDistanceThreshold;

	vector<vec2i> mappingPairs; // (source node ID, target node ID) : # of valid mappings
	vector<Eigen::Affine3d> mappingXforms; // xform from source node to target node : # of valid mappings

	for (int matchID = 0; matchID < (int)mMatchPairs.size(); matchID++) {
		int sourceID = mMatchPairs[matchID][0];
		int targetID = mMatchPairs[matchID][1];
		Eigen::Affine3d &xform = mMatchTransformations[matchID];
		mappingPairs.push_back(vec2i(sourceID, targetID));
		mappingXforms.push_back(xform);
	}

	for (int targetID = 0; targetID < (int)mpTargetGraph->mAllNodes.size(); targetID++) {
		vec2i mapping = mpPreviousSolution->mReplaceMapping[targetID];
		if (mapping[1] >= 0) {
			int sourceID = mapping[0];
			Eigen::Affine3d &xform = mpPreviousSolution->mReplaceTransformation[targetID];
			mappingPairs.push_back(vec2i(sourceID, targetID));
			mappingXforms.push_back(xform);
		}
	}

	// update slots with source slots

	set<vec2i> reusedSlotsSet;
	vector<pair<vec2i, vec2i>> reusedSlotsPairs;

	for (int pairID = 0; pairID < (int)mappingPairs.size(); pairID++) {

		int sourceID = mappingPairs[pairID][0];
		int targetID = mappingPairs[pairID][1];
		Eigen::Affine3d xform = mappingXforms[pairID];

		map<int, int> sourceNodeMap; // source slot adjacent node ID => source slot ID
		for (TSlot &slot : mSourceNodeSlots[sourceID]) {
			if (slot.adjacentIndex[0] >= 0) {
				sourceNodeMap[slot.adjacentIndex[0]] = slot.selfIndex[1];
			}
		}

		for (TSlot &slot : mAlignedNodeSlots[targetID]) {
			if (slot.isVirtualSlot) continue; // skip virtual slot
			if (reusedSlotsSet.find(slot.selfIndex) != reusedSlotsSet.end()) continue; // already re-used
			int slotID = slot.selfIndex[1];
			int nbTargetID = slot.adjacentIndex[0];
			int nbSlotID = slot.adjacentIndex[1];
			if (nbTargetID < 0) continue;
			TSlot &nbSlot = mAlignedNodeSlots[nbTargetID][nbSlotID];
			if (nbSlot.isVirtualSlot) continue; // skip virtual slot
			if (reusedSlotsSet.find(nbSlot.selfIndex) != reusedSlotsSet.end()) continue; // already re-used

			int nbSourceID = mpPreviousSolution->mReplaceMapping[nbTargetID][0];
			auto &it = sourceNodeMap.find(nbSourceID);
			if (it != sourceNodeMap.end()) { // found slots connecting those two nodes on source shape
				int sourceSlotID = it->second;
				TSlot &sourceSlot = mSourceNodeSlots[sourceID][sourceSlotID];
				int nbSourceSlotID = sourceSlot.adjacentIndex[1];
				TSlot &nbSourceSlot = mSourceNodeSlots[nbSourceID][nbSourceSlotID];

				// re-use previous slots
				TSlot newSlot = slot;
				TSlot newNbSlot = nbSlot;
				Eigen::Affine3d nbXform = mpPreviousSolution->mReplaceTransformation[nbTargetID];
				newSlot.samples = xform * sourceSlot.samples;
				newSlot.center = newSlot.samples.rowwise().mean();
				newNbSlot.samples = nbXform * nbSourceSlot.samples;
				newNbSlot.center = newNbSlot.samples.rowwise().mean();
				newSlot.isVirtualSlot = false;
				newSlot.isReUsedSlot = true;
				newNbSlot.isVirtualSlot = false;
				newNbSlot.isReUsedSlot = true;
				if ((newSlot.center - newNbSlot.center).norm() < reusingDistanceThreshold) {
					slot = newSlot;
					nbSlot = newNbSlot;
					reusedSlotsSet.insert(slot.selfIndex);
					reusedSlotsSet.insert(nbSlot.selfIndex);
					reusedSlotsPairs.push_back(make_pair(slot.selfIndex, nbSlot.selfIndex));
				}
			}
		}
	}

	if (TUtil::gDebugVisualization) {
		if (!visualizeSlots(
			"reusedSlots.ply",
			mTargetWorkingNodes,
			mAlignedNodeSlots)) return false;
		system("pause");
	}
	/*
	// check connectivity from initial solution

	vector<vec2i> unprocessedNodePairs;
	if (true) {

		int numTargetNodes = (int)mpTargetGraph->mAllNodes.size();

		// mark current working target nodes
		vector<bool> workingNodeFlag(numTargetNodes, false);
		for (int nodeID : mTargetWorkingNodes) workingNodeFlag[nodeID] = true;

		// mark processed target node pairs
		set<vec2i> processedNodePairs;
		for (int nodeID = 0; nodeID < numTargetNodes; nodeID++) {
			for (TSlot &slot : mAlignedNodeSlots[nodeID]) {
				int neighborID = slot.adjacentIndex[0];
				if (neighborID < 0) continue;
				vec2i key(nodeID, neighborID);
				if (key[0] > key[1]) swap(key[0], key[1]);
				processedNodePairs.insert(key);
			}
		}

		// map target nodes
		map<int, int> nodeMapping; // initial target node ID => current target node ID
		for (int nodeID = 0; nodeID < numTargetNodes; nodeID++) {
			int initNodeID = mpPreviousSolution->mReplaceMapping[nodeID][1];
			if(initNodeID >= 0) nodeMapping[initNodeID] = nodeID;
		}
		for (int initNodeID = 0; initNodeID < (int)mpPreviousSolution->mNodeMapping.size(); initNodeID++) {
			int mappedNodeID = mpPreviousSolution->mNodeMapping[initNodeID];
			if (mappedNodeID >= 0) nodeMapping[initNodeID] = mappedNodeID;
		}

		// check "unprocessed" adjacent node pairs on initial target graph
		for (int nodeID = 0; nodeID < (int)mpInitialSolution->mGraph.mAllNodes.size(); nodeID++) {
			int mappedNodeID = nodeMapping[nodeID];
			if (mappedNodeID < 0 || !workingNodeFlag[mappedNodeID]) continue;
			for (TNode *nbNode : mpInitialSolution->mGraph.mAllNodes[nodeID]->mAdjacent) {
				int neighborID = nbNode->mUID;
				int mappedNeighborID = nodeMapping[neighborID];
				if (mappedNeighborID < 0 || !workingNodeFlag[mappedNeighborID]) continue;
				vec2i key(mappedNodeID, mappedNeighborID);
				if (key[0] > key[1]) swap(key[0], key[1]);
				if (processedNodePairs.find(key) != processedNodePairs.end()) continue;
				unprocessedNodePairs.push_back(key);
				processedNodePairs.insert(key);
			}
		}
	}

	// build virtual nodes for "unprocessed" adjacent node pairs from initial solution

	if (true) {
		for (vec2i nodePair : unprocessedNodePairs) {
			int sourceID = nodePair[0]; // one node in candidate graph
			int targetID = nodePair[1]; // another node in candidate graph
			Eigen::Matrix3Xd sourceSampleMat, targetSampleMat;
			if (!SampleUtil::buildMatrices(mTargetNodeSamples[sourceID], sourceSampleMat)) return false;
			if (!SampleUtil::buildMatrices(mTargetNodeSamples[targetID], targetSampleMat)) return false;
			SKDTree sourceTree, targetTree;
			SKDTreeData sourceTreeData, targetTreeData;
			if (!SampleUtil::buildKdTree(sourceSampleMat, sourceTree, sourceTreeData)) return false;
			if (!SampleUtil::buildKdTree(targetSampleMat, targetTree, targetTreeData)) return false;
			Eigen::VectorXi sourceIndices, targetIndices;
			if (!SampleUtil::findNearestNeighbors(sourceTree, targetSampleMat, targetIndices)) return false;
			if (!SampleUtil::findNearestNeighbors(targetTree, sourceSampleMat, sourceIndices)) return false;
			vector<vec2i> adjacentPointList;
			for (int pointID = 0; pointID < (int)sourceIndices.size(); pointID++) {
				int nearestID = sourceIndices[pointID];
				if (targetIndices[nearestID] == pointID) {
					adjacentPointList.push_back(vec2i(pointID, nearestID));
				}
			}
			if (!adjacentPointList.empty()) {
				// add slots
				int numSourceSlots = (int)mAlignedNodeSlots[sourceID].size();
				int numTargetSlots = (int)mAlignedNodeSlots[targetID].size();
				if ((int)mTargetNodeSlots[sourceID].size() != numSourceSlots ||
					(int)mTargetNodeSlots[targetID].size() != numTargetSlots)
				{
					cout << "Error: inconsistent number of slots between target graph and aligned graph" << endl;
					return false;
				}
				TSlot sourceSlot, targetSlot;
				sourceSlot.samples.resize(3, (int)adjacentPointList.size());
				targetSlot.samples.resize(3, (int)adjacentPointList.size());
				for (int id = 0; id < (int)adjacentPointList.size(); id++) {
					sourceSlot.samples.col(id) = sourceSampleMat.col(adjacentPointList[id][0]);
					targetSlot.samples.col(id) = targetSampleMat.col(adjacentPointList[id][1]);
				}
				sourceSlot.center = sourceSlot.samples.rowwise().mean();
				targetSlot.center = targetSlot.samples.rowwise().mean();
				sourceSlot.selfIndex = vec2i(sourceID, numSourceSlots);
				targetSlot.selfIndex = vec2i(targetID, numTargetSlots);
				sourceSlot.adjacentIndex = targetSlot.selfIndex;
				targetSlot.adjacentIndex = sourceSlot.selfIndex;
				sourceSlot.isVirtualSlot = true;
				sourceSlot.isReUsedSlot = false;
				targetSlot.isVirtualSlot = true;
				targetSlot.isReUsedSlot = false;

				mAlignedNodeSlots[sourceID].push_back(sourceSlot);
				mAlignedNodeSlots[targetID].push_back(targetSlot);
				mTargetNodeSlots[sourceID].push_back(sourceSlot);
				mTargetNodeSlots[targetID].push_back(targetSlot);
			}
		}
	}
	*/
	return true;
}

bool ContextPartGraphAssemble::alignSlots(double &outError) {

	int numMatchPairs = (int)mMatchPairs.size();
	int numSourceNodes = (int)mpSourceGraph->mAllNodes.size();
	int numTargetNodes = (int)mpTargetGraph->mAllNodes.size();

	// initialize node transformations

	mNodeTransformations.resize(numTargetNodes);
	for (int nodeID = 0; nodeID < numTargetNodes; nodeID++) {
		mNodeTransformations[nodeID].setIdentity();
	}

	outError = 0;
	if (numMatchPairs == 0) return true; // no matching
	double errorThreshold = StyleSynthesisConfig::mAssemble_SlotsAlignmentErrorThreshold;

	// determine global alignment groups & gather alignment info

	vector<vector<int>> alignmentGroups(0);
	vector<int> alignmentPrimitive(0);
	vector<Eigen::Vector3d> alignmentOrientation(0);
	if (true) {
		vector<int> replaceGroup(0);
		for (int pairID = 0; pairID < numMatchPairs; pairID++) {
			int targetID = mMatchPairs[pairID][1];
			replaceGroup.push_back(targetID);
		}
		vector<bool> visited(numTargetNodes, false);
		for (int pairID = 0; pairID < numMatchPairs; pairID++) {
			int sourceID = mMatchPairs[pairID][0];
			int targetID = mMatchPairs[pairID][1];
			visited[targetID] = true;

			TNode *srcNode = mpSourceGraph->mAllNodes[sourceID];
			int primitive = srcNode->mNodeDescriptors.mPrimitive;
			Eigen::Affine3d xform = mMatchTransformations[pairID];
			Eigen::Vector3d orientation = xform.rotation() * srcNode->mNodeDescriptors.mMajorOrientation;

			alignmentGroups.push_back(vector<int>(1, targetID));
			alignmentPrimitive.push_back(primitive);
			alignmentOrientation.push_back(orientation);
		}
		for (int id = 0; id < (int)mTargetWorkingNodes.size(); id++) {

			int nodeID = mTargetWorkingNodes[id];
			if (visited[nodeID]) continue;
			visited[nodeID] = true;

			vector<int> queue(1, nodeID);
			int head = 0;
			//head = 1; // HACK: skip extracting connected components
			while (head < (int)queue.size()) {
				int nodeID = queue[head];
				for (TSlot &slot : mAlignedNodeSlots[nodeID]) {
					int otherID = slot.adjacentIndex[0];
					int otherSlotID = slot.adjacentIndex[1];
					if (otherID < 0 || otherSlotID < 0 || visited[otherID]) continue;
					TSlot &otherSlot = mAlignedNodeSlots[otherID][otherSlotID];
					if (slot.isVirtualSlot || otherSlot.isVirtualSlot) continue;
					visited[otherID] = true;
					queue.push_back(otherID);
				}
				head++;
			}

			alignmentGroups.push_back(queue);

			// allow global alignment for all nodes
			//TNode *tgtNode = mpTargetGraph->mAllNodes[queue[0]];
			//alignmentPrimitive.push_back(tgtNode->mNodeDescriptors.mPrimitive);
			//alignmentOrientation.push_back(tgtNode->mNodeDescriptors.mMajorOrientation);

			// special handling for non-replaced nodes
			alignmentPrimitive.push_back(-1);
			alignmentOrientation.push_back(Eigen::Vector3d::UnitY());
		}
	}

	// compute baseline label error and initial error before global alignment

	Eigen::VectorXd baselineError;
	if (true) {
		if (!TUtil::checkSlotsAlignment(
			alignmentGroups,
			mTargetNodeSlots,
			baselineError)) return false;

		Eigen::VectorXd initialError;
		if (!TUtil::checkSlotsAlignment(
			alignmentGroups,
			mAlignedNodeSlots,
			initialError)) return false;
		Eigen::VectorXd realError = (initialError.cwiseSqrt() - baselineError.cwiseSqrt()).cwiseMax(0.0);
		outError = realError.norm() / cml::sqrt_safe((double)realError.size());
#ifdef OUTPUT_PROGRESS
		cout << "Slot error: " << outError << endl;
		if (TUtil::gDebugVisualization) {
			cout << "Baseline: " << baselineError.cwiseSqrt().transpose() << endl;
			cout << "Current: " << initialError.cwiseSqrt().transpose() << endl;
			cout << "Delta: " << realError.transpose() << endl;
			cout << endl;
		}
#endif

		//if (outError < errorThreshold) return true; // early quit
	}

	// iterative global slots alignment

	int numIterations = StyleSynthesisConfig::mAssemble_SlotsAlignmentIteration;
	//if (alignmentPrimitive[0] == 2) numIterations = 0; // HACK: ...
	//if (numIterations) outError = DBL_MAX; // always do global alignment

	vector<Eigen::Affine3d> currentTransformations = mNodeTransformations;
	if (TUtil::gDebugVisualization) {
		if (!visualizeNodes("before-global.ply")) return false;
		if (!visualizeSlots(
			"alignedSlots.ply",
			mTargetWorkingNodes,
			mAlignedNodeSlots)) return false;
		system("pause");
	}

	for (int iterID = 0; iterID < numIterations; iterID++) {
#ifdef OUTPUT_PROGRESS
		cout << "================ Iteration " << iterID << " ================" << endl;
#endif
		
		int globalMode = StyleSynthesisConfig::mAssemble_GlobalAlignmentMode;
		if (globalMode < 0) globalMode = mAlignmentMode;
		vector<Eigen::Affine3d> newTransformations(numTargetNodes);
		if (!TUtil::alignGlobalSlots(
			globalMode,
			alignmentGroups,
			alignmentPrimitive,
			alignmentOrientation,
			mAlignedNodeSlots,
			newTransformations)) return false;

		for (int nodeID : mTargetWorkingNodes) {
			currentTransformations[nodeID] = newTransformations[nodeID] * currentTransformations[nodeID];
			//cout << "Global alignment for node " << nodeID << endl;
			//cout << newTransformations[nodeID].matrix() << endl;
		}

		// check slots alignment after each iteration

		double currentError;
		if (true) {
			Eigen::VectorXd labelError;
			if (!TUtil::checkSlotsAlignment(
				alignmentGroups,
				mAlignedNodeSlots,
				labelError)) return false;
			Eigen::VectorXd realError = (labelError.cwiseSqrt() - baselineError.cwiseSqrt()).cwiseMax(0.0);
			currentError = realError.norm() / cml::sqrt_safe((double)realError.size());

#ifdef OUTPUT_PROGRESS
			cout << "Slot error: " << currentError << endl;
			if (TUtil::gDebugVisualization) {
				cout << "Baseline: " << baselineError.cwiseSqrt().transpose() << endl;
				cout << "Current: " << labelError.cwiseSqrt().transpose() << endl;
				cout << "Delta: " << realError.transpose() << endl;
				cout << endl;
			}
#endif
		}

		if (currentError < outError) {
			outError = currentError;
			mNodeTransformations = currentTransformations;
		}

		// update virtual slots

		if (true) {
			vector<Eigen::Matrix3Xd> xformedSourcePoints(numMatchPairs);
			for (int pairID = 0; pairID < numMatchPairs; pairID++) {
				int sourceID = mMatchPairs[pairID][0];
				int targetID = mMatchPairs[pairID][1];
				Eigen::Matrix3Xd sourceMat;
				if (!SampleUtil::buildMatrices(mSourceNodeSamples[sourceID], sourceMat)) return false;
				Eigen::Affine3d nodeXform = currentTransformations[targetID];
				Eigen::Affine3d matchXform = mMatchTransformations[pairID];
				xformedSourcePoints[pairID] = nodeXform * matchXform * sourceMat;
			}

			if (!TUtil::updateVirtualSlots(xformedSourcePoints, mMatchPairs, mAlignedNodeSlots)) return false;
		}

		if (TUtil::gDebugVisualization) {
			mNodeTransformations.swap(currentTransformations);
			if (!visualizeNodes("after-global.ply")) return false;
			mNodeTransformations.swap(currentTransformations);
			if (!visualizeSlots(
				"alignedSlots.ply",
				mTargetWorkingNodes,
				mAlignedNodeSlots)) return false;
			system("pause");
		}

		if (StyleSynthesisConfig::mAssemble_SlotsAlignmentAllowEarlyQuit) {
			if (outError < errorThreshold*0.5) return true; // early quit
		}
	}

	return true;
}

bool ContextPartGraphAssemble::postAlignSlots(double &outError) {

	int numSourceNodes = (int)mpSourceGraph->mAllNodes.size();
	int numTargetNodes = (int)mpTargetGraph->mAllNodes.size();

	// regularize replacing transformation

	double regEps = StyleSynthesisConfig::mAssemble_PostRegularizationEpsilon;

	mNodeTransformations.resize(numTargetNodes);
	for (int nodeID = 0; nodeID < numTargetNodes; nodeID++) {
		mNodeTransformations[nodeID].setIdentity();
	}

	vector<Eigen::Affine3d> regReplaceXforms = mpPreviousSolution->mReplaceTransformation;
	for (int id = 0; id < (int)mTargetWorkingNodes.size(); id++) {
		int targetID = mTargetWorkingNodes[id];
		int sourceID = mpPreviousSolution->mReplaceMapping[targetID][0];
		Eigen::Affine3d repXform = mpPreviousSolution->mReplaceTransformation[targetID];
		Eigen::Affine3d &regXform = regReplaceXforms[targetID];
		Eigen::Matrix3Xd nodeMat;
		if (!SampleUtil::buildMatrices(mSourceNodeSamples[sourceID], nodeMat)) return false;
		Eigen::Vector3d nodeCenter = nodeMat.rowwise().mean();
#ifdef OUTPUT_PROGRESS
		cout << "Regularizing node " << targetID << endl;
#endif
		if (!TUtil::regularizeNodeTransformation(repXform, regXform, nodeCenter, regEps)) return false;
		mNodeTransformations[targetID] = regXform * repXform.inverse();
	}

	// align exemplar slots

	if (!TUtil::alignExemplarSlots(
		mpPreviousSolution->mReplaceMapping,
		regReplaceXforms,
		mNodeTransformations,
		mSourceNodeSlots,
		mTargetNodeSlots,
		mAlignedNodeSlots)) return false;

	if (TUtil::gDebugVisualization) {
		if (!visualizeNodes("regularized.ply")) return false;
		if (!visualizeSlots(
			"regSlots.ply",
			mTargetWorkingNodes,
			mAlignedNodeSlots)) return false;
		system("pause");
	}

	outError = 0;
	double errorThreshold = StyleSynthesisConfig::mAssemble_SlotsAlignmentErrorThreshold;

	// compute initial error before global alignment

	vector<vector<int>> alignmentGroups(mTargetWorkingNodes.size());
	for (int id = 0; id < (int)mTargetWorkingNodes.size(); id++) {
		alignmentGroups[id].assign(1, mTargetWorkingNodes[id]);
	}

	if (true) {
		Eigen::VectorXd initialError;
		if (!TUtil::checkSlotsAlignment(
			alignmentGroups,
			mAlignedNodeSlots,
			initialError)) return false;
		outError = cml::sqrt_safe(initialError.mean());
#ifdef OUTPUT_PROGRESS
		cout << "Slot error: " << outError << endl;
		if (TUtil::gDebugVisualization) {
			cout << "Current: " << initialError.cwiseSqrt().transpose() << endl;
		}
#endif

		//if (outError < errorThreshold) return true; // early quit
	}

	// iterative global slots alignment

	int numIterations = StyleSynthesisConfig::mAssemble_SlotsAlignmentIteration;
	//if (numIterations) outError = DBL_MAX; // always do global alignment

	vector<Eigen::Affine3d> currentTransformations = mNodeTransformations;
	if (TUtil::gDebugVisualization) {
		if (!visualizeNodes("before-global.ply")) return false;
		if (!visualizeSlots(
			"alignedSlots.ply",
			mTargetWorkingNodes,
			mAlignedNodeSlots)) return false;
		system("pause");
	}

	// determine alignment order

	vector<int> orderedWorkingNodes;
	if (true) {
		vector<vector<TSlot>> &allSlots = mAlignedNodeSlots;
		vector<int> &workingNodes = mTargetWorkingNodes;
		int numAllNodes = (int)mAlignedNodeSlots.size();
		int numWorkingNodes = (int)workingNodes.size();
		// build node graph
		vector<set<int>> nodeGraph(numAllNodes);
		for (int nodeID : workingNodes) {
			for (TSlot &slot : allSlots[nodeID]) {
				int nbNodeID = slot.adjacentIndex[0];
				if (nbNodeID >= 0) {
					nodeGraph[nodeID].insert(nbNodeID);
					nodeGraph[nbNodeID].insert(nodeID);
				}
			}
		}
		// compute all pair distance (Floyd algorithm)
		Eigen::MatrixXi distMat(numAllNodes, numAllNodes);
		distMat.setZero();
		for (int i : workingNodes) for (int j : workingNodes) distMat(i, j) = numWorkingNodes; // use this as max distance
		for (int k : workingNodes) distMat(k, k) = 0;
		for (int i : workingNodes) for (int j : nodeGraph[i]) distMat(i, j) = 1;
		for (int k : workingNodes) {
			for (int i : workingNodes) {
				if (k == i) continue;
				for (int j : workingNodes) {
					if (k == j || i == j) continue;
					int dist = distMat(i, k) + distMat(k, j);
					if (distMat(i, j) > dist) distMat(i, j) = dist;
				}
			}
		}
		// determine order by node radius
		Eigen::VectorXi radiusVec = distMat.rowwise().sum();
		vector<int> nodeScore(0);
		for (int nodeID : workingNodes) nodeScore.push_back(radiusVec[nodeID]);
		vector<int> orderList(numWorkingNodes);
		for (int k = 0; k < numWorkingNodes; k++) orderList[k] = k;
		sort(orderList.begin(), orderList.end(),
			[&nodeScore](int &lhs, int &rhs) { return nodeScore[lhs] < nodeScore[rhs]; });
		orderedWorkingNodes.resize(numWorkingNodes);
		for (int k = 0; k < numWorkingNodes; k++) orderedWorkingNodes[k] = workingNodes[orderList[k]];
	}

	for (int iterID = 0; iterID < numIterations; iterID++) {
#ifdef OUTPUT_PROGRESS
		cout << "================ Iteration " << iterID << " ================" << endl;
#endif

		vector<Eigen::Affine3d> newTransformations(numTargetNodes);
		if (!TUtil::alignGlobalSlots(
			orderedWorkingNodes,
			mAlignedNodeSlots,
			newTransformations)) return false;

		for (int nodeID : mTargetWorkingNodes) {
			currentTransformations[nodeID] = newTransformations[nodeID] * currentTransformations[nodeID];
			//cout << "Global alignment for node " << nodeID << endl;
			//cout << newTransformations[nodeID].matrix() << endl;
		}

		// check slots alignment after each iteration

		double currentError;
		if (true) {
			Eigen::VectorXd labelError;
			if (!TUtil::checkSlotsAlignment(
				alignmentGroups,
				mAlignedNodeSlots,
				labelError)) return false;
			currentError = cml::sqrt_safe(labelError.mean());

#ifdef OUTPUT_PROGRESS
			cout << "Slot error: " << currentError << endl;
			if (TUtil::gDebugVisualization) {
				cout << "Current: " << labelError.cwiseSqrt().transpose() << endl;
			}
#endif
		}

		if (currentError < outError) {
			outError = currentError;
			mNodeTransformations = currentTransformations;
		}

		if (TUtil::gDebugVisualization) {
			mNodeTransformations.swap(currentTransformations);
			if (!visualizeNodes("after-global.ply")) return false;
			mNodeTransformations.swap(currentTransformations);
			if (!visualizeSlots(
				"alignedSlots.ply",
				mTargetWorkingNodes,
				mAlignedNodeSlots)) return false;
			system("pause");
		}

		if (StyleSynthesisConfig::mAssemble_SlotsAlignmentAllowEarlyQuit) {
			if (outError < errorThreshold*0.5) return true; // early quit
		}
	}

	return true;
}

bool ContextPartGraphAssemble::alignAddings(double &outError) {

	int numTargetNodes = (int)mpTargetGraph->mAllNodes.size();

	mAddingNodeIndices.clear();
	mAddingTransformations.clear();

	outError = 0;
	int numTotalSlots = 0;

	for (int nodeID : mAddingNodes) {
		TNode *node = mpSourceGraph->mAllNodes[nodeID];

		int primitive = node->mNodeDescriptors.mPrimitive;
		Eigen::Vector3d orientation = node->mNodeDescriptors.mMajorOrientation;

		int numSlots = (int)mAddingNodeSlots[nodeID].size();
		if (numSlots <= 0) continue; // no reliable adjacent part to determine transfer destination
		numTotalSlots += numSlots;

		// get slot transforms list

		vector<Eigen::Affine3d> slotXformList(0);
		if (true) {
			vector<int> slotReplacedCount(numSlots, 0);
			for (int slotID = 0; slotID < numSlots; slotID++) {
				TSlot &slot = mAddingNodeSlots[nodeID][slotID];
				int adjNodeID = slot.adjacentIndex[0];
				for (int tgtNodeID = 0; tgtNodeID < numTargetNodes; tgtNodeID++) {
					if (mpPreviousSolution->mReplaceMapping[tgtNodeID][0] == adjNodeID &&
						mpPreviousSolution->mReplaceMapping[tgtNodeID][1] >= 0)
					{
						slotReplacedCount[slotID]++;
					}
				}
			}
			int minCount = slotReplacedCount[0];
			int minSlotID = 0;
			for (int slotID = 0; slotID < numSlots; slotID++) {
				int count = slotReplacedCount[slotID];
				if (count && count < minCount) {
					minCount = count;
					minSlotID = slotID;
				}
			}

			TSlot &slot = mAddingNodeSlots[nodeID][minSlotID];
			int adjNodeID = slot.adjacentIndex[0];
			for (int tgtNodeID = 0; tgtNodeID < numTargetNodes; tgtNodeID++) {
				if (mpPreviousSolution->mReplaceMapping[tgtNodeID][0] == adjNodeID &&
					mpPreviousSolution->mReplaceMapping[tgtNodeID][1] >= 0)
				{
					Eigen::Affine3d replaceXform = mpPreviousSolution->mReplaceTransformation[tgtNodeID];
					replaceXform = mNodeTransformations[tgtNodeID] * replaceXform;
					Eigen::Affine3d slotXform;
					slotXform.setIdentity();
					slotXform.pretranslate(-slot.center);
					slotXform.prerotate(replaceXform.rotation());
					slotXform.pretranslate(replaceXform * slot.center);
					slotXformList.push_back(slotXform);
				}
			}
		}

		// process each transformation in the list

		for (Eigen::Affine3d slotXform : slotXformList) {

			vector<TSlot> xformedSourceSlots;
			vector<TSlot> xformedTargetSlots;

			for (TSlot &slot : mAddingNodeSlots[nodeID]) {
				int adjNodeID = slot.adjacentIndex[0];
				int adjSlotID = slot.adjacentIndex[1];
				if (adjNodeID < 0) continue;
				TSlot &adjSlot = mAddingNodeSlots[adjNodeID][adjSlotID];

				Eigen::Affine3d replaceXform;
				double minDist = DBL_MAX;
				for (int tgtNodeID = 0; tgtNodeID < numTargetNodes; tgtNodeID++) {
					if (mpPreviousSolution->mReplaceMapping[tgtNodeID][0] == adjNodeID &&
						mpPreviousSolution->mReplaceMapping[tgtNodeID][1] >= 0)
					{
						Eigen::Affine3d currentXform = mpPreviousSolution->mReplaceTransformation[tgtNodeID];
						currentXform = mNodeTransformations[tgtNodeID] * currentXform;

						double dist = (currentXform * adjSlot.center - slotXform * slot.center).norm();
						if (dist < minDist) {
							minDist = dist;
							replaceXform = currentXform;
						}
					}
				}

				TSlot xformedSlot = slot;
				xformedSlot.center = slotXform * xformedSlot.center;
				xformedSlot.samples = slotXform * xformedSlot.samples;

				TSlot xformedAdjSlot = adjSlot;
				xformedAdjSlot.center = replaceXform * xformedAdjSlot.center;
				xformedAdjSlot.samples = replaceXform * xformedAdjSlot.samples;

				xformedSourceSlots.push_back(xformedSlot);
				xformedTargetSlots.push_back(xformedAdjSlot);
			}

			orientation = slotXform.rotation() * orientation;

			int globalMode = StyleSynthesisConfig::mAssemble_GlobalAlignmentMode;
			if (globalMode < 0) globalMode = mAlignmentMode;

			Eigen::Affine3d alignXform;
			double alignError;
			if (!TUtil::alignSlotPairs(
				globalMode,
				primitive, orientation,
				xformedSourceSlots, xformedTargetSlots,
				alignXform, alignError)) return false;

			Eigen::Affine3d addingXform = alignXform * slotXform;

			mAddingNodeIndices.push_back(nodeID);
			mAddingTransformations.push_back(addingXform);

			outError = max(outError, alignError);
		}
	}

	//if (numTotalSlots == 0) outError = DBL_MAX;

	if (TUtil::gDebugVisualization) {
		if (!visualizeNodes("after-adding.ply")) return false;
	}

	return true;
}

bool ContextPartGraphAssemble::visualize(string folderName) {

	string matchVName = folderName + "matching.ply";

	if (!visualizeMatchings(matchVName)) return false;

	return true;
}

bool ContextPartGraphAssemble::visualizeMatchings(string fileName) {
	
	TTriangleMesh *srcMesh = mpSourceGraph->mRootNode->mpGraphMesh;
	TTriangleMesh *tgtMesh = mpTargetGraph->mRootNode->mpGraphMesh;

	vector<vec3i> srcColor(srcMesh->indices.size(), vec3i(127, 127, 127));
	vector<vec3i> tgtColor(tgtMesh->indices.size(), vec3i(127, 127, 127));

	int numMatchPairs = (int)mMatchPairs.size();
	for (int matchID = 0; matchID < numMatchPairs; matchID++) {

		int sourceID = mMatchPairs[matchID][0];
		int targetID = mMatchPairs[matchID][1];
		TNode *srcNode = mpSourceGraph->mAllNodes[sourceID];
		TNode *tgtNode = mpTargetGraph->mAllNodes[targetID];
		vector<int> &srcSegment = (*mpSourceGraph->mRootNode->mpGraphSegments)[srcNode->mPartLevelID][srcNode->mPartSegmentID];
		vector<int> &tgtSegment = (*mpTargetGraph->mRootNode->mpGraphSegments)[tgtNode->mPartLevelID][tgtNode->mPartSegmentID];

		for (int faceID : srcSegment) srcColor[faceID] = vec3i(255, 0, 0);
		for (int faceID : tgtSegment) tgtColor[faceID] = vec3i(0, 0, 255);
	}

	vec3 srcBBMin, srcBBMax;
	vec3 tgtBBMin, tgtBBMax;
	if (!MeshUtil::computeAABB(*srcMesh, srcBBMin, srcBBMax)) return false;
	if (!MeshUtil::computeAABB(*tgtMesh, tgtBBMin, tgtBBMax)) return false;
	vec3 offset(0.0f, -(tgtBBMax[1] - srcBBMin[1]) * 1.2f, 0.0f);

	PlyExporter pe;
	if (!pe.addMesh(&srcMesh->positions, &srcMesh->normals, &srcMesh->indices, &srcColor)) return false;
	if (!pe.addMesh(&tgtMesh->positions, &tgtMesh->normals, &tgtMesh->indices, &tgtColor, offset)) return false;
	if (!pe.output(fileName)) return false;

	return true;
}

bool ContextPartGraphAssemble::visualizeNodes(string fileName) {

	int numTargetNodes = (int)mpTargetGraph->mAllNodes.size();
	vector<bool> replacedFlags(numTargetNodes, false);

	PlyExporter pe;

	int numMatchPairs = (int)mMatchPairs.size();
	for (int pairID = 0; pairID < numMatchPairs; pairID++) {
		int sourceID = mMatchPairs[pairID][0];
		int targetID = mMatchPairs[pairID][1];
		replacedFlags[targetID] = true;
		Eigen::Affine3d xform = mNodeTransformations[targetID] * mMatchTransformations[pairID];
		Eigen::Matrix4d eigenMat = xform.matrix();
		matrix4f cmlMat;
		for (int col = 0; col < 4; col++) {
			for (int row = 0; row < 4; row++) {
				cmlMat(row, col) = (float)eigenMat(row, col);
			}
		}
		TSampleSet &points = mSourceNodeSamples[sourceID];
		if (!pe.addPoint(&points.positions, &points.normals, cmlMat, SegmentUtil::colorMapping(pairID))) return false;
	}

	for (int nodeID : mRemovingNodes) {
		replacedFlags[nodeID] = true;
	}

	for (int nodeID : mTargetWorkingNodes) {
		if (replacedFlags[nodeID]) continue;
		Eigen::Affine3d xform = mNodeTransformations[nodeID];
		Eigen::Matrix4d eigenMat = xform.matrix();
		matrix4f cmlMat;
		for (int col = 0; col < 4; col++) {
			for (int row = 0; row < 4; row++) {
				cmlMat(row, col) = (float)eigenMat(row, col);
			}
		}
		TSampleSet &points = mTargetNodeSamples[nodeID];
		if (!pe.addPoint(&points.positions, &points.normals, cmlMat, vec3i(127, 127, 127))) return false;
	}

	for (int addingID = 0; addingID < (int)mAddingNodeIndices.size(); addingID++) {
		int nodeID = mAddingNodeIndices[addingID];
		Eigen::Affine3d xform = mAddingTransformations[addingID];
		Eigen::Matrix4d eigenMat = xform.matrix();
		matrix4f cmlMat;
		for (int col = 0; col < 4; col++) {
			for (int row = 0; row < 4; row++) {
				cmlMat(row, col) = (float)eigenMat(row, col);
			}
		}
		TSampleSet &points = mSourceNodeSamples[nodeID];
		if (!pe.addPoint(&points.positions, &points.normals, cmlMat, SegmentUtil::colorMapping(addingID))) return false;
	}

	if (!pe.output(fileName)) return false;

	return true;
}

bool ContextPartGraphAssemble::visualizeSlots(
	string fileName,
	vector<int> &workingNodes,
	vector<vector<TSlot>> &nodeSlots)
{

	PlyExporter pe;
	for (int nodeID : workingNodes) {
		int numSlots = (int)nodeSlots[nodeID].size();
		vec3i color = SegmentUtil::colorMapping(nodeID);
		for (int slotID = 0; slotID < numSlots; slotID++) {
			TSlot &slot = nodeSlots[nodeID][slotID];
			vector<vec3> slotPoints;
			for (int id = 0; id < (int)slot.samples.cols(); id++) {
				slotPoints.push_back(vec3d(slot.samples(0, id), slot.samples(1, id), slot.samples(2, id)));
			}
			vec3i slotColor = slot.isVirtualSlot ? vec3i(255, 255, 255) : color;
			if (!pe.addPoint(&slotPoints, 0, cml::identity_4x4(), slotColor)) return false;
			color = color / 2;
		}
	}
	if (!pe.output(fileName)) return false;
	
	return true;
}
