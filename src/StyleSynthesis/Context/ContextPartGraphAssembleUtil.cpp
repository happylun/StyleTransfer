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

#include "ContextPartGraphAssembleUtil.h"

#include <iostream>
#include <fstream>
#include <set>
#include <unordered_set>

#include "ContextPartGraph.h"

#include "Mesh/MeshUtil.h"

#include "Sample/SampleSimplePoissonDisk.h"
#include "Sample/SampleUtil.h"

#include "Match/MatchLabeledICP.h"
#include "Match/MatchSimpleICP.h"
#include "Match/MatchPrimitiveICP.h"

#include "Utility/PlyExporter.h"

#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

#define OUTPUT_PROGRESS

bool ContextPartGraphAssembleUtil::gDebugVisualization = false;

bool ContextPartGraphAssembleUtil::markNodeTypes(
	TGraph *inGraph,
	vector<int> &inReplaceNodeList,
	vector<int> &inRemoveNodeList,
	vector<int> &inPreviousNodeList,
	vector<TNodeType> &outTypes)
{

	int numNodes = (int)inGraph->mAllNodes.size();
	outTypes.assign(numNodes, TNodeType_ACTIVE);

	// process replacing nodes

	for (int nodeID : inReplaceNodeList) {

		// mark this node as "replaced"
		outTypes[nodeID] = TNodeType_REPLACED;

		// mark all ancester nodes as "modified"
		int parentID = inGraph->mAllNodes[nodeID]->mParent->mUID;
		while (parentID >= 0 && outTypes[parentID] == TNodeType_ACTIVE) {
			outTypes[parentID] = TNodeType_MODIFIED;
			parentID = inGraph->mAllNodes[parentID]->mParent->mUID;
		}

		// mark all descendant nodes as "removed"
		vector<int> queue(1, nodeID);
		int head = 0;
		while (head < (int)queue.size()) {
			int currentID = queue[head];
			for (TNode *child : inGraph->mAllNodes[currentID]->mChildren) {
				int childID = child->mUID;
				if (outTypes[childID] == TNodeType_ACTIVE) {
					queue.push_back(childID);
					outTypes[childID] = TNodeType_REMOVED;
				}
			}
			head++;
		}
	}

	// process removing nodes

	for (int nodeID : inRemoveNodeList) {

		// mark this node as "removed"
		outTypes[nodeID] = TNodeType_REMOVED;

		// mark all ancester nodes as "modified"
		int parentID = inGraph->mAllNodes[nodeID]->mParent->mUID;
		while (parentID >= 0) {
			outTypes[parentID] = TNodeType_MODIFIED;
			parentID = inGraph->mAllNodes[parentID]->mParent->mUID;
		}

		// mark all descendant nodes as "removed"
		vector<int> queue(1, nodeID);
		int head = 0;
		while (head < (int)queue.size()) {
			int currentID = queue[head];
			for (TNode *child : inGraph->mAllNodes[currentID]->mChildren) {
				int childID = child->mUID;
				queue.push_back(childID);
				outTypes[childID] = TNodeType_REMOVED;
			}
			head++;
		}
	}

	// process previously removed nodes

	for (int nodeID : inPreviousNodeList) {

		// implicitly this node is marked as "active" and all its descendants will be marked as "inactive"

		// mark all ancester nodes as "modified"
		int parentID = inGraph->mAllNodes[nodeID]->mParent->mUID;
		while (parentID >= 0) {
			outTypes[parentID] = TNodeType_MODIFIED;
			parentID = inGraph->mAllNodes[parentID]->mParent->mUID;
		}
	}

	// mark active nodes as "inactive" if its parent is active and is not the root

	for (int nodeID = numNodes - 1; nodeID >= 0; nodeID--) {
		int parentID = inGraph->mAllNodes[nodeID]->mParent->mUID;
		if (parentID >= 0 && outTypes[parentID] == TNodeType_ACTIVE) {
			outTypes[nodeID] = TNodeType_INACTIVE;
		}
	}

	return true;
}

bool ContextPartGraphAssembleUtil::extractNodeMesh(TNode *inNode, TTriangleMesh &outMesh) {

	if (outMesh.amount) return true; // already extracted

	vector<int> &nodeSegment = (*(inNode->mpGraphSegments))[inNode->mPartLevelID][inNode->mPartSegmentID];
	if (!MeshUtil::extractSubMesh(*inNode->mpGraphMesh, nodeSegment, outMesh)) return false;

	return true;
}

bool ContextPartGraphAssembleUtil::extractMeshSamples(TTriangleMesh &inMesh, TSampleSet &outSamples) {

	if (outSamples.amount) return true; // already extracted

	if (inMesh.indices.empty()) {
		// empty mesh
		outSamples.positions.clear();
		outSamples.normals.clear();
		outSamples.indices.clear();
		outSamples.amount = 0;
		return true;
	}

	double sampleRadius = StyleSynthesisConfig::mAssemble_SlotsSamplingRadius;

	SampleSimplePoissonDisk sspd(&inMesh);
	if (!sspd.runSampling(sampleRadius)) return false;
	if (!sspd.exportSample(outSamples)) return false;

	if (outSamples.amount < 10) { // redo sampling if not enough sample points
		if (!sspd.runSampling(20)) return false;
		if (!sspd.exportSample(outSamples)) return false;
	}

	return true;
}

bool ContextPartGraphAssembleUtil::extractNodeSlots(
	TGraph *inGraph,
	vector<TSampleSet> &inSamples,
	vector<int> &inWorkingNodes,
	vector<vector<TSlot>> &outSlots)
{

	// thresholds

	double distEps = StyleSynthesisConfig::mAssemble_SlotsSamplingRadius;
	double neighborDistance = distEps * 1.0; // UNDONE: param neighbor distance for determining slot points
	double wallDistance = distEps * 2.0;
	double virtualDistance = distEps * 2.0;

	// build KD trees for working nodes

	vector<SKDTree*> trees(inGraph->mAllNodes.size(), 0);
	vector<SKDTreeData> treesData(inGraph->mAllNodes.size());
	for (int nodeID : inWorkingNodes) {
		trees[nodeID] = new SKDTree();
		if (!SampleUtil::buildKdTree(inSamples[nodeID].positions, *(trees[nodeID]), treesData[nodeID])) return false;
	}

	// find adjacent node pairs

	vector<vec2i> adjacentNodePairs;
	if (true) {
		set<vec2i> nodePairSet;
		set<int> workingNodeSet(inWorkingNodes.begin(), inWorkingNodes.end());
		for (int nodeID : inWorkingNodes) {
			for (TNode *nbNode : inGraph->mAllNodes[nodeID]->mAdjacent) {
				int nbID = nbNode->mUID;
				if (workingNodeSet.find(nbID) != workingNodeSet.end()) {
					vec2i key(nodeID, nbID);
					if (key[0] > key[1]) swap(key[0], key[1]);
					nodePairSet.insert(key);
				}
			}
		}
		adjacentNodePairs.assign(nodePairSet.begin(), nodePairSet.end());
	}
	int numNodePairs = (int)adjacentNodePairs.size();

	// process each pair of adjacent nodes

	outSlots.assign(inGraph->mAllNodes.size(), vector<TSlot>(0));

	for (int pairID = 0; pairID < numNodePairs; pairID++) {
		int sourceID = adjacentNodePairs[pairID][0];
		int targetID = adjacentNodePairs[pairID][1];

		TSampleSet &srcSamples = inSamples[sourceID];
		TSampleSet &tgtSamples = inSamples[targetID];
		SKDTree &srcTree = *(trees[sourceID]);
		SKDTree &tgtTree = *(trees[targetID]);

		// find slot points

		vector<double> srcPointFlags(srcSamples.amount, -1);
		vector<double> tgtPointFlags(tgtSamples.amount, -1);
#pragma omp parallel for
		for (int srcPointID = 0; srcPointID < srcSamples.amount; srcPointID++) {
			vec3 point = srcSamples.positions[srcPointID];
			SKDT::NamedPoint queryPoint(point[0], point[1], point[2]);
			Thea::BoundedSortedArray<SKDTree::Neighbor> queryResult(20);
			tgtTree.kClosestElements<Thea::MetricL2>(queryPoint, queryResult, neighborDistance);
			for (int qID = 0; qID < queryResult.size(); qID++) {
				int tgtPointID = (int)tgtTree.getElements()[queryResult[qID].getIndex()].id;
				double dist = (tgtSamples.positions[tgtPointID] - point).length();
				if (tgtPointFlags[tgtPointID] < 0) {
					tgtPointFlags[tgtPointID] = dist;
				} else {
					tgtPointFlags[tgtPointID] = min(dist, tgtPointFlags[tgtPointID]);
				}
			}
		}
#pragma omp parallel for
		for (int tgtPointID = 0; tgtPointID < tgtSamples.amount; tgtPointID++) {
			vec3 point = tgtSamples.positions[tgtPointID];
			SKDT::NamedPoint queryPoint(point[0], point[1], point[2]);
			Thea::BoundedSortedArray<SKDTree::Neighbor> queryResult(20);
			srcTree.kClosestElements<Thea::MetricL2>(queryPoint, queryResult, neighborDistance);
			for (int qID = 0; qID < queryResult.size(); qID++) {
				int srcPointID = (int)srcTree.getElements()[queryResult[qID].getIndex()].id;
				double dist = (srcSamples.positions[srcPointID] - point).length();
				if (srcPointFlags[srcPointID] < 0) {
					srcPointFlags[srcPointID] = dist;
				} else {
					srcPointFlags[srcPointID] = min(dist, srcPointFlags[srcPointID]);
				}
			}
		}

		// build point set

		TPointSet srcPoints;
		TPointSet tgtPoints;
		double srcAvgDist = 0;
		double tgtAvgDist = 0;
		if (true) {
			srcPoints.positions.clear();
			srcPoints.normals.clear();
			for (int pointID = 0; pointID < srcSamples.amount; pointID++) {
				if (srcPointFlags[pointID] >= 0) {
					srcPoints.positions.push_back(srcSamples.positions[pointID]);
					srcPoints.normals.push_back(srcSamples.normals[pointID]);
					srcAvgDist += srcPointFlags[pointID];
				}
			}
			srcPoints.amount = (int)srcPoints.positions.size();
			if (srcPoints.amount) srcAvgDist /= srcPoints.amount;
		}
		if (true) {
			tgtPoints.positions.clear();
			tgtPoints.normals.clear();
			for (int pointID = 0; pointID < tgtSamples.amount; pointID++) {
				if (tgtPointFlags[pointID] >= 0) {
					tgtPoints.positions.push_back(tgtSamples.positions[pointID]);
					tgtPoints.normals.push_back(tgtSamples.normals[pointID]);
					tgtAvgDist += tgtPointFlags[pointID];
				}
			}
			tgtPoints.amount = (int)tgtPoints.positions.size();
			if (tgtPoints.amount) tgtAvgDist /= tgtPoints.amount;
		}
		if (srcPoints.amount == 0 || tgtPoints.amount == 0) continue; // no slot points extracted

		// output slots
		int srcSlotID = (int)outSlots[sourceID].size();
		int tgtSlotID = (int)outSlots[targetID].size();
		if (true) {
			outSlots[sourceID].push_back(TSlot());
			TSlot &srcSlot = outSlots[sourceID].back();
			if (!SampleUtil::buildMatrices(srcPoints, srcSlot.samples)) return false;
			srcSlot.center = srcSlot.samples.rowwise().mean();
			srcSlot.selfIndex = vec2i(sourceID, srcSlotID);
			srcSlot.adjacentIndex = vec2i(targetID, tgtSlotID);
			srcSlot.isVirtualSlot = (srcAvgDist > virtualDistance);
			srcSlot.isReUsedSlot = false;
		}
		if (true) {
			outSlots[targetID].push_back(TSlot());
			TSlot &tgtSlot = outSlots[targetID].back();
			if (!SampleUtil::buildMatrices(tgtPoints, tgtSlot.samples)) return false;
			tgtSlot.center = tgtSlot.samples.rowwise().mean();
			tgtSlot.selfIndex = vec2i(targetID, tgtSlotID);
			tgtSlot.adjacentIndex = vec2i(sourceID, srcSlotID);
			tgtSlot.isVirtualSlot = (tgtAvgDist > virtualDistance);
			tgtSlot.isReUsedSlot = false;
		}
	}

	// add floor/ceiling slots

	double floorHeight, ceilingHeight;
	if (true) {
		vec3 bbMin, bbMax;
		if (!MeshUtil::computeAABB(*inGraph->mRootNode->mpGraphMesh, bbMin, bbMax)) return false;
		floorHeight = bbMin[1];
		ceilingHeight = bbMax[1];
	}

	bool addFloor = StyleSynthesisConfig::mAssemble_AddFloorSlot;
	bool addCeiling = StyleSynthesisConfig::mAssemble_AddCeilingSlot;

	for (int nodeID : inWorkingNodes) {

		vector<vec3> floorPoints;
		vector<vec3> ceilingPoints;
		for (vec3 &point : inSamples[nodeID].positions) {
			if (fabs(point[1] - floorHeight) < wallDistance) {
				floorPoints.push_back(point);
			}
			if (fabs(point[1] - ceilingHeight) < wallDistance) {
				ceilingPoints.push_back(point);
			}
		}
		if (addFloor && !floorPoints.empty()) {
			// add floor slot
			int slotID = (int)outSlots[nodeID].size();
			outSlots[nodeID].push_back(TSlot());
			TSlot &floorSlot = outSlots[nodeID].back();
			if (!SampleUtil::buildMatrices(floorPoints, floorSlot.samples)) return false;
			floorSlot.center = floorSlot.samples.rowwise().mean();
			floorSlot.selfIndex = vec2i(nodeID, slotID);
			floorSlot.adjacentIndex = vec2i(-1, -1);
			floorSlot.isVirtualSlot = false;
			floorSlot.isReUsedSlot = false;
		}
		if (addCeiling && !ceilingPoints.empty()) {
			// add ceiling slot
			int slotID = (int)outSlots[nodeID].size();
			outSlots[nodeID].push_back(TSlot());
			TSlot &ceilingSlot = outSlots[nodeID].back();
			if (!SampleUtil::buildMatrices(ceilingPoints, ceilingSlot.samples)) return false;
			ceilingSlot.center = ceilingSlot.samples.rowwise().mean();
			ceilingSlot.selfIndex = vec2i(nodeID, slotID);
			ceilingSlot.adjacentIndex = vec2i(-2, -2);
			ceilingSlot.isVirtualSlot = false;
			ceilingSlot.isReUsedSlot = false;
		}
	}

	// don't forget to release memory from trees

	for (SKDTree* &tree : trees) {
		if (tree) delete tree;
		tree = 0;
	}

	return true;
}

bool ContextPartGraphAssembleUtil::computeSlotDistance(
	TSlot &sourceSlot,
	TSlot &targetSlot,
	double &nearestDistance,
	double &farthestDistance)
{

	//distance = (sourceSlot.center - targetSlot.center).norm();
	//return true;

	SKDTree srcTree, tgtTree;
	SKDTreeData srcTreeData, tgtTreeData;
	if (!SampleUtil::buildKdTree(sourceSlot.samples, srcTree, srcTreeData)) return false;
	if (!SampleUtil::buildKdTree(targetSlot.samples, tgtTree, tgtTreeData)) return false;

	Eigen::VectorXi srcSliceIndices, tgtSliceIndices;
	if (!SampleUtil::findNearestNeighbors(srcTree, targetSlot.samples, srcSliceIndices)) return false;
	if (!SampleUtil::findNearestNeighbors(tgtTree, sourceSlot.samples, tgtSliceIndices)) return false;

	Eigen::Matrix3Xd slicedSourcePoints, slicedTargetPoints;
	if (!SampleUtil::sliceMatrices(sourceSlot.samples, srcSliceIndices, slicedSourcePoints)) return false;
	if (!SampleUtil::sliceMatrices(targetSlot.samples, tgtSliceIndices, slicedTargetPoints)) return false;

	Eigen::VectorXd srcDist = (sourceSlot.samples - slicedTargetPoints).colwise().norm();
	Eigen::VectorXd tgtDist = (targetSlot.samples - slicedSourcePoints).colwise().norm();

	nearestDistance = min(srcDist.minCoeff(), tgtDist.minCoeff());
	farthestDistance = max(srcDist.maxCoeff(), tgtDist.maxCoeff());

	return true;
}

bool ContextPartGraphAssembleUtil::alignMatchedParts(
	int alignmentMode,
	TNode *sourceNode,
	TNode *targetNode,
	TSampleSet &sourceSamples,
	TSampleSet &targetSamples,
	Eigen::Affine3d &matchTransform,
	Eigen::Affine3d &fittingTransform)
{

	// get alignment info

	int srcPrimitive = sourceNode->mNodeDescriptors.mPrimitive;
	int tgtPrimitive = targetNode->mNodeDescriptors.mPrimitive;
	Eigen::Vector3d srcOrientation = sourceNode->mNodeDescriptors.mMajorOrientation;
	Eigen::Vector3d tgtOrientation = targetNode->mNodeDescriptors.mMajorOrientation;
	if (StyleSynthesisConfig::mContext_HandleRotationalSymmetry) {
		Eigen::Vector3d srcMassOffset = sourceNode->mNodeDescriptors.mBoundingBox.center(); // NOTE: don't use mass center -- may be flipped
		Eigen::Vector3d tgtMassOffset = targetNode->mNodeDescriptors.mBoundingBox.center();
		// special handling for rotational symmetric parts...
		if (srcPrimitive == 1 && srcMassOffset.squaredNorm()) {// && fabs(srcOrientation.dot(srcMassOffset.normalized())) < 0.1) {
			srcMassOffset = Eigen::Vector3d::UnitY().cross(srcMassOffset);
		}
		if (tgtPrimitive == 1 && tgtMassOffset.squaredNorm()) {// && fabs(tgtOrientation.dot(tgtMassOffset.normalized())) < 0.1) {
			tgtMassOffset = Eigen::Vector3d::UnitY().cross(tgtMassOffset);
		}
		if (srcOrientation.dot(srcMassOffset) < 0) srcOrientation = -srcOrientation;
		if (tgtOrientation.dot(tgtMassOffset) < 0) tgtOrientation = -tgtOrientation;
	}
	if (!MatchPrimitiveICP::regularizeAxis(srcOrientation)) return false;
	if (!MatchPrimitiveICP::regularizeAxis(tgtOrientation)) return false;
#ifdef OUTPUT_PROGRESS
	if (gDebugVisualization) {
		cout << "Source primitive: " << srcPrimitive << endl;
		cout << "Target primitive: " << tgtPrimitive << endl;
		cout << "Source orientation: " << srcOrientation.transpose() << endl;
		cout << "Target orientation: " << tgtOrientation.transpose() << endl;
	}
#endif

	// build matrices

	Eigen::Matrix3Xd sourceMat, targetMat;
	if (!SampleUtil::buildMatrices(sourceSamples, sourceMat)) return false;
	if (!SampleUtil::buildMatrices(targetSamples, targetMat)) return false;

	// pre-alignment

	Eigen::Affine3d preTransform = Eigen::Affine3d::Identity();
	if (true) {
		double angleThreshold = StyleSynthesisConfig::mAssemble_PrimitiveAngleThreshold;
		double orienThreshold = StyleSynthesisConfig::mAssemble_OrientationAngleThreshold;
		double angle = cml::deg(cml::acos_safe(fabs(srcOrientation.dot(tgtOrientation))));
		if (srcPrimitive == tgtPrimitive && (srcPrimitive != 2 || StyleSynthesisConfig::mAssemble_AllowSphereAlignment)) {
			// same primitive type
			if (angle > 90.0 - angleThreshold) {
				// align major axis
				//if (srcOrientation.dot(tgtOrientation) < 0) srcOrientation = -srcOrientation;
				Eigen::Vector3d rotAxis = srcOrientation.cross(tgtOrientation);
				if (rotAxis.squaredNorm()) {
					rotAxis.normalize();
				} else {
					if (fabs(srcOrientation[0]) || fabs(srcOrientation[2])) rotAxis = Eigen::Vector3d::UnitY();
					else rotAxis = Eigen::Vector3d::UnitZ();
				}
				double rotAngle = cml::acos_safe(srcOrientation.dot(tgtOrientation));
				if (StyleSynthesisConfig::mAssemble_AllowMinimalAlignment) {
					// applied rotation with minimal effort
					if (fabs(rotAngle) > cml::constantsd::pi_over_2()) {
						if (rotAngle > 0) rotAngle -= cml::constantsd::pi();
						else rotAngle += cml::constantsd::pi();
					}
				}
				Eigen::AngleAxisd rotation(rotAngle, rotAxis);
				if(rotAngle) preTransform.prerotate(rotation);
			}
		} else if(srcPrimitive != 2 && tgtPrimitive != 2) {
			// one is stick while other is plane
			if (angle < angleThreshold) {
				// push stick onto plane
				//Eigen::Vector3d rotAxis = (srcOrientation.cross(tgtOrientation)).normalized();
				//double rotAngle = cml::acos_safe(srcOrientation.dot(tgtOrientation)) - cml::rad(90.0);
				//Eigen::AngleAxisd rotation(rotAngle, rotAxis);
				//preTransform.prerotate(rotation);

				// align local CS
				Eigen::Matrix3Xd matS = sourceMat.colwise() - sourceMat.rowwise().mean();
				Eigen::Matrix3Xd matT = targetMat.colwise() - targetMat.rowwise().mean();
				Eigen::JacobiSVD< Eigen::Matrix3Xd > svdS(matS, Eigen::ComputeThinU);
				Eigen::JacobiSVD< Eigen::Matrix3Xd > svdT(matT, Eigen::ComputeThinU);
				Eigen::Matrix3d localCSS = svdS.matrixU();
				Eigen::Matrix3d localCST = svdT.matrixU();
				if (localCSS.col(0).cross(localCSS.col(1)).dot(localCSS.col(2)) < 0) {
					localCSS.col(2) = -localCSS.col(2);
				}
				if (localCST.col(0).cross(localCST.col(1)).dot(localCST.col(2)) < 0) {
					localCST.col(2) = -localCST.col(2);
				}
				if (localCSS.col(0).dot(localCST.col(0)) < 0) {
					localCSS.col(0) = -localCSS.col(0);
					localCSS.col(1) = -localCSS.col(1);
				}
				if (localCSS.col(1).dot(localCST.col(1)) < 0) {
					localCSS.col(1) = -localCSS.col(1);
					localCSS.col(2) = -localCSS.col(2);
				}
				Eigen::Matrix3d csXform = localCST * localCSS.transpose();
				preTransform.linear() = csXform;
				if (gDebugVisualization) {
					if (!MatchPrimitiveICP::visualize("stickPlane-1.ply", sourceMat, targetMat, preTransform)) return false;
					if (!MatchPrimitiveICP::visualize("stickPlane-2.ply", targetMat, sourceMat, preTransform.inverse())) return false;
				}
			}
		}
	}
//#ifdef OUTPUT_PROGRESS
//	cout << "pre-aligned orientation: " << srcOrientation.transpose() << endl;
//#endif
	sourceMat = preTransform * sourceMat;
	srcOrientation = preTransform * srcOrientation; // preTransform is guaranteed to be a rotation matrix
//#ifdef OUTPUT_PROGRESS
//	cout << "aligned orientation: " << srcOrientation.transpose() << endl;
//#endif

	// get transformation to place source part

	matchTransform.setIdentity();
	if (gDebugVisualization) {
		if (!MatchPrimitiveICP::visualize("prealignedPart.ply", sourceMat, targetMat, matchTransform)) return false;
	}
//#ifdef OUTPUT_PROGRESS
//	cout << "ICP setting:" << endl;
//	cout << "  mode = " << alignmentMode << endl;
//	cout << "  primitive = " << srcPrimitive << endl;
//	cout << "  orientation = " << srcOrientation.transpose() << endl;
//#endif
	if (!MatchPrimitiveICP::run(10,
		alignmentMode, srcPrimitive,
		sourceMat, targetMat,
		srcOrientation, matchTransform)) return false;
	if (gDebugVisualization) {
		if (!MatchPrimitiveICP::visualize("alignedPart.ply", sourceMat, targetMat, matchTransform)) return false;
	}
	matchTransform = matchTransform * preTransform;

	//fittingTransform = matchTransform; // HACK: skip fitting
	//return true;

	// get transformation to fit source part to target part (for matching slots)

	int fittingMode = alignmentMode;
	if (srcPrimitive == 0) { // stick
		fittingMode = max(2, fittingMode); // translation + scaling along major axis + uniform scaling of cross section
	} else if (srcPrimitive == 1) { // plane
		fittingMode = max(4, fittingMode); // translation + non-uniform scaling of cross section
	} else if (srcPrimitive == 2) { // sphere
		fittingMode = max(2, fittingMode); // translation + uniform scaling
	} else {
		return error("incorrect primitive type " + srcPrimitive);
	}

	if (!MatchPrimitiveICP::fit(
		fittingMode, srcPrimitive,
		sourceMat, targetMat,
		srcOrientation, fittingTransform)) return false;

	if (gDebugVisualization) {
		if (!MatchPrimitiveICP::visualize("fittedPart.ply", sourceMat, targetMat, fittingTransform)) return false;
		system("pause");
	}

	fittingTransform = fittingTransform * preTransform;

	return true;
}

bool ContextPartGraphAssembleUtil::alignMatchedSlots(
	Eigen::Matrix3Xd &inSourceFittingPoints,
	Eigen::Matrix3Xd &inSourceXformedPoints,
	vector<TSlot> &inSourceFittingSlots,
	vector<TSlot> &inSourceXformedSlots,
	vector<TSlot> &inTargetSlots,
	vector<TSlot> &outAlignedSlots)
{

	int numSourceSlots = (int)inSourceXformedSlots.size();
	int numTargetSlots = (int)inTargetSlots.size();

	double maxSlotNearestDist = StyleSynthesisConfig::mAssemble_SlotsMatchingNearestDistance;
	double maxSlotFarthestDist = StyleSynthesisConfig::mAssemble_SlotsMatchingFarthestDistance;
	double maxVirtualSlotDistance = StyleSynthesisConfig::mAssemble_SlotsSamplingRadius * 1.0;

	// initialize data for use in adding virtual slots

	SKDTree tree;
	SKDTreeData treeData;
	if (!SampleUtil::buildKdTree(inSourceFittingPoints, tree, treeData)) return false;

	double floorHeight = inSourceFittingPoints.row(1).minCoeff();
	double ceilingHeight = inSourceFittingPoints.row(1).maxCoeff();
	
	// match slots one by one

	outAlignedSlots.resize(numTargetSlots);

	for (int targetSlotID = 0; targetSlotID < numTargetSlots; targetSlotID++) {

		TSlot &targetSlot = inTargetSlots[targetSlotID];
		TSlot &alignedSlot = outAlignedSlots[targetSlotID];

		// get slot position

		if (targetSlot.adjacentIndex[0] < 0) {
			// wall slot
			bool foundSlot = false;
			for (int sourceSlotID = 0; sourceSlotID < numSourceSlots; sourceSlotID++) {
				TSlot &sourceSlot = inSourceFittingSlots[sourceSlotID];
				if (sourceSlot.adjacentIndex[0] == targetSlot.adjacentIndex[0]) {
					// use existing wall slot
					alignedSlot = inSourceXformedSlots[sourceSlotID];
					foundSlot = true;
					break;
				}
			}
			if (!foundSlot) {
				// add wall slot
				vector<Eigen::Vector3d> wallPoints;
				bool isFloor = (targetSlot.adjacentIndex[0] == -1);
				for (int pointID = 0; pointID < (int)inSourceFittingPoints.cols(); pointID++) {
					Eigen::Vector3d fittingPoint = inSourceFittingPoints.col(pointID);
					Eigen::Vector3d xformedPoint = inSourceXformedPoints.col(pointID);
					if (isFloor) {
						if (fabs(fittingPoint[1] - floorHeight) < maxVirtualSlotDistance) wallPoints.push_back(xformedPoint);
					} else {
						if (fabs(fittingPoint[1] - ceilingHeight) < maxVirtualSlotDistance) wallPoints.push_back(xformedPoint);
					}
				}
				alignedSlot.samples.resize(3, (int)wallPoints.size());
				for (int pointID = 0; pointID < (int)wallPoints.size(); pointID++) {
					alignedSlot.samples.col(pointID) = wallPoints[pointID];
				}
				alignedSlot.center = alignedSlot.samples.rowwise().mean();
				alignedSlot.isVirtualSlot = false;
				alignedSlot.isReUsedSlot = false;
			}
		} else {
			// find closest source slot
			double minDist = maxSlotNearestDist;
			int minSlotID = -1;
			for (int sourceSlotID = 0; sourceSlotID < numSourceSlots; sourceSlotID++) {
				TSlot &sourceSlot = inSourceFittingSlots[sourceSlotID];
				double slotMinDist, slotMaxDist;
				if (!computeSlotDistance(sourceSlot, targetSlot, slotMinDist, slotMaxDist)) return false;
				if (slotMinDist < minDist && slotMaxDist < maxSlotFarthestDist) {
					minDist = slotMinDist;
					minSlotID = sourceSlotID;
				}
			}

			if (minSlotID >= 0) {
				// use closest source slot
				alignedSlot = inSourceXformedSlots[minSlotID];
			} else {
				// add virtual slot by collecting nearest points
				set<int> nnSet;
				if (!SampleUtil::findNearestNeighbors(tree, targetSlot.samples, nnSet, maxVirtualSlotDistance)) return false;
				if (nnSet.empty()) {
					Eigen::VectorXi nnIndices;
					if (!SampleUtil::findNearestNeighbors(tree, targetSlot.samples, nnIndices)) return false;
					nnSet.insert(nnIndices.data(), nnIndices.data() + nnIndices.size());
				}
				vector<int> nnList(nnSet.begin(), nnSet.end());
				alignedSlot.samples.resize(3, (int)nnList.size());
				for (int nnID = 0; nnID < (int)nnList.size(); nnID++) {
					alignedSlot.samples.col(nnID) = inSourceXformedPoints.col(nnList[nnID]);
				}
				alignedSlot.center = alignedSlot.samples.rowwise().mean();
				alignedSlot.isVirtualSlot = true;
				alignedSlot.isReUsedSlot = false;
			}
		}

		// copy slot ID info

		alignedSlot.selfIndex = targetSlot.selfIndex;
		alignedSlot.adjacentIndex = targetSlot.adjacentIndex;
	}

	return true;
}

bool ContextPartGraphAssembleUtil::alignGlobalSlots(
	int inAlignmentMode,
	vector<vector<int>> &inAlignmentGroups,
	vector<int> &inAlignmentPrimitive,
	vector<Eigen::Vector3d> &inoutAlignmentOrientation,
	vector<vector<TSlot>> &inoutWorkingSlots,
	vector<Eigen::Affine3d> &outNodeTransformation)
{

	outNodeTransformation.resize(inoutWorkingSlots.size(), Eigen::Affine3d::Identity());

	int numGroups = (int)inAlignmentGroups.size();
	for (int groupID = 0; groupID < numGroups; groupID++) {

		//cout << "------------------------" << endl;

		// fix all other nodes not in group, transform nodes in the group to align slots

		vector<int> &group = inAlignmentGroups[groupID];
		unordered_set<int> groupSet(group.begin(), group.end());

		// get alignment info

		int mode = inAlignmentMode;
		int primitive = inAlignmentPrimitive[groupID];
		Eigen::Vector3d orientation = inoutAlignmentOrientation[groupID];
		if (primitive < 0) {
			// non-replaced nodes
			if (StyleSynthesisConfig::mAssemble_AllowGlobalAlignment) {
				// allow translation only
				mode = 0;
				primitive = 0;
				orientation = Eigen::Vector3d::UnitY();
			} else {
				// skip any transformation
				for (int nodeID : group) {
					outNodeTransformation[nodeID].setIdentity();
				}
				continue;
			}
		}
		/*
		// use best guess for each node according to primitives
		if (StyleSynthesisConfig::mAssemble_GlobalAlignmentMode == -1) {
			mode = StyleSynthesisConfig::mAssemble_BestGuessAlignmentMode.values[primitive];
		}
		*/

		// get slot data

		vector<Eigen::Vector3d> srcSlotPoints;
		vector<Eigen::Vector3d> tgtSlotPoints;
		vector<int> srcSlotLabels;
		vector<int> tgtSlotLabels;
		vector<double> labelWeights;

		int numSlotPairs = 0;
		for (int nodeID : group) {
			int numSlots = (int)inoutWorkingSlots[nodeID].size();
			bool alignVirtualSlot = true;
			// check whether should align virtual slots
			// if any "real" slot exists, skip aligning virtual slots
			for (int slotID = 0; slotID < numSlots; slotID++) {
				TSlot &slot = inoutWorkingSlots[nodeID][slotID];
				int nbNodeID = slot.adjacentIndex[0];
				int nbSlotID = slot.adjacentIndex[1];
				if (nbNodeID < 0) continue;
				TSlot &nbSlot = inoutWorkingSlots[nbNodeID][nbSlotID];
				if (!slot.isVirtualSlot && !nbSlot.isVirtualSlot) {
					alignVirtualSlot = false;
					break;
				}
			}
			for (int slotID = 0; slotID < numSlots; slotID++) {
				TSlot &slot = inoutWorkingSlots[nodeID][slotID];
				int nbNodeID = slot.adjacentIndex[0];
				int nbSlotID = slot.adjacentIndex[1];
				if (nbNodeID < 0) continue; // skip wall slots for global alignment
				if (groupSet.find(nbNodeID) != groupSet.end()) continue; // both nodes are in the group
				TSlot &nbSlot = inoutWorkingSlots[nbNodeID][nbSlotID];
				if (false) { // HACK: always align any kind of slots
				//if (slot.isVirtualSlot || nbSlot.isVirtualSlot) {
					if (alignVirtualSlot) {
						// align center for virtual slots
						srcSlotPoints.push_back(slot.center);
						srcSlotLabels.push_back(numSlotPairs);
						tgtSlotPoints.push_back(nbSlot.center);
						tgtSlotLabels.push_back(numSlotPairs);
					} else {
						continue;
					}
				} else {
					for (int id = 0; id < (int)slot.samples.cols(); id++) {
						srcSlotPoints.push_back(slot.samples.col(id));
						srcSlotLabels.push_back(numSlotPairs);
					}
					for (int id = 0; id < (int)nbSlot.samples.cols(); id++) {
						tgtSlotPoints.push_back(nbSlot.samples.col(id));
						tgtSlotLabels.push_back(numSlotPairs);
					}
				}
				if (slot.isReUsedSlot || nbSlot.isReUsedSlot) {
					labelWeights.push_back(StyleSynthesisConfig::mAssemble_SlotsReusingWeightFactor);
				} else {
					labelWeights.push_back(1.0);
				}
				//cout << "Slots: " << slot.selfIndex << " " << (slot.isVirtualSlot ? "T" : "F") << (slot.isReUsedSlot ? "T" : "F");
				//cout << " ; " << nbSlot.selfIndex << " " << (nbSlot.isVirtualSlot ? "T" : "F") << (nbSlot.isReUsedSlot ? "T" : "F") << endl;
				numSlotPairs++;
			}
		}

		// build matrices for labeled ICP

		Eigen::Matrix3Xd matSourcePoints;
		Eigen::Matrix3Xd matTargetPoints;
		Eigen::VectorXi vecSourceLabels;
		Eigen::VectorXi vecTargetLabels;
		Eigen::VectorXd vecLabelWeights;
		if (!MatchLabeledICP::buildMatrices(srcSlotPoints, srcSlotLabels, matSourcePoints, vecSourceLabels)) return false;
		if (!MatchLabeledICP::buildMatrices(tgtSlotPoints, tgtSlotLabels, matTargetPoints, vecTargetLabels)) return false;
		vecLabelWeights = Eigen::Map<Eigen::VectorXd>(labelWeights.data(), labelWeights.size());

		// apply labeled ICP

		Eigen::Affine3d xform;
		xform.setIdentity();
		if (gDebugVisualization) {
			if (!MatchLabeledICP::visualize(
				"ICP-before.ply",
				matSourcePoints, matTargetPoints,
				vecSourceLabels, vecTargetLabels,
				xform)) return false;
		}
		if (!MatchLabeledICP::run(10, // UNDONE: param number of ICP iterations
			mode, primitive,
			matSourcePoints, matTargetPoints,
			vecSourceLabels, vecTargetLabels,
			orientation, xform, &vecLabelWeights)) return false;
		if (gDebugVisualization) {
			cout << "Orientation: " << orientation.transpose() << endl;
			cout << "Xform:\n" << xform.matrix() << endl;
			if (!MatchLabeledICP::visualize(
				"ICP-after.ply",
				matSourcePoints, matTargetPoints,
				vecSourceLabels, vecTargetLabels,
				xform)) return false;
			system("pause");
		}

		// apply transformation on this group

		for (int nodeID : group) {
			for (TSlot &slot : inoutWorkingSlots[nodeID]) {
				slot.samples = xform * slot.samples;
				slot.center = xform * slot.center;
			}
			outNodeTransformation[nodeID] = xform;
		}

		if (inAlignmentPrimitive[groupID] >= 0) {
			inoutAlignmentOrientation[groupID] = xform.rotation() * inoutAlignmentOrientation[groupID];
		}
	}

	return true;
}

bool ContextPartGraphAssembleUtil::alignGlobalSlots(
	vector<int> &inWorkingNodes,
	vector<vector<TSlot>> &inoutWorkingSlots,
	vector<Eigen::Affine3d> &outNodeTransformation)
{
	int numAllNodes = (int)inoutWorkingSlots.size();

	outNodeTransformation.resize(inoutWorkingSlots.size());
	for (int nodeID = 0; nodeID < numAllNodes; nodeID++) {
		outNodeTransformation[nodeID].setIdentity();
	}

	for (int nodeID : inWorkingNodes) {

		// fix all other nodes, transform current node to align slots

		// get slot data

		vector<Eigen::Vector3d> srcSlotPoints;
		vector<Eigen::Vector3d> tgtSlotPoints;
		vector<int> srcSlotLabels;
		vector<int> tgtSlotLabels;
		vector<double> labelWeights;

		int numSlotPairs = 0;
		for (int slotID = 0; slotID < (int)inoutWorkingSlots[nodeID].size(); slotID++) {
			TSlot &slot = inoutWorkingSlots[nodeID][slotID];
			int nbNodeID = slot.adjacentIndex[0];
			int nbSlotID = slot.adjacentIndex[1];
			if (nbNodeID < 0) continue; // skip wall slots for global alignment
			Eigen::Affine3d nbXform = outNodeTransformation[nbNodeID];
			TSlot &nbSlot = inoutWorkingSlots[nbNodeID][nbSlotID];
			for (int id = 0; id < (int)slot.samples.cols(); id++) {
				srcSlotPoints.push_back(slot.samples.col(id));
				srcSlotLabels.push_back(numSlotPairs);
			}
			Eigen::Matrix3Xd nbSlotSamples = nbXform * nbSlot.samples;
			for (int id = 0; id < (int)nbSlot.samples.cols(); id++) {
				tgtSlotPoints.push_back(nbSlotSamples.col(id));
				tgtSlotLabels.push_back(numSlotPairs);
			}
			if (slot.isReUsedSlot || nbSlot.isReUsedSlot) {
				labelWeights.push_back(5.0); // UNDONE: param re-used slots label weight
			} else {
				labelWeights.push_back(1.0);
			}
			//cout << "Slots: " << slot.selfIndex << "; " << nbSlot.selfIndex << endl;
			numSlotPairs++;
		}

		// build matrices for labeled ICP

		Eigen::Matrix3Xd matSourcePoints;
		Eigen::Matrix3Xd matTargetPoints;
		Eigen::VectorXi vecSourceLabels;
		Eigen::VectorXi vecTargetLabels;
		Eigen::VectorXd vecLabelWeights;
		if (!MatchLabeledICP::buildMatrices(srcSlotPoints, srcSlotLabels, matSourcePoints, vecSourceLabels)) return false;
		if (!MatchLabeledICP::buildMatrices(tgtSlotPoints, tgtSlotLabels, matTargetPoints, vecTargetLabels)) return false;
		vecLabelWeights = Eigen::Map<Eigen::VectorXd>(labelWeights.data(), labelWeights.size());

		// apply labeled ICP

		int mode = 0; // only allow translation & uniform scaling
		int primitive = 2; // sphere
		Eigen::Vector3d orientation = Eigen::Vector3d::UnitY();

		Eigen::Affine3d xform;
		xform.setIdentity();
		if (gDebugVisualization) {
			if (!MatchLabeledICP::visualize(
				"ICP-before.ply",
				matSourcePoints, matTargetPoints,
				vecSourceLabels, vecTargetLabels,
				xform)) return false;
		}
		if (!MatchLabeledICP::run(10, // UNDONE: param number of ICP iterations
			mode, primitive,
			matSourcePoints, matTargetPoints,
			vecSourceLabels, vecTargetLabels,
			orientation, xform, &vecLabelWeights)) return false;
		if (gDebugVisualization) {
			if (!MatchLabeledICP::visualize(
				"ICP-after.ply",
				matSourcePoints, matTargetPoints,
				vecSourceLabels, vecTargetLabels,
				xform)) return false;
			system("pause");
		}

		// apply transformation on this group

		for (TSlot &slot : inoutWorkingSlots[nodeID]) {
			slot.samples = xform * slot.samples;
			slot.center = xform * slot.center;
		}
		outNodeTransformation[nodeID] = xform;
	}

	return true;
}

bool ContextPartGraphAssembleUtil::alignExemplarSlots(
	vector<vec2i> &replaceMapping,
	vector<Eigen::Affine3d> &exemplarXform,
	vector<Eigen::Affine3d> &candidateXform,
	vector<vector<TSlot>> &exemplarSlots,
	vector<vector<TSlot>> &candidateSlots,
	vector<vector<TSlot>> &alignedSlots)
{
	int numSourceNodes = (int)exemplarSlots.size();
	int numTargetNodes = (int)replaceMapping.size();

	vector<vector<int>> srcNodeMappings(numSourceNodes, vector<int>(0)); // target ID : # of target nodes replaced by this source node : # of source nodes
	for (int targetID = 0; targetID < numTargetNodes; targetID++) {
		int sourceID = replaceMapping[targetID][0];
		if (sourceID < 0) continue;
		srcNodeMappings[sourceID].push_back(targetID);
	}

	set<vec2i> processedNodePairs;

	alignedSlots.assign(numTargetNodes, vector<TSlot>(0));

	// process existing exemplar slots

	for (int targetID = 0; targetID < numTargetNodes; targetID++) {
		int sourceID = replaceMapping[targetID][0];
		if (sourceID < 0) continue;

		int numSlots = (int)exemplarSlots[sourceID].size();
		for (int srcSlotID = 0; srcSlotID < numSlots; srcSlotID++) {
			TSlot &srcSlot = exemplarSlots[sourceID][srcSlotID];
			int nbSourceID = srcSlot.adjacentIndex[0];
			if (nbSourceID < 0) continue;
			TSlot &nbSrcSlot = exemplarSlots[nbSourceID][srcSlot.adjacentIndex[1]];
			for (int nbTargetID : srcNodeMappings[nbSourceID]) {
				vec2i key(targetID, nbTargetID);
				if (key[0] > key[1]) swap(key[0], key[1]);
				if (processedNodePairs.find(key) != processedNodePairs.end()) continue;

				// transform slots

				TSlot tgtSlot, nbTgtSlot;

				tgtSlot.samples = exemplarXform[targetID] * srcSlot.samples;
				tgtSlot.center = tgtSlot.samples.rowwise().mean();
				tgtSlot.isVirtualSlot = false;
				tgtSlot.isReUsedSlot = true;
				tgtSlot.selfIndex[0] = targetID;
				tgtSlot.selfIndex[1] = (int)alignedSlots[targetID].size();

				nbTgtSlot.samples = exemplarXform[nbTargetID] * nbSrcSlot.samples;
				nbTgtSlot.center = nbTgtSlot.samples.rowwise().mean();
				nbTgtSlot.isVirtualSlot = false;
				nbTgtSlot.isReUsedSlot = true;
				nbTgtSlot.selfIndex[0] = nbTargetID;
				nbTgtSlot.selfIndex[1] = (int)alignedSlots[nbTargetID].size();

				// check adjacency

				if (true) {
					double distThr = StyleSynthesisConfig::mAssemble_SlotsReusingDistanceThreshold;

					SKDTree tgtTree, nbTgtTree;
					SKDTreeData tgtTreeData, nbTgtTreeData;
					if (!SampleUtil::buildKdTree(tgtSlot.samples, tgtTree, tgtTreeData)) return false;
					if (!SampleUtil::buildKdTree(nbTgtSlot.samples, nbTgtTree, nbTgtTreeData)) return false;
					Eigen::VectorXi tgtNNIndices, nbTgtNNIndices;
					if (!SampleUtil::findNearestNeighbors(nbTgtTree, tgtSlot.samples, tgtNNIndices)) return false;
					if (!SampleUtil::findNearestNeighbors(tgtTree, nbTgtSlot.samples, nbTgtNNIndices)) return false;

					vector<int> tgtSet, nbTgtSet;
					double minDist = DBL_MAX;
					for (int id = 0; id < (int)tgtNNIndices.size(); id++) {
						int nbID = tgtNNIndices[id];
						if (nbTgtNNIndices[nbID] == id) { // REAL nearest neighbors
							tgtSet.push_back(id);
							nbTgtSet.push_back(nbID);
							double dist = (tgtSlot.samples.col(id) - nbTgtSlot.samples.col(nbID)).norm();
							minDist = min(minDist, dist);
						}
					}
					if (minDist > distThr) continue; // slots are too far away

					double centerDist = (tgtSlot.center - nbTgtSlot.center).norm();
					if (centerDist < distThr) {
						// just use these slots
					} else {
						// keep nearest points only
						Eigen::Matrix3Xd tgtMat(3, (int)tgtSet.size());
						Eigen::Matrix3Xd nbTgtMat(3, (int)nbTgtSet.size());
						for (int id = 0; id < (int)tgtSet.size(); id++) {
							tgtMat.col(id) = tgtSlot.samples.col(tgtSet[id]);
						}
						for (int id = 0; id < (int)nbTgtSet.size(); id++) {
							nbTgtMat.col(id) = nbTgtSlot.samples.col(nbTgtSet[id]);
						}
						tgtSlot.samples.swap(tgtMat);
						nbTgtSlot.samples.swap(nbTgtMat);
						tgtSlot.center = tgtSlot.samples.rowwise().mean();
						nbTgtSlot.center = nbTgtSlot.samples.rowwise().mean();

						tgtSlot.isReUsedSlot = false;
						tgtSlot.isVirtualSlot = true;

						nbTgtSlot.isReUsedSlot = false;
						nbTgtSlot.isVirtualSlot = true;
					}
				}

				// add slots

				tgtSlot.adjacentIndex = nbTgtSlot.selfIndex;
				nbTgtSlot.adjacentIndex = tgtSlot.selfIndex;

				alignedSlots[targetID].push_back(tgtSlot);
				alignedSlots[nbTargetID].push_back(nbTgtSlot);

				processedNodePairs.insert(key);
			}
		}
	}

	if (StyleSynthesisConfig::mAssemble_PostAlignCandidateSlots) {

		// process "unprocessed" candidate slots

		for (int targetID = 0; targetID < numTargetNodes; targetID++) {
			int numSlots = (int)candidateSlots[targetID].size();
			for (int tgtSlotID = 0; tgtSlotID < numSlots; tgtSlotID++) {
				TSlot &tgtSlot = candidateSlots[targetID][tgtSlotID];
				int nbTargetID = tgtSlot.adjacentIndex[0];
				if (nbTargetID < 0) continue;
				TSlot &nbTgtSlot = candidateSlots[nbTargetID][tgtSlot.adjacentIndex[1]];

				vec2i key(targetID, nbTargetID);
				if (key[0] > key[1]) swap(key[0], key[1]);
				if (processedNodePairs.find(key) != processedNodePairs.end()) continue;
				processedNodePairs.insert(key);

				// transform slots

				TSlot newTgtSlot, newNbTgtSlot;

				newTgtSlot.samples = candidateXform[targetID] * tgtSlot.samples;
				newTgtSlot.center = newTgtSlot.samples.rowwise().mean();
				newTgtSlot.isVirtualSlot = true;
				newTgtSlot.isReUsedSlot = false;
				newTgtSlot.selfIndex[0] = targetID;
				newTgtSlot.selfIndex[1] = (int)alignedSlots[targetID].size();

				newNbTgtSlot.samples = candidateXform[nbTargetID] * nbTgtSlot.samples;
				newNbTgtSlot.center = newNbTgtSlot.samples.rowwise().mean();
				newNbTgtSlot.isVirtualSlot = true;
				newNbTgtSlot.isReUsedSlot = false;
				newNbTgtSlot.selfIndex[0] = nbTargetID;
				newNbTgtSlot.selfIndex[1] = (int)alignedSlots[nbTargetID].size();

				// add slots

				newTgtSlot.adjacentIndex = newNbTgtSlot.selfIndex;
				newNbTgtSlot.adjacentIndex = newTgtSlot.selfIndex;

				alignedSlots[targetID].push_back(newTgtSlot);
				alignedSlots[nbTargetID].push_back(newNbTgtSlot);
			}
		}
	}
	
	return true;
}

bool ContextPartGraphAssembleUtil::alignSlotPairs(
	int inAlignmentMode,
	int inAlignmentPrimitive,
	Eigen::Vector3d &inAlignmentOrientation,
	vector<TSlot> &inSourceSlots,
	vector<TSlot> &inTargetSlots,
	Eigen::Affine3d &outTransformation,
	double &outError)
{

	// gather slots data

	vector<Eigen::Vector3d> srcSlotPoints;
	vector<Eigen::Vector3d> tgtSlotPoints;
	vector<int> srcSlotLabels;
	vector<int> tgtSlotLabels;

	int numSlots = (int)inTargetSlots.size();
	for (int slotID = 0; slotID < numSlots; slotID++) {
		for (int pointID = 0; pointID < (int)inSourceSlots[slotID].samples.cols(); pointID++) {
			srcSlotPoints.push_back(inSourceSlots[slotID].samples.col(pointID));
			srcSlotLabels.push_back(slotID);
		}
		for (int pointID = 0; pointID < (int)inTargetSlots[slotID].samples.cols(); pointID++) {
			tgtSlotPoints.push_back(inTargetSlots[slotID].samples.col(pointID));
			tgtSlotLabels.push_back(slotID);
		}
	}

	// build matrices

	Eigen::Matrix3Xd matSourcePoints;
	Eigen::Matrix3Xd matTargetPoints;
	Eigen::VectorXi vecSourceLabels;
	Eigen::VectorXi vecTargetLabels;
	if (!MatchLabeledICP::buildMatrices(srcSlotPoints, srcSlotLabels, matSourcePoints, vecSourceLabels)) return false;
	if (!MatchLabeledICP::buildMatrices(tgtSlotPoints, tgtSlotLabels, matTargetPoints, vecTargetLabels)) return false;

	// apply labeled ICP

	Eigen::Affine3d xform;
	xform.setIdentity();
	if (gDebugVisualization) {
		if (!MatchLabeledICP::visualize(
			"ICP-before.ply",
			matSourcePoints, matTargetPoints,
			vecSourceLabels, vecTargetLabels,
			xform)) return false;
	}
	if (!MatchLabeledICP::run(10, // UNDONE: param number of ICP iterations
		inAlignmentMode, inAlignmentPrimitive,
		matSourcePoints, matTargetPoints,
		vecSourceLabels, vecTargetLabels,
		inAlignmentOrientation, xform)) return false;
	if (gDebugVisualization) {
		if (!MatchLabeledICP::visualize(
			"ICP-after.ply",
			matSourcePoints, matTargetPoints,
			vecSourceLabels, vecTargetLabels,
			xform)) return false;
		//system("pause");
	}

	outTransformation = xform;

	// compute error

	matSourcePoints = xform * matSourcePoints;

	Eigen::VectorXd labelErrorST, labelErrorTS;
	if (!MatchLabeledICP::error(
		matSourcePoints, matTargetPoints,
		vecSourceLabels, vecTargetLabels,
		labelErrorST)) return false;
	if (!MatchLabeledICP::error(
		matTargetPoints, matSourcePoints,
		vecTargetLabels, vecSourceLabels,
		labelErrorTS)) return false;
	Eigen::VectorXd labelError = (labelErrorST + labelErrorTS) / 2;

	outError = labelError.mean();

#ifdef OUTPUT_PROGRESS
	cout << "Slots error: " << outError << endl;
#endif

	return true;
}

bool ContextPartGraphAssembleUtil::updateVirtualSlots(
	vector<Eigen::Matrix3Xd> &inSourcePoints,
	vector<vec2i> &inMatchNodePairs,
	vector<vector<TSlot>> &inWorkingSlots)
{

	int numMatchPairs = (int)inMatchNodePairs.size();

	double maxVirtualSlotDistance = StyleSynthesisConfig::mAssemble_SlotsSamplingRadius * 1.0;

	for (int pairID = 0; pairID < numMatchPairs; pairID++) {
		int sourceID = inMatchNodePairs[pairID][0];
		int targetID = inMatchNodePairs[pairID][1];

		Eigen::Matrix3Xd &sourceMat = inSourcePoints[pairID];
		SKDTree tree;
		SKDTreeData treeData;
		if (!SampleUtil::buildKdTree(sourceMat, tree, treeData)) return false;

		for (TSlot &slot : inWorkingSlots[targetID]) {
			if (!slot.isVirtualSlot) continue;
			TSlot &nbSlot = inWorkingSlots[slot.adjacentIndex[0]][slot.adjacentIndex[1]];
			//if (nbSlot.isVirtualSlot) continue; // if both are virtual slots, do not update

			// update virtual slot with nearest points
			set<int> nnSet;
			if (!SampleUtil::findNearestNeighbors(tree, nbSlot.samples, nnSet, maxVirtualSlotDistance)) return false;
			if (nnSet.empty()) {
				Eigen::VectorXi nnIndices;
				if (!SampleUtil::findNearestNeighbors(tree, nbSlot.samples, nnIndices)) return false;
				nnSet.insert(nnIndices.data(), nnIndices.data() + nnIndices.size());
			}
			vector<int> nnList(nnSet.begin(), nnSet.end());
			slot.samples.resize(3, (int)nnList.size());
			for (int nnID = 0; nnID < (int)nnList.size(); nnID++) {
				slot.samples.col(nnID) = sourceMat.col(nnList[nnID]);
			}
			slot.center = slot.samples.rowwise().mean();
		}
	}

	return true;
}

bool ContextPartGraphAssembleUtil::checkSlotsAlignment(
	vector<vector<int>> &inAlignmentGroups,
	vector<vector<TSlot>> &inWorkingSlots,
	Eigen::VectorXd &outError)
{

	// check slots alignment (only check slot pairs across groups)

	int numGroups = (int)inAlignmentGroups.size();
	int numNodes = (int)inWorkingSlots.size();

	// compute group mapping

	vector<int> nodeGroupMapping(numNodes, -1);
	vector<int> workingNodes(0);
	for (int groupID = 0; groupID < numGroups; groupID++) {
		vector<int> &group = inAlignmentGroups[groupID];
		for (int nodeID : group) {
			nodeGroupMapping[nodeID] = groupID;
			workingNodes.push_back(nodeID);
		}
	}

	// get floor / ceiling height

	Eigen::AlignedBox3d slotsBB;
	for (int nodeID : workingNodes) {
		for (TSlot &slot : inWorkingSlots[nodeID]) {
			Eigen::Vector3d minPoint = slot.samples.rowwise().minCoeff();
			Eigen::Vector3d maxPoint = slot.samples.rowwise().maxCoeff();
			slotsBB.extend(minPoint);
			slotsBB.extend(maxPoint);
		}
	}
	double floorHeight = slotsBB.min()[1];
	double ceilingHeight = slotsBB.max()[1];

	// get all slot pairs

	vector<vector<bool>> slotFlags(numNodes); // visited flag to prevent counting each slot twice
	for (int nodeID : workingNodes) {
		int numSlots = (int)inWorkingSlots[nodeID].size();
		slotFlags[nodeID].assign(numSlots, false);
	}

	vector<Eigen::Vector3d> srcSlotPoints;
	vector<Eigen::Vector3d> tgtSlotPoints;
	vector<int> srcSlotLabels;
	vector<int> tgtSlotLabels;
	vector<double> wallErrorList;

	int numSlotPairs = 0;
	vector<vec2i> slotPairs(0);
	for (int nodeID : workingNodes) {
		int numSlots = (int)inWorkingSlots[nodeID].size();
		for (int slotID = 0; slotID < numSlots; slotID++) {
			if (slotFlags[nodeID][slotID]) continue;
			TSlot &slot = inWorkingSlots[nodeID][slotID];
			int nbNodeID = slot.adjacentIndex[0];
			int nbSlotID = slot.adjacentIndex[1];
			if (nbNodeID == -1) {
				double floorError = fabs(slot.samples.row(1).minCoeff() - floorHeight);
				wallErrorList.push_back(cml::sqr(floorError));
			} else if (nbNodeID == -2) {
				double ceilingError = fabs(slot.samples.row(1).maxCoeff() - ceilingHeight);
				wallErrorList.push_back(cml::sqr(ceilingError));
			} else {
				if (nodeGroupMapping[nodeID] == nodeGroupMapping[nbNodeID]) continue; // both nodes are in the same group
				TSlot &nbSlot = inWorkingSlots[nbNodeID][nbSlotID];
				slotFlags[nbNodeID][nbSlotID] = true;
				if (true) {
					for (int pointID = 0; pointID < (int)slot.samples.cols(); pointID++) {
						srcSlotPoints.push_back(slot.samples.col(pointID));
						srcSlotLabels.push_back(numSlotPairs);
					}
					for (int pointID = 0; pointID < (int)nbSlot.samples.cols(); pointID++) {
						tgtSlotPoints.push_back(nbSlot.samples.col(pointID));
						tgtSlotLabels.push_back(numSlotPairs);
					}
				}
				slotPairs.push_back(vec2i(nodeID, slotID));
				numSlotPairs++;
			}
		}
	}
	if (wallErrorList.empty()) wallErrorList.push_back(0.0);
	Eigen::VectorXd wallError = Eigen::Map<Eigen::VectorXd>(wallErrorList.data(), wallErrorList.size());

	if (srcSlotPoints.empty() || tgtSlotPoints.empty()) {
		outError.resize(1);
		outError[0] = wallError.maxCoeff();
		return true;
	}

	// compute labeled ICP error

	Eigen::VectorXd labelError;
	if (true) {
		Eigen::Matrix3Xd matSourcePoints;
		Eigen::Matrix3Xd matTargetPoints;
		Eigen::VectorXi vecSourceLabels;
		Eigen::VectorXi vecTargetLabels;
		if (!MatchLabeledICP::buildMatrices(srcSlotPoints, srcSlotLabels, matSourcePoints, vecSourceLabels)) return false;
		if (!MatchLabeledICP::buildMatrices(tgtSlotPoints, tgtSlotLabels, matTargetPoints, vecTargetLabels)) return false;

		Eigen::VectorXd labelErrorST, labelErrorTS;
		if (!MatchLabeledICP::error(
			matSourcePoints, matTargetPoints,
			vecSourceLabels, vecTargetLabels,
			labelErrorST)) return false;
		if (!MatchLabeledICP::error(
			matTargetPoints, matSourcePoints,
			vecTargetLabels, vecSourceLabels,
			labelErrorTS)) return false;
		labelError = (labelErrorST + labelErrorTS) / 2;
		//labelError = labelErrorST.cwiseMin(labelErrorTS);
	}

	// compute final error

	//outError.resize(labelError.size() + wallError.size());
	//outError << labelError, wallError;
	outError.resize(labelError.size() + 1);
	outError << labelError, wallError.maxCoeff();

	return true;
}

bool ContextPartGraphAssembleUtil::regularizeTransformation(
	Eigen::Affine3d &inTransform,
	Eigen::Affine3d &outTransform,
	double eps)
{

	// decompose affine transformation

	Eigen::Vector3d vecT = inTransform.translation();
	Eigen::Matrix3d matLinear = inTransform.linear();
	Eigen::JacobiSVD<Eigen::Matrix3d> svd(matLinear, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Vector3d vecS = svd.singularValues();

	/*
	cout << "====================" << endl;
	cout << "before: " << endl;
	cout << "\t T: " << vecT.transpose() << endl;
	cout << "\t S: " << vecS.transpose() << endl;
	*/

	// regularize

	for (int dim = 0; dim < 3; dim++) {
		if (fabs(vecT[dim]) < eps) vecT[dim] = 0.0;
		if (fabs(max(vecS[dim], 1.0 / vecS[dim]) - 1.0) < 1e-3) vecS[dim] = 1.0; // UNDONE: param scale factor threshold
	}

	/*
	cout << "after: " << endl;
	cout << "\t T: " << vecT.transpose() << endl;
	cout << "\t S: " << vecS.transpose() << endl;
	*/

	// output regularized transformation

	outTransform.setIdentity();
	outTransform.linear() = svd.matrixU() * vecS.asDiagonal() * svd.matrixV().adjoint();
	outTransform.pretranslate(vecT);

	return true;
}

bool ContextPartGraphAssembleUtil::regularizeNodeTransformation(
	Eigen::Affine3d &inTransform,
	Eigen::Affine3d &outTransform,
	Eigen::Vector3d &inCenter, double eps)
{
	double maxRegScale = 1.0 + eps;
	double minRegScale = 1.0 / maxRegScale;
	double maxRegAngle = 30.0;

	outTransform.setIdentity();
	
	Eigen::Matrix3d matLinear = inTransform.linear();
	Eigen::JacobiSVD<Eigen::Matrix3d> svd(matLinear, Eigen::ComputeFullU | Eigen::ComputeFullV);

	// regularize scaling
	Eigen::Vector3d scale = svd.singularValues();
	for (int dim = 0; dim < 3; dim++) {
		double oldScale = scale[dim];
		if (oldScale < maxRegScale && oldScale > minRegScale) {
			scale[dim] = 1.0;
		} else {
			if (StyleSynthesisConfig::mAssemble_PostRegularizationRounding) {
				if (oldScale >= maxRegScale) {
					scale[dim] = ceil(oldScale - eps);
				} else {
					scale[dim] = 1.0 / ceil(1.0 / oldScale - eps);
				}
			}
		}
#ifdef OUTPUT_PROGRESS
		cout << "Regularized scaling " << dim << " : " << oldScale << " => " << scale[dim] << endl;
#endif
	}

	// regularize rotation
	Eigen::Matrix3d rotation = svd.matrixU() * svd.matrixV().adjoint();
	Eigen::AngleAxisd aa(rotation);
	if (cml::deg(aa.angle()) < maxRegAngle) {
#ifdef OUTPUT_PROGRESS
		cout << "Regularized rotation : " << cml::deg(aa.angle()) << " => " << "0.0" << endl;
#endif
		outTransform.linear() = svd.matrixU() * scale.asDiagonal() * svd.matrixU().adjoint();
	} else {
#ifdef OUTPUT_PROGRESS
		cout << "Kept rotation : " << cml::deg(aa.angle()) << endl;
#endif
		outTransform.linear() = svd.matrixU() * scale.asDiagonal() * svd.matrixV().adjoint();
	}

	// get translation
	Eigen::Vector3d srcCenter = inCenter;
	Eigen::Vector3d tgtCenter = inTransform * srcCenter;
	outTransform.translate(-srcCenter);
	outTransform.pretranslate(tgtCenter);

	return true;
}