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

#include "ContextPartGraphTrain.h"

#include <iostream>
#include <fstream>
#include <set>

#include "Mesh/MeshUtil.h"
#include "Sample/SampleUtil.h"
#include "Segment/SegmentUtil.h"

#include "Match/MatchSimpleICP.h"

#include "Context/ContextPartGraphMatch.h"

#include "Utility/PlyExporter.h"

#include "Data/DataUtil.h"
#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

bool gUseFirstLevelNodesOnly = true; // UNDONE: param use first node only for learning

ContextPartGraphTrain::ContextPartGraphTrain() {
}

ContextPartGraphTrain::~ContextPartGraphTrain() {
}

bool ContextPartGraphTrain::loadSceneMesh(vector<TTriangleMesh> &sceneMeshes) {
	mpSceneMeshes = &sceneMeshes;
	return true;
}
bool ContextPartGraphTrain::loadSceneSegments(vector<vector<vector<vector<int>>>> &sceneSegments) {
	mpSceneSegments = &sceneSegments;
	return true;
}
bool ContextPartGraphTrain::loadSceneGraph(vector<TGraph> &sceneGraphs) {
	mpSceneGraphs = &sceneGraphs;
	return true;
}

bool ContextPartGraphTrain::process() {

	if (!generateElementGroups()) return false;
	if (!generateSymmetryGroups()) return false;
	if (!extractTrainPartPairs()) return false;

	return true;
}

bool ContextPartGraphTrain::generateElementGroups() {

	int numMeshes = (int)mpSceneMeshes->size();

	vector<vector<bool>> nodeVisitedFlags(numMeshes);
	for (int meshID = 0; meshID < numMeshes; meshID++) {
		int numNodes = (int)(*mpSceneGraphs)[meshID].mAllNodes.size();
		nodeVisitedFlags[meshID].assign(numNodes, false);
	}

	mElementGroup.clear();

	for (int meshID = 0; meshID < numMeshes; meshID++) {
		
		int numNodes;
		if (!gUseFirstLevelNodesOnly) {
			numNodes = (int)(*mpSceneGraphs)[meshID].mAllNodes.size();
		} else {
			numNodes = (int)(*mpSceneGraphs)[meshID].mRootNode->mChildren.size();
		}

		for (int nodeID = 0; nodeID < numNodes; nodeID++) {
			if (nodeVisitedFlags[meshID][nodeID]) continue;
			nodeVisitedFlags[meshID][nodeID] = true;

			cout << "\rChecking mesh " << (meshID + 1) << "/" << numMeshes << " node " << (nodeID + 1) << "/" << numNodes << "            ";

			vector<vec2i> elementGroup(1, vec2i(meshID, nodeID));
			TNode *thisNode = (*mpSceneGraphs)[meshID].mAllNodes[nodeID];

			// check symmetric nodes

			for (TNode *otherNode : thisNode->mSymmetry) {
				int otherNodeID = otherNode->mUID;
				elementGroup.push_back(vec2i(meshID, otherNodeID));
				nodeVisitedFlags[meshID][otherNodeID] = true;
			}

			// check nodes in other meshes

			for (int otherMeshID = meshID + 1; otherMeshID < numMeshes; otherMeshID++) {
				vector<bool> checkedFlags = nodeVisitedFlags[otherMeshID];

				int numOtherNodes;
				if (!gUseFirstLevelNodesOnly) {
					numOtherNodes = (int)(*mpSceneGraphs)[otherMeshID].mAllNodes.size();
				} else {
					numOtherNodes = (int)(*mpSceneGraphs)[otherMeshID].mRootNode->mChildren.size();
				}

				for (int otherNodeID = 0; otherNodeID < numOtherNodes; otherNodeID++) {
					if (checkedFlags[otherNodeID]) continue;
					checkedFlags[otherNodeID] = true;

					// check node
					TNode *otherNode = (*mpSceneGraphs)[otherMeshID].mAllNodes[otherNodeID];
					bool isElement;
					if (!checkNodes(thisNode, otherNode, isElement)) return false;

					// add element
					if (isElement) {
						elementGroup.push_back(vec2i(otherMeshID, otherNodeID));
						nodeVisitedFlags[otherMeshID][otherNodeID] = true;
					}

					// process symmetric nodes
					for (TNode *otherSymNode : otherNode->mSymmetry) {
						int otherSymNodeID = otherSymNode->mUID;
						if (isElement) {
							elementGroup.push_back(vec2i(otherMeshID, otherSymNodeID));
							nodeVisitedFlags[otherMeshID][otherSymNodeID] = true;
						}
						checkedFlags[otherSymNodeID] = true;
					}
				}
			}

			// add element group

			if ((int)elementGroup.size() > 1) {
				int anyMeshID = elementGroup[0][0];
				bool acrossMeshes = false;
				for (vec2i element : elementGroup) {
					if (element[0] != anyMeshID) acrossMeshes = true;
				}
				if (acrossMeshes) mElementGroup.push_back(elementGroup);
			}
		}
	}
	cout << endl;

	return true;
}

bool ContextPartGraphTrain::generateSymmetryGroups() {

	int numMeshes = (int)mpSceneMeshes->size();

	vector<vector<bool>> nodeVisitedFlags(numMeshes);
	for (int meshID = 0; meshID < numMeshes; meshID++) {
		int numNodes = (int)(*mpSceneGraphs)[meshID].mAllNodes.size();
		nodeVisitedFlags[meshID].assign(numNodes, false);
	}

	mSymmetryGroup.clear();

	for (int meshID = 0; meshID < numMeshes; meshID++) {

		int numNodes;
		if (!gUseFirstLevelNodesOnly) {
			numNodes = (int)(*mpSceneGraphs)[meshID].mAllNodes.size();
		} else {
			numNodes = (int)(*mpSceneGraphs)[meshID].mRootNode->mChildren.size();
		}

		for (int nodeID = 0; nodeID < numNodes; nodeID++) {
			if (nodeVisitedFlags[meshID][nodeID]) continue;
			nodeVisitedFlags[meshID][nodeID] = true;

			cout << "\rChecking mesh " << (meshID + 1) << "/" << numMeshes << " node " << (nodeID + 1) << "/" << numNodes << "            ";

			vector<vec2i> elementGroup(1, vec2i(meshID, nodeID));
			TNode *thisNode = (*mpSceneGraphs)[meshID].mAllNodes[nodeID];

			// check symmetric nodes

			for (TNode *otherNode : thisNode->mSymmetry) {
				int otherNodeID = otherNode->mUID;
				elementGroup.push_back(vec2i(meshID, otherNodeID));
				nodeVisitedFlags[meshID][otherNodeID] = true;
			}

			// add symmetry group

			if ((int)elementGroup.size() > 1) {
				mSymmetryGroup.push_back(elementGroup);
			}
		}
	}
	cout << endl;

	return true;
}

bool ContextPartGraphTrain::extractTrainPartPairs() {

	// get mesh pairs

	cout << "getting mesh pairs" << endl;
	set<vec2i> meshPairSet;
	for (auto &group : mElementGroup) {
		set<int> meshSet;
		for (vec2i element : group) meshSet.insert(element[0]);
		vector<int> meshList(meshSet.begin(), meshSet.end());
		for (int i = 0; i < (int)meshList.size(); i++) {
			for (int j = i + 1; j < (int)meshList.size(); j++) {
				vec2i meshPair(meshList[i], meshList[j]);
				if (meshPair[0] > meshPair[1]) swap(meshPair[0], meshPair[1]);
				meshPairSet.insert(meshPair);
			}
		}
	}
	mTrainMeshPairs.assign(meshPairSet.begin(), meshPairSet.end());

	// match mesh pairs

	string weightsFolder = StyleSynthesisConfig::mData_DataSetRootFolder + "weights/";

	cout << "matching mesh pairs" << endl;

	int numMeshPairs = (int)mTrainMeshPairs.size();
	vector<Eigen::MatrixXd> similarityMatrices(numMeshPairs);

	vector<double> accumNodeSigma(0);
	vector<double> accumEdgeSigma(0);
	int dimNodeSigma = 0;
	int dimEdgeSigma = 0;

	for (int meshPairID = 0; meshPairID < numMeshPairs; meshPairID++) {
		vec2i meshPair = mTrainMeshPairs[meshPairID];

		ContextPartGraphMatch cpgm;
		if (!cpgm.loadWeights(weightsFolder)) return false;
		if (!cpgm.loadGraph((*mpSceneGraphs)[meshPair[0]], (*mpSceneGraphs)[meshPair[1]])) return false;
		if (!cpgm.computeNodeDistance()) return false;
		if (!cpgm.computeNodeSigma()) return false;
		if (!cpgm.evaluateNodeSimilarity()) return false;
		similarityMatrices[meshPairID] = cpgm.mNodeSimilarityMatrix;

		if (accumNodeSigma.empty()) {
			accumNodeSigma = cpgm.mNodeSigma;
			dimNodeSigma = (int)accumNodeSigma.size();
		} else {
			for (int dim = 0; dim < dimNodeSigma; dim++) {
				accumNodeSigma[dim] += cpgm.mNodeSigma[dim];
			}
		}

		if (!cpgm.computeEdgeDistance()) {
			// special handling: huge graph
			cpgm.mEdgeSigma = cpgm.mNodeSigma;
			cout << "Note: huge graph detected" << endl;
			//system("pause");
		} else {
			if (!cpgm.computeEdgeSigma()) return false;
		}

		if (accumEdgeSigma.empty()) {
			accumEdgeSigma = cpgm.mEdgeSigma;
			dimEdgeSigma = (int)accumEdgeSigma.size();
		} else {
			for (int dim = 0; dim < dimEdgeSigma; dim++) {
				accumEdgeSigma[dim] += cpgm.mEdgeSigma[dim];
			}
		}
	}

	// get distance sigma

	if (true) {
		mSceneNodeSigma.assign(dimNodeSigma, 0);
		mSceneEdgeSigma.assign(dimEdgeSigma, 0);
		for (int dim = 0; dim < dimNodeSigma; dim++) {
			mSceneNodeSigma[dim] = accumNodeSigma[dim] / numMeshPairs;
		}
		for (int dim = 0; dim < dimEdgeSigma; dim++) {
			mSceneEdgeSigma[dim] = accumEdgeSigma[dim] / numMeshPairs;
		}
	}

	// extract part pairs

	cout << "extracting part pairs" << endl;
	mTrainPartPairs.assign(numMeshPairs, vector<vec2i>(0));
	for (auto &group : mElementGroup) {
		for (int meshPairID = 0; meshPairID < numMeshPairs; meshPairID++) {
			Eigen::MatrixXd &simMat = similarityMatrices[meshPairID];
			vector<int> srcParts(0);
			vector<int> tgtParts(0);
			for (vec2i element : group) {
				if (element[0] == mTrainMeshPairs[meshPairID][0]) srcParts.push_back(element[1]);
				if (element[0] == mTrainMeshPairs[meshPairID][1]) tgtParts.push_back(element[1]);
			}
			if (srcParts.empty() || tgtParts.empty()) continue;
			Eigen::MatrixXd slicedMat(srcParts.size(), tgtParts.size());
			for (int row = 0; row < (int)srcParts.size(); row++) {
				for (int col = 0; col < (int)tgtParts.size(); col++) {
					slicedMat(row, col) = simMat(srcParts[row], tgtParts[col]);
				}
			}
			while (true) {
				int maxRowID, maxColID;
				double maxSim = slicedMat.maxCoeff(&maxRowID, &maxColID);
				if (maxSim == 0) break;
				mTrainPartPairs[meshPairID].push_back(vec2i(srcParts[maxRowID], tgtParts[maxColID]));
				slicedMat.row(maxRowID).setZero();
				slicedMat.col(maxColID).setZero();
			}
		}
	}

	return true;
}

bool ContextPartGraphTrain::checkNodes(TNode *sourceNode, TNode *targetNode, bool &isElement) {

	double errorFactor = StyleSynthesisConfig::mContext_SymmetryCheckingThreshold;

	vector<int> &srcSegment = (*(sourceNode->mpGraphSegments))[sourceNode->mPartLevelID][sourceNode->mPartSegmentID];
	vector<int> &tgtSegment = (*(targetNode->mpGraphSegments))[targetNode->mPartLevelID][targetNode->mPartSegmentID];

	Eigen::Matrix3Xd matSP, matSN, matTP, matTN;
	if (StyleSynthesisConfig::mContext_SymmetryCheckingByVertex) {
		// compare mesh vertices
		TTriangleMesh srcPartMesh, tgtPartMesh;
		if (!MeshUtil::extractSubMesh(*(sourceNode->mpGraphMesh), srcSegment, srcPartMesh)) return false;
		if (!MeshUtil::extractSubMesh(*(targetNode->mpGraphMesh), tgtSegment, tgtPartMesh)) return false;
		if (!MeshUtil::recomputeNormals(srcPartMesh)) return false;
		if (!MeshUtil::recomputeNormals(tgtPartMesh)) return false;
		if (!SampleUtil::buildMatrices(srcPartMesh, matSP, matSN)) return false;
		if (!SampleUtil::buildMatrices(tgtPartMesh, matTP, matTN)) return false;
	} else {
		// compare sample points
		// pre-process node data (will skip if already pre-processed)
		if (!sourceNode->mNodeDescriptors.preprocess(*(sourceNode->mpGraphMesh), srcSegment)) return false;
		if (!targetNode->mNodeDescriptors.preprocess(*(targetNode->mpGraphMesh), tgtSegment)) return false;
		matSP = sourceNode->mNodeDescriptors.mSampleMatP;
		matSN = sourceNode->mNodeDescriptors.mSampleMatN;
		matTP = targetNode->mNodeDescriptors.mSampleMatP;
		matTN = targetNode->mNodeDescriptors.mSampleMatN;
	}

	// check similarity by ICP

	Eigen::Affine3d xform;
	if (!MatchSimpleICP::runRegularShape(20, matSP, matSN, matTP, matTN, xform)) return false; // UNDONE: param match iterations

	// check scaling

	if (true) {
		Eigen::Matrix3d linearMat = xform.linear();
		Eigen::JacobiSVD<Eigen::Matrix3d> svd(linearMat);
		Eigen::Vector3d scaleFactor = svd.singularValues();
		if (scaleFactor.maxCoeff() > 5.0 || scaleFactor.minCoeff() < 0.2) { // weird scaling
			isElement = false;
			return true;
		}
	}

	// check point to point error

	double errorST;
	Eigen::Matrix3Xd matXSP = xform * matSP;
	if (!MatchSimpleICP::error(matXSP, matTP, errorST)) return false;

	double errorTS;
	Eigen::Matrix3Xd matXTP = xform.inverse() * matTP;
	if (!MatchSimpleICP::error(matXTP, matSP, errorTS)) return false;

	isElement = false;
	if (errorST < cml::sqr(targetNode->mNodeDescriptors.mSizeEpsilon * errorFactor) &&
		errorTS < cml::sqr(sourceNode->mNodeDescriptors.mSizeEpsilon * errorFactor))
	{
		isElement = true;
	}

	return true;
}

bool ContextPartGraphTrain::visualizeElements(string fileName) {

	int numMesh = (int)mpSceneMeshes->size();

	vector<vector<vec3i>> faceColors(numMesh);
	for (int meshID = 0; meshID < numMesh; meshID++) {
		faceColors[meshID].assign((*mpSceneMeshes)[meshID].indices.size(), vec3i(127, 127, 127));
	}

	int numElements = (int)mElementGroup.size();
	for (int elementID = 0; elementID < numElements; elementID++) {

		vec3i color = SegmentUtil::colorMapping(elementID);

		for (vec2i element : mElementGroup[elementID]) {
			int meshID = element[0];
			int nodeID = element[1];
			TNode *node = (*mpSceneGraphs)[meshID].mAllNodes[nodeID];
			vector<int> &segment = (*mpSceneSegments)[meshID][node->mPartLevelID][node->mPartSegmentID];
			for (int faceID : segment) {
				faceColors[meshID][faceID] = color;
			}
		}
	}

	PlyExporter pe;

	vec3 offset(0.0f, 0.0f, 0.0f);
	for (int meshID = 0; meshID < numMesh; meshID++) {
		TTriangleMesh &mesh = (*mpSceneMeshes)[meshID];

		vec3 bbMin, bbMax;
		if (!MeshUtil::computeAABB(mesh, bbMin, bbMax)) return false;
		vec3 bbSize = bbMax - bbMin;
		offset[0] += bbSize[0] * 0.6f;

		if (!pe.addMesh(&mesh.positions, &mesh.normals, &mesh.indices, &faceColors[meshID], offset)) return false;

		offset[0] += bbSize[0] * 0.6f;
	}

	if (!pe.output(fileName)) return false;

	return true;
}

bool ContextPartGraphTrain::visualizeSymmetries(string fileName) {

	int numMesh = (int)mpSceneMeshes->size();

	vector<vector<vec3i>> faceColors(numMesh);
	for (int meshID = 0; meshID < numMesh; meshID++) {
		faceColors[meshID].assign((*mpSceneMeshes)[meshID].indices.size(), vec3i(127, 127, 127));
	}

	int numGroups = (int)mSymmetryGroup.size();
	for (int groupID = 0; groupID < numGroups; groupID++) {

		vec3i color = SegmentUtil::colorMapping(groupID);

		for (vec2i element : mSymmetryGroup[groupID]) {
			int meshID = element[0];
			int nodeID = element[1];
			if (nodeID >= (int)(*mpSceneGraphs)[meshID].mRootNode->mChildren.size()) continue; // skip non-group nodes
			TNode *node = (*mpSceneGraphs)[meshID].mAllNodes[nodeID];
			vector<int> &segment = (*mpSceneSegments)[meshID][node->mPartLevelID][node->mPartSegmentID];
			for (int faceID : segment) {
				faceColors[meshID][faceID] = color;
			}
		}
	}

	PlyExporter pe;

	vec3 offset(0.0f, 0.0f, 0.0f);
	for (int meshID = 0; meshID < numMesh; meshID++) {
		TTriangleMesh &mesh = (*mpSceneMeshes)[meshID];

		vec3 bbMin, bbMax;
		if (!MeshUtil::computeAABB(mesh, bbMin, bbMax)) return false;
		vec3 bbSize = bbMax - bbMin;
		offset[0] += bbSize[0] * 0.6f;

		if (!pe.addMesh(&mesh.positions, &mesh.normals, &mesh.indices, &faceColors[meshID], offset)) return false;

		offset[0] += bbSize[0] * 0.6f;
	}

	if (!pe.output(fileName)) return false;

	return true;
}

bool ContextPartGraphTrain::visualizePartPairs(string fileName) {

	PlyExporter pe;
	vec3 offset(0.0f, 0.0f, 0.0f);

	int numMeshPairs = (int)mTrainMeshPairs.size();
	for (int meshPairID = 0; meshPairID < numMeshPairs; meshPairID++) {
		vec2i meshPair = mTrainMeshPairs[meshPairID];
		TTriangleMesh &srcMesh = (*mpSceneMeshes)[meshPair[0]];
		TTriangleMesh &tgtMesh = (*mpSceneMeshes)[meshPair[1]];

		vector<vec3i> srcFaceColors((*mpSceneMeshes)[meshPair[0]].indices.size(), vec3i(127, 127, 127));
		vector<vec3i> tgtFaceColors((*mpSceneMeshes)[meshPair[1]].indices.size(), vec3i(127, 127, 127));

		for (int partPairID = 0; partPairID < (int)mTrainPartPairs[meshPairID].size(); partPairID++) {
			vec2i partPair = mTrainPartPairs[meshPairID][partPairID];
			TNode *srcNode = (*mpSceneGraphs)[meshPair[0]].mAllNodes[partPair[0]];
			TNode *tgtNode = (*mpSceneGraphs)[meshPair[1]].mAllNodes[partPair[1]];
			vector<int> &srcSegment = (*mpSceneSegments)[meshPair[0]][srcNode->mPartLevelID][srcNode->mPartSegmentID];
			vector<int> &tgtSegment = (*mpSceneSegments)[meshPair[1]][tgtNode->mPartLevelID][tgtNode->mPartSegmentID];
			vec3i color = SegmentUtil::colorMapping(partPairID);
			for (int faceID : srcSegment) srcFaceColors[faceID] = color;
			for (int faceID : tgtSegment) tgtFaceColors[faceID] = color;
		}

		vec3 srcBBMin, srcBBMax;
		vec3 tgtBBMin, tgtBBMax;
		if (!MeshUtil::computeAABB(srcMesh, srcBBMin, srcBBMax)) return false;
		if (!MeshUtil::computeAABB(tgtMesh, tgtBBMin, tgtBBMax)) return false;
		vec3 srcBBSize = srcBBMax - srcBBMin;
		vec3 tgtBBSize = tgtBBMax - tgtBBMin;
		vec3 maxBBSize = srcBBSize;
		maxBBSize.maximize(tgtBBSize);
		offset[0] += maxBBSize[0] * 0.6f;
		vec3 stOffset(0.0f, maxBBSize[1] * 1.2f, 0.0f);

		if (!pe.addMesh(&srcMesh.positions, &srcMesh.normals, &srcMesh.indices, &srcFaceColors, offset)) return false;
		if (!pe.addMesh(&tgtMesh.positions, &tgtMesh.normals, &tgtMesh.indices, &tgtFaceColors, offset+stOffset)) return false;

		offset[0] += maxBBSize[0] * 0.6f;
	}

	if (!pe.output(fileName)) return false;

	return true;
}


bool ContextPartGraphTrain::outputElements(string fileName) {

	ofstream file(fileName);
	if (!file.is_open()) return error("cannot write to file " + fileName);

	int numGroups = (int)mElementGroup.size();
	file << numGroups << endl;

	for (int groupID = 0; groupID < numGroups; groupID++) {
		int numElements = (int)mElementGroup[groupID].size();
		file << numElements << endl;
		for (int elementID = 0; elementID < numElements; elementID++) {
			vec2i index = mElementGroup[groupID][elementID];
			file << index << endl;
		}
	}

	file.close();

	return true;
}

bool ContextPartGraphTrain::outputPartPairs(string fileName) {

	ofstream file(fileName);
	if (!file.is_open()) return error("cannot write to file " + fileName);

	int numMeshPairs = (int)mTrainMeshPairs.size();
	file << numMeshPairs << endl;
	for (int meshPairID = 0; meshPairID < numMeshPairs; meshPairID++) {
		vec2i meshPair = mTrainMeshPairs[meshPairID];
		file << meshPair[0] << " " << meshPair[1];
		for (vec2i partPair : mTrainPartPairs[meshPairID]) {
			file << " " << partPair[0] << " " << partPair[1];
		}
		file << endl;
	}

	file.close();

	return true;
}

bool ContextPartGraphTrain::outputSigma(string nodeSigmaName, string edgeSigmaName) {

	if (!DataUtil::saveValueListASCII(nodeSigmaName, mSceneNodeSigma)) return false;
	if (!DataUtil::saveValueListASCII(edgeSigmaName, mSceneEdgeSigma)) return false;

	return true;
}