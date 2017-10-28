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

#include "ContextPartGraphAssembleResult.h"

#include <iostream>
#include <fstream>
#include <set>

#include "Mesh/MeshUtil.h"

#include "Segment/SegmentGroupApxCvx.h"
#include "Segment/SegmentUtil.h"

#include "Utility/PlyExporter.h"

#include "ContextPartGraphAssemble.h"

#include "Data/DataUtil.h"

using namespace StyleSynthesis;

ContextPartGraphAssembleResult::ContextPartGraphAssembleResult() {

}

ContextPartGraphAssembleResult::~ContextPartGraphAssembleResult() {

}

bool ContextPartGraphAssembleResult::identity(TGraph *graph) {

	mGraph = *graph;
	mMesh = *graph->mRootNode->mpGraphMesh;
	mSegment = *graph->mRootNode->mpGraphSegments;

	int numNodes = (int)graph->mAllNodes.size();

	mNodeMapping.resize(numNodes);
	for (int k = 0; k < numNodes; k++) mNodeMapping[k] = k;
	mNodeTransformation.assign(numNodes, Eigen::Affine3d::Identity());

	mReplaceMapping.assign(numNodes, vec2i(-1, -1));
	mReplaceTransformation.assign(numNodes, Eigen::Affine3d::Identity());

	return true;
}

bool ContextPartGraphAssembleResult::generate(TAssemble *assemble) {

	int numMatchPairs = (int)assemble->mMatchPairs.size();
	int numAddingParts = (int)assemble->mAddingNodeIndices.size();

	TTriangleMesh &sourceMesh = *assemble->mpSourceGraph->mRootNode->mpGraphMesh;
	TTriangleMesh &targetMesh = *assemble->mpTargetGraph->mRootNode->mpGraphMesh;
	int numSourceMeshFaces = (int)sourceMesh.indices.size();
	int numTargetMeshFaces = (int)targetMesh.indices.size();

	////////////////////////////////// assemble mesh //////////////////////////////////

	// NOTE: it's necessary to use separate face mappings for source / adding
	// one source node may be used to "replace several target nodes" / "be added to target shape" under different transformations
	vector<vector<int>> sourceFaceMapping(numMatchPairs, vector<int>(numSourceMeshFaces, -1)); // face ID in assemble mesh : # of faces in source mesh : # of matches
	vector<vector<int>> addingFaceMapping(numAddingParts, vector<int>(numSourceMeshFaces, -1)); // face ID in assemble mesh : # of faces in source mesh : # of adding operations
	vector<int> targetFaceMapping(numTargetMeshFaces, -1); // face ID in assemble mesh : # of faces in target mesh

	int numTargetNodes = (int)assemble->mpTargetGraph->mAllNodes.size();
	vector<bool> targetNodeActiveFlags(numTargetNodes, false);

	if (true) {

		PlyExporter collector;
		int numAssembleFaces = 0;

		// mark active nodes

		for (int targetID : assemble->mTargetWorkingNodes) {
			targetNodeActiveFlags[targetID] = true;
		}
		for (int matchID = 0; matchID < numMatchPairs; matchID++) {
			int targetID = assemble->mMatchPairs[matchID][1];
			targetNodeActiveFlags[targetID] = false;
		}

		// process source replacement mesh

		for (int matchID = 0; matchID < numMatchPairs; matchID++) {
			int sourceID = assemble->mMatchPairs[matchID][0];
			int targetID = assemble->mMatchPairs[matchID][1];
			TNode *sourceNode = assemble->mpSourceGraph->mAllNodes[sourceID];
			vector<int> &sourceSegment = (*sourceNode->mpGraphSegments)[sourceNode->mPartLevelID][sourceNode->mPartSegmentID];

			// compute face mapping
			for (int faceID : sourceSegment) {
				sourceFaceMapping[matchID][faceID] = numAssembleFaces;
				numAssembleFaces++;
			}

			// transform replacement mesh
			TTriangleMesh replaceMesh;
			if (!MeshUtil::extractSubMesh(sourceMesh, sourceSegment, replaceMesh)) return false;
			Eigen::Affine3d matchXform = assemble->mMatchTransformations[matchID];
			Eigen::Affine3d nodeXform = assemble->mNodeTransformations[targetID];
			Eigen::Affine3d transform = nodeXform * matchXform;
			Eigen::Matrix3d rotation = transform.rotation();
			for (vec3 &position : replaceMesh.positions) {
				Eigen::Vector3d vec(vec3d(position).data());
				vec = transform * vec;
				position = vec3d(vec[0], vec[1], vec[2]);
			}
			for (vec3 &normal : replaceMesh.normals) {
				Eigen::Vector3d vec(vec3d(normal).data());
				vec = rotation * vec;
				normal = vec3d(vec[0], vec[1], vec[2]);
			}

			// collect replacement mesh
			if (!collector.addMesh(&replaceMesh.indices, &replaceMesh.positions, &replaceMesh.normals)) return false;
		}

		// process target active node mesh

		for (int targetID : assemble->mTargetWorkingNodes) {
			if (!targetNodeActiveFlags[targetID]) continue; // process active nodes only

			TNode *targetNode = assemble->mpTargetGraph->mAllNodes[targetID];
			vector<int> &targetSegment = (*targetNode->mpGraphSegments)[targetNode->mPartLevelID][targetNode->mPartSegmentID];

			// compute face mapping
			for (int faceID : targetSegment) {
				targetFaceMapping[faceID] = numAssembleFaces;
				numAssembleFaces++;
			}

			// transform node mesh
			TTriangleMesh nodeMesh;
			if (!MeshUtil::extractSubMesh(targetMesh, targetSegment, nodeMesh)) return false;
			Eigen::Affine3d transform = assemble->mNodeTransformations[targetID];
			Eigen::Matrix3d rotation = transform.rotation();
			for (vec3 &position : nodeMesh.positions) {
				Eigen::Vector3d vec(vec3d(position).data());
				vec = transform * vec;
				position = vec3d(vec[0], vec[1], vec[2]);
			}
			for (vec3 &normal : nodeMesh.normals) {
				Eigen::Vector3d vec(vec3d(normal).data());
				vec = rotation * vec;
				normal = vec3d(vec[0], vec[1], vec[2]);
			}

			// collect node mesh
			if (!collector.addMesh(&nodeMesh.indices, &nodeMesh.positions, &nodeMesh.normals)) return false;
		}

		// process source adding mesh

		for (int addingID = 0; addingID < numAddingParts; addingID++) {
			int sourceID = assemble->mAddingNodeIndices[addingID];
			TNode *sourceNode = assemble->mpSourceGraph->mAllNodes[sourceID];
			vector<int> &sourceSegment = (*sourceNode->mpGraphSegments)[sourceNode->mPartLevelID][sourceNode->mPartSegmentID];

			// compute face mapping
			for (int faceID : sourceSegment) {
				addingFaceMapping[addingID][faceID] = numAssembleFaces;
				numAssembleFaces++;
			}

			// transfer adding mesh
			TTriangleMesh addingMesh;
			if (!MeshUtil::extractSubMesh(sourceMesh, sourceSegment, addingMesh)) return false;
			Eigen::Affine3d transform = assemble->mAddingTransformations[addingID];
			Eigen::Matrix3d rotation = transform.rotation();
			for (vec3 &position : addingMesh.positions) {
				Eigen::Vector3d vec(vec3d(position).data());
				vec = transform * vec;
				position = vec3d(vec[0], vec[1], vec[2]);
			}
			for (vec3 &normal : addingMesh.normals) {
				Eigen::Vector3d vec(vec3d(normal).data());
				vec = rotation * vec;
				normal = vec3d(vec[0], vec[1], vec[2]);
			}

			// collect adding mesh
			if (!collector.addMesh(&addingMesh.indices, &addingMesh.positions, &addingMesh.normals)) return false;
		}

		// assemble mesh

		if (true) {
			mMesh.positions = collector.mVertices;
			mMesh.normals = collector.mNormals;
			mMesh.indices = collector.mFaceIndices;
			mMesh.amount = (int)mMesh.positions.size();
			vector<int> vertexIndices;
			MeshUtil::removeDuplicateVertices(mMesh, mMesh, vertexIndices);
		}
	}

	////////////////////////////////// construct hierarchy //////////////////////////////////

	typedef ContextPartGraphAssembleUtil::TNodeType TNodeType;

	vector<vector<int>> assembleSegmentNode; // target node ID (-1 for replaced/modified/added node) : # of nodes in this level : # of levels
	vector<vector<int>> assembleSegmentMatch; // match pair ID (-1 for non-replaced node) : # of nodes in this level : # of levels
	vector<vector<int>> assembleSegmentAdding; // adding operation ID (-1 for non-added node) : # of nodes in this level : # of levels

	if (true) {
		int numTargetNodes = (int)assemble->mpTargetGraph->mAllNodes.size();

		// compute node level

		vector<int> nodeLevel(numTargetNodes, -1);
		int numMaxLevels = 3;
		for (int nodeID = 0; nodeID < numTargetNodes; nodeID++) {
			int parentID = assemble->mpTargetGraph->mAllNodes[nodeID]->mParent->mUID;
			if (parentID < 0) nodeLevel[nodeID] = 0;
			else nodeLevel[nodeID] = nodeLevel[parentID] + 1;
			numMaxLevels = max(numMaxLevels, nodeLevel[nodeID] + 1);
		}

		mSegment.assign(numMaxLevels, vector<vector<int>>(0));
		assembleSegmentNode.assign(numMaxLevels, vector<int>(0)); // make sure it's synchronized with mSegment
		assembleSegmentMatch.assign(numMaxLevels, vector<int>(0)); // make sure it's synchronized with mSegment
		assembleSegmentAdding.assign(numMaxLevels, vector<int>(0)); // make sure it's synchronized with mSegment

		// push replaced segments from source nodes into hierarchy (should be placed before target node segments)

		vector<vector<int>> nodeNewSegment(numTargetNodes, vector<int>(0));

		for (int matchID = 0; matchID < numMatchPairs; matchID++) {
			int sourceID = assemble->mMatchPairs[matchID][0];
			int targetID = assemble->mMatchPairs[matchID][1];
			int levelID = nodeLevel[targetID];
			vector<vec2i> queue(1, vec2i(sourceID, levelID));
			int head = 0;
			while (head < (int)queue.size()) {
				int nodeID = queue[head][0];
				int levelID = queue[head][1];
				TNode *node = assemble->mpSourceGraph->mAllNodes[nodeID];

				vector<int> &segment = (*node->mpGraphSegments)[node->mPartLevelID][node->mPartSegmentID];
				vector<int> newSegment(0);
				for (int faceID : segment) newSegment.push_back(sourceFaceMapping[matchID][faceID]);
				if (!newSegment.empty()) {
					mSegment[levelID].push_back(newSegment);
					assembleSegmentNode[levelID].push_back(-1);
					assembleSegmentMatch[levelID].push_back(head == 0 ? matchID : -1);
					assembleSegmentAdding[levelID].push_back(-1);

					if (head == 0) {
						nodeNewSegment[targetID] = newSegment;
					}
				}

				int nextLevelID = levelID + 1;
				if (nextLevelID < numMaxLevels) {
					for (TNode *child : node->mChildren) {
						queue.push_back(vec2i(child->mUID, nextLevelID));
					}
				}

				head++;
			}
		}

		// get target node segments

		for (int nodeID = numTargetNodes - 1; nodeID >= 0; nodeID--) {
			int levelID = nodeLevel[nodeID];
			vector<int> &newSegment = nodeNewSegment[nodeID];

			TNodeType nodeType = assemble->mTargetNodeTypes[nodeID];
			if (nodeType == TNodeType::TNodeType_ACTIVE || nodeType == TNodeType::TNodeType_INACTIVE) {
				// keep segment as is
				TNode *node = assemble->mpTargetGraph->mAllNodes[nodeID];
				vector<int> &segment = (*node->mpGraphSegments)[node->mPartLevelID][node->mPartSegmentID];
				for (int faceID : segment) newSegment.push_back(targetFaceMapping[faceID]);
				if (!newSegment.empty()) {
					mSegment[levelID].push_back(newSegment);
					assembleSegmentNode[levelID].push_back(nodeID);
					assembleSegmentMatch[levelID].push_back(-1);
					assembleSegmentAdding[levelID].push_back(-1);
				}
			} else if (nodeType == TNodeType::TNodeType_MODIFIED) {
				// use all children's segments
				for (TNode *child : assemble->mpTargetGraph->mAllNodes[nodeID]->mChildren) {
					int childID = child->mUID;
					vector<int> &childSegment = nodeNewSegment[childID];
					newSegment.insert(newSegment.end(), childSegment.begin(), childSegment.end());
				}
				if (!newSegment.empty()) {
					mSegment[levelID].push_back(newSegment);
					assembleSegmentNode[levelID].push_back(-1);
					assembleSegmentMatch[levelID].push_back(-1);
					assembleSegmentAdding[levelID].push_back(-1);
				}
			} else {
				// skip other types of nodes (segments in replaced node should already be ready)
			}
		}

		// push added segments from source nodes into hierarchy

		for (int addingID = 0; addingID < numAddingParts; addingID++) {
			int sourceID = assemble->mAddingNodeIndices[addingID];
			int levelID = 0; // add to top level (just below root)
			vector<vec2i> queue(1, vec2i(sourceID, levelID));
			int head = 0;
			while (head < (int)queue.size()) {
				int nodeID = queue[head][0];
				int levelID = queue[head][1];
				TNode *node = assemble->mpSourceGraph->mAllNodes[nodeID];

				vector<int> &segment = (*node->mpGraphSegments)[node->mPartLevelID][node->mPartSegmentID];
				vector<int> newSegment(0);
				for (int faceID : segment) newSegment.push_back(addingFaceMapping[addingID][faceID]);
				if (!newSegment.empty()) {
					mSegment[levelID].push_back(newSegment);
					assembleSegmentNode[levelID].push_back(-1);
					assembleSegmentMatch[levelID].push_back(-1);
					assembleSegmentAdding[levelID].push_back(head == 0 ? addingID : -1);
				}

				int nextLevelID = levelID + 1;
				if (nextLevelID < numMaxLevels) {
					for (TNode *child : node->mChildren) {
						queue.push_back(vec2i(child->mUID, nextLevelID));
					}
				}

				head++;
			}
		}
	}

	////////////////////////////////// construct graph //////////////////////////////////

	cout << "segment level: ";
	for (auto &segment : mSegment) {
		cout << segment.size() << " ";
	}
	cout << endl;
	if (!mGraph.buildGraphHierarchy(*assemble->mpNodeGenerator, &mMesh, &mSegment)) return false;
	if (!mGraph.buildGraphDescriptor()) return false;
	if (!mGraph.buildGraphContext()) return false;

	int numNewNodes = (int)mGraph.mAllNodes.size();

	////////////////////////////////// build node mapping //////////////////////////////////

	if (true) {
		mNodeMapping.assign(numTargetNodes, -1);
		for (int nodeID = 0; nodeID < numNewNodes; nodeID++) {
			TNode *node = mGraph.mAllNodes[nodeID];
			int mappedNodeID = assembleSegmentNode[node->mPartLevelID][node->mPartSegmentID];
			if (mappedNodeID >= 0) mNodeMapping[mappedNodeID] = nodeID;
		}

		//cout << "Node mapping:";
		//for (int nodeID = 0; nodeID < numTargetNodes; nodeID++) {
		//	cout << " " << mNodeMapping[nodeID];
		//}
		//cout << endl;
	}

	////////////////////////////////// build node transformation //////////////////////////////////

	if (true) {
		vector<Eigen::Affine3d> targetNodeTransformation = assemble->mNodeTransformations;
		vector<bool> targetNodeValidFlags = targetNodeActiveFlags;
		for (int nodeID = 0; nodeID < numTargetNodes; nodeID++) {
			int parentID = assemble->mpTargetGraph->mAllNodes[nodeID]->mParent->mUID;
			if (parentID >= 0 && targetNodeValidFlags[parentID]) {
				targetNodeTransformation[nodeID] = targetNodeTransformation[parentID]; // copy parent's transformation
				targetNodeValidFlags[nodeID] = true;
			}
		}
		mNodeTransformation.assign(numTargetNodes, Eigen::Affine3d::Identity());
		for (int nodeID = 0; nodeID < numTargetNodes; nodeID++) {
			if (targetNodeValidFlags[nodeID]) {
				mNodeTransformation[nodeID] = targetNodeTransformation[nodeID];
			}
		}
	}

	////////////////////////////////// build replacement info //////////////////////////////////

	if (true) {
		mReplaceMapping.assign(numNewNodes, vec2i(-1, -1));
		mReplaceTransformation.assign(numNewNodes, Eigen::Affine3d::Identity());
		for (int nodeID = 0; nodeID < numNewNodes; nodeID++) {
			TNode *node = mGraph.mAllNodes[nodeID];
			int mappedMatchID = assembleSegmentMatch[node->mPartLevelID][node->mPartSegmentID];
			if (mappedMatchID >= 0) {
				vec2i matchPair = assemble->mMatchPairs[mappedMatchID];
				Eigen::Affine3d matchXform = assemble->mMatchTransformations[mappedMatchID];
				Eigen::Affine3d nodeXform = assemble->mNodeTransformations[matchPair[1]];
				mReplaceMapping[nodeID] = matchPair;
				mReplaceTransformation[nodeID] = nodeXform * matchXform;
			}
			int mappedAddingID = assembleSegmentAdding[node->mPartLevelID][node->mPartSegmentID];
			if (mappedAddingID >= 0) {
				int sourceID = assemble->mAddingNodeIndices[mappedAddingID];
				Eigen::Affine3d addingXform = assemble->mAddingTransformations[mappedAddingID];
				mReplaceMapping[nodeID] = vec2i(sourceID, -1);
				mReplaceTransformation[nodeID] = addingXform;
			}
		}
	}

	return true;
}

bool ContextPartGraphAssembleResult::chain(ContextPartGraphAssembleResult *previousResult) {

	// chain replacement info (should before chaining nodes info)

	vector<vec2i> chainedReplaceMapping = mReplaceMapping;
	vector<Eigen::Affine3d> chainedReplaceTransformation = mReplaceTransformation;
	if (true) {
		vector<int> targetNodeMapping(mNodeMapping.size(), -1); // original target node ID : # of current target nodes
		for (int origTargetID = 0; origTargetID < (int)previousResult->mNodeMapping.size(); origTargetID++) {
			int targetID = previousResult->mNodeMapping[origTargetID];
			if (targetID >= 0) targetNodeMapping[targetID] = origTargetID;
		}
		for (vec2i &replacePair : chainedReplaceMapping) {
			if (replacePair[1] >= 0) replacePair[1] = targetNodeMapping[replacePair[1]];
		}

		for (int nodeID = 0; nodeID < (int)previousResult->mReplaceMapping.size(); nodeID++) {
			vec2i replacePair = previousResult->mReplaceMapping[nodeID];
			Eigen::Affine3d replaceXform = previousResult->mReplaceTransformation[nodeID];
			if (replacePair[0] >= 0) {
				int newNodeID = mNodeMapping[nodeID];
				chainedReplaceMapping[newNodeID] = replacePair;
				chainedReplaceTransformation[newNodeID] = mNodeTransformation[nodeID] * replaceXform;
			}
		}
	}

	// chain nodes info

	vector<int> chainedNodeMapping = previousResult->mNodeMapping;
	vector<Eigen::Affine3d> chainedNodeTransformation = previousResult->mNodeTransformation;
	if (true) {
		int numNodes = (int)chainedNodeMapping.size();
		for (int nodeID = 0; nodeID < numNodes; nodeID++) {
			int chainedID = chainedNodeMapping[nodeID];
			Eigen::Affine3d chainedXform = chainedNodeTransformation[nodeID];
			if (chainedID >= 0) {
				chainedNodeMapping[nodeID] = mNodeMapping[chainedID];
				chainedNodeTransformation[nodeID] = mNodeTransformation[chainedID] * chainedXform;
			}
		}
	}

	// update data

	mNodeMapping.swap(chainedNodeMapping);
	mNodeTransformation.swap(chainedNodeTransformation);
	mReplaceMapping.swap(chainedReplaceMapping);
	mReplaceTransformation.swap(chainedReplaceTransformation);

	return true;
}

bool ContextPartGraphAssembleResult::visualize(string fileName) {

	// find replaced nodes

	int numNodes = (int)mGraph.mAllNodes.size();
	vector<bool> nodeReplacedFlags(numNodes, true);
	for (int nodeID : mNodeMapping) {
		if (nodeID >= 0) {
			nodeReplacedFlags[nodeID] = false;
			int parentID = mGraph.mAllNodes[nodeID]->mParent->mUID;
			while (parentID >= 0) {
				nodeReplacedFlags[parentID] = false;
				parentID = mGraph.mAllNodes[parentID]->mParent->mUID;
			}
		}
	}
	for (int nodeID = numNodes - 1; nodeID >= 0; nodeID--) {
		if (!nodeReplacedFlags[nodeID]) continue;
		int parentID = mGraph.mAllNodes[nodeID]->mParent->mUID;
		if (parentID >= 0 && nodeReplacedFlags[parentID]) nodeReplacedFlags[nodeID] = false;
	}

	// mark face colors

	int numFaces = (int)mMesh.indices.size();
	vector<vec3i> faceColors(numFaces, vec3i(127, 127, 127));
	for (int nodeID = 0; nodeID < numNodes; nodeID++) {
		if (!nodeReplacedFlags[nodeID]) continue;
		TNode *node = mGraph.mAllNodes[nodeID];
		vector<int> &nodeSegment = mSegment[node->mPartLevelID][node->mPartSegmentID];
		vec3i nodeColor = SegmentUtil::colorMapping(nodeID);
		for (int faceID : nodeSegment) {
			faceColors[faceID] = nodeColor;
		}
	}

	// visualize assemble mesh

	PlyExporter pe;
	if (!pe.addMesh(&mMesh.positions, &mMesh.normals, &mMesh.indices, &faceColors)) return false;
	if (!pe.output(fileName)) return false;

	return true;
}

bool ContextPartGraphAssembleResult::visualize(string fileName, ContextPartGraphAssembleResult *prevSolution) {
	
	int numNodes = (int)mGraph.mAllNodes.size();
	vector<bool> nodeReplacedFlags(numNodes, false);

	// find replaced nodes

	set<int> prevRepSet;
	for (vec2i rm : prevSolution->mReplaceMapping) {
		if (rm[1] >= 0) prevRepSet.insert(rm[1]);
	}
	for (int nodeID = 0; nodeID < numNodes; nodeID++) {
		vec2i rm = mReplaceMapping[nodeID];
		if (rm[0] >= 0 && rm[1] < 0) {
			nodeReplacedFlags[nodeID] = true;
		}
		if (rm[1] >= 0 && prevRepSet.find(rm[1]) == prevRepSet.end()) {
			nodeReplacedFlags[nodeID] = true;
		}
	}

	/*

	for (int nodeID : mNodeMapping) {
		if (nodeID >= 0) {
			nodeReplacedFlags[nodeID] = false;
			int parentID = mGraph.mAllNodes[nodeID]->mParent->mUID;
			while (parentID >= 0) {
				nodeReplacedFlags[parentID] = false;
				parentID = mGraph.mAllNodes[parentID]->mParent->mUID;
			}
		}
	}
	for (int nodeID = numNodes - 1; nodeID >= 0; nodeID--) {
		if (!nodeReplacedFlags[nodeID]) continue;
		int parentID = mGraph.mAllNodes[nodeID]->mParent->mUID;
		if (parentID >= 0 && nodeReplacedFlags[parentID]) nodeReplacedFlags[nodeID] = false;
	}
	*/

	// mark face colors

	int numFaces = (int)mMesh.indices.size();
	//vector<vec3i> faceColors(numFaces, vec3i(127, 127, 127));
	vector<vec3i> faceColors(numFaces, vec3i(255, 255, 255));
	for (int nodeID = 0; nodeID < numNodes; nodeID++) {
		if (!nodeReplacedFlags[nodeID]) continue;
		TNode *node = mGraph.mAllNodes[nodeID];
		vector<int> &nodeSegment = mSegment[node->mPartLevelID][node->mPartSegmentID];
		//vec3i nodeColor = SegmentUtil::colorMapping(nodeID);
		vec3i nodeColor(128, 255, 128);
		for (int faceID : nodeSegment) {
			faceColors[faceID] = nodeColor;
		}
	}

	// visualize assemble mesh

	PlyExporter pe;
	if (!pe.addMesh(&mMesh.positions, &mMesh.normals, &mMesh.indices, &faceColors)) return false;
	if (!pe.output(fileName)) return false;

	return true;
}

bool ContextPartGraphAssembleResult::saveData(string folderName) {

	string meshName = folderName + "mesh.ply";
	string segmentName = folderName + "segment.txt";
	string graphHName = folderName + "graph-hierarchy.txt";
	string graphDName = folderName + "graph-descriptor.txt";
	string graphCName = folderName + "graph-context.txt";
	string nodeMappingName = folderName + "node-mapping.txt";
	string nodeTransformationName = folderName + "node-transformation.txt";
	string replaceMappingName = folderName + "replace-mapping.txt";
	string replaceTransformationName = folderName + "replace-transformation.txt";

	if (!MeshUtil::saveMesh(meshName, mMesh)) return false;
	if (!SegmentGroupApxCvx::saveSegments(segmentName, mSegment)) return false;
	if (!mGraph.saveGraphHierarchy(graphHName)) return false;
	if (!mGraph.saveGraphDescriptor(graphDName)) return false;
	if (!mGraph.saveGraphContext(graphCName)) return false;

	if (!DataUtil::saveIndexListASCII(nodeMappingName, mNodeMapping)) return false;
	if (true) {
		int numNodes = (int)mNodeMapping.size();
		int numValidNodes = 0;
		for (int mapping : mNodeMapping) {
			if (mapping >= 0) numValidNodes++;
		}
		Eigen::MatrixXd allXformMat(4, 4 * numValidNodes);
		int validID = 0;
		for (int nodeID = 0; nodeID < numNodes; nodeID++) {
			if (mNodeMapping[nodeID] < 0) continue;
			allXformMat.middleCols(validID * 4, 4) = mNodeTransformation[nodeID].matrix();
			validID++;
		}
		if (!DataUtil::saveMatrixASCII(nodeTransformationName, allXformMat)) return false;
	}

	if (!DataUtil::savePairListASCII(replaceMappingName, mReplaceMapping)) return false;
	if (true) {
		int numNodes = (int)mReplaceMapping.size();
		int numReplacedNodes = 0;
		for (vec2i index : mReplaceMapping) {
			if (index[0] >= 0) numReplacedNodes++;
		}
		Eigen::MatrixXd allXformMat(4, 4 * numReplacedNodes);
		int replaceID = 0;
		for (int nodeID = 0; nodeID < numNodes; nodeID++) {
			if (mReplaceMapping[nodeID][0] < 0) continue;
			allXformMat.middleCols(replaceID * 4, 4) = mReplaceTransformation[nodeID].matrix();
			replaceID++;
		}
		if (!DataUtil::saveMatrixASCII(replaceTransformationName, allXformMat)) return false;
	}

	return true;
}

bool ContextPartGraphAssembleResult::loadData(string folderName, TNodeGen &nodeGen) {

	string meshName = folderName + "mesh.ply";
	string segmentName = folderName + "segment.txt";
	string graphHName = folderName + "graph-hierarchy.txt";
	string graphDName = folderName + "graph-descriptor.txt";
	string graphCName = folderName + "graph-context.txt";
	string nodeMappingName = folderName + "node-mapping.txt";
	string nodeTransformationName = folderName + "node-transformation.txt";
	string replaceMappingName = folderName + "replace-mapping.txt";
	string replaceTransformationName = folderName + "replace-transformation.txt";

	if (!MeshUtil::loadMesh(meshName, mMesh)) return false;
	if (!SegmentGroupApxCvx::loadSegments(segmentName, mSegment)) return false;
	if (!mGraph.loadGraphHierarchy(graphHName, nodeGen, &mMesh, &mSegment)) return false;
	if (!mGraph.loadGraphDescriptor(graphDName)) return false;
	if (!mGraph.loadGraphContext(graphCName)) return false;

	if (!DataUtil::loadIndexListASCII(nodeMappingName, mNodeMapping)) return false;
	if (true) {
		Eigen::MatrixXd allXformMat;
		if (!DataUtil::loadMatrixASCII(nodeTransformationName, allXformMat)) return false;
		int numNodes = (int)mNodeMapping.size();
		mNodeTransformation.assign(numNodes, Eigen::Affine3d::Identity());
		int validID = 0;
		for (int nodeID = 0; nodeID < numNodes; nodeID++) {
			if (mNodeMapping[nodeID] < 0) continue;
			mNodeTransformation[nodeID].matrix() = allXformMat.middleCols(validID * 4, 4);
			validID++;
			if (validID * 4 > (int)allXformMat.cols()) {
				return error("incorrect node transformation file");
			}
		}
	}

	if (!DataUtil::loadPairListASCII(replaceMappingName, mReplaceMapping)) return false;
	if (true) {
		Eigen::MatrixXd allXformMat;
		if (!DataUtil::loadMatrixASCII(replaceTransformationName, allXformMat)) return false;
		int numNodes = (int)mReplaceMapping.size();
		mReplaceTransformation.assign(numNodes, Eigen::Affine3d::Identity());
		int validID = 0;
		for (int nodeID = 0; nodeID < numNodes; nodeID++) {
			if (mReplaceMapping[nodeID][0] < 0) continue;
			mReplaceTransformation[nodeID].matrix() = allXformMat.middleCols(validID * 4, 4);
			validID++;
			if (validID * 4 > (int)allXformMat.cols()) {
				return error("incorrect replace transformation file");
			}
		}
	}

	return true;
}