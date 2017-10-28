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

#include "ContextPartGraph.h"

#include <iostream>
#include <fstream>

#include "Match/MatchRigidICP.h"

#include "Mesh/MeshUtil.h"
#include "Sample/SampleUtil.h"
#include "Segment/SegmentUtil.h"

#include "Sample/SampleSimplePoissonDisk.h"
#include "Segment/SegmentMeshApxCvx.h"

#include "Utility/PlyExporter.h"

#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

ContextPartGraph::ContextPartGraph() {
}

ContextPartGraph::~ContextPartGraph() {
}

bool ContextPartGraph::buildGraphHierarchy(
	TNodeGen &nodeGenerator,
	TTriangleMesh *mesh,
	vector<vector<vector<int>>> *segments)
{

	cout << "Building graph hierarchy..." << endl;

	mRootNode = nodeGenerator.generateNode();
	mRootNode->clearNode();
	mRootNode->mpGraphMesh = mesh;
	mRootNode->mpGraphSegments = segments;
	mRootNode->mUID = -1;

	mAllNodes.clear();
	
	int numFaces = (int)mesh->indices.size();
	vector<int> currentLevelFaceMap(numFaces, -1); // node ID : # of faces
	vector<int> lastLevelFaceMap(numFaces, -1); // initially the parent node is the root node
	vector<int> nodeFaceAmount; // face amount (used for detecting identical segment) : # of nodes

	int numNodes = 0;
	int numLevels = min(3, (int)segments->size()); // UNDONE: param number of levels used for building graph
	for (int levelID = 0; levelID < numLevels; levelID++) {

		vector<vector<int>> &currentLevelSegments = (*segments)[levelID];
		int numSegments = (int)currentLevelSegments.size();

		for (int segmentID = 0; segmentID < numSegments; segmentID++) {

			vector<int> &currentSegment = currentLevelSegments[segmentID];
			TTriangleMesh segmentMesh;
			if (!MeshUtil::extractSubMesh(*mesh, currentSegment, segmentMesh)) return false;

			// check flat segment & small segment
			if (levelID > 0) {// HACK: always keep all segments in first level (part group level)
				double volume;
				if (!MeshUtil::computeVolume(segmentMesh, volume)) return false;
				if (volume < 1e-6) continue; // UNDONE: param flat segment volume threshold

				double sampleRadius = StyleSynthesisConfig::mAssemble_SlotsSamplingRadius;
				TSampleSet segmentSamples;
				SampleSimplePoissonDisk sspd(&segmentMesh);
				if (!sspd.runSampling(sampleRadius)) return false;
				if (!sspd.exportSample(segmentSamples)) return false;
				if (segmentSamples.amount < 20) continue; // UNDONE: param small segment sampling amount threshold
			}

			int parentNodeID = lastLevelFaceMap[currentSegment[0]];
			TNode *parentNode = parentNodeID >= 0 ? mAllNodes[parentNodeID] : mRootNode;

			// check identical segment
			if (parentNodeID >= 0 && nodeFaceAmount[parentNodeID] == (int)currentSegment.size()) {
				continue;
			}

			// add node
			TNode *node = nodeGenerator.generateNode();
			node->mpGraphMesh = mesh;
			node->mpGraphSegments = segments;
			node->mPartLevelID = levelID;
			node->mPartSegmentID = segmentID;
			node->mUID = numNodes; // will be re-assigned later
			mAllNodes.push_back(node);
			nodeFaceAmount.push_back((int)currentSegment.size());

			// update face map
			for (int faceID : currentSegment) {
				currentLevelFaceMap[faceID] = numNodes; // current node ID in mAllNodes
			}

			// link to parent node
			node->mParent = parentNode;
			parentNode->mChildren.push_back(node);

			numNodes++;
		}

		// update face map

		lastLevelFaceMap = currentLevelFaceMap;
	}

	// prune graph

	for (TNode* node : mAllNodes) {
		int numChildren = (int)node->mChildren.size();
		if (numChildren == 1) { // collapse single bridge
			TNode *child = node->mChildren[0];
			for (TNode *grandchild : child->mChildren) grandchild->mParent = node;
			node->mChildren.swap(child->mChildren);
			child->mParent = 0;
			child->mChildren.clear();
		}
	}

	// re-organize mAllNodes to be strictly ordered by layer (nodes in upper layer always precede nodes in lower layer)

	if (true) {
		vector<vector<TNode*>> graphLevels;
		if (!extractGraphLevels(graphLevels)) return false;

		vector<TNode*> newAllNodes(0);
		int levelID = 0;
		for (auto &level : graphLevels) {
			for (TNode *node : level) {
				node->mUID = (int)newAllNodes.size();
				newAllNodes.push_back(node);
			}
			cout << "Level " << levelID << " : " << level.size() << " nodes" << endl;
			levelID++;
		}

		mAllNodes.swap(newAllNodes);
	}

	return true;
}

bool ContextPartGraph::buildGraphDescriptor(int mode) {

	cout << "Building graph descriptor..." << endl;

	int numNodes = (int)mAllNodes.size();

#pragma omp parallel for
	for (int nodeID = 0; nodeID < numNodes; nodeID++) {
		if (!mAllNodes[nodeID]->computeDescriptor(mode)) error("compute node descriptor");
	}

	return true;
}

bool ContextPartGraph::buildGraphContext(int mode) {

	cout << "Building graph context..." << endl;

	for (TNode *node : mAllNodes) {
		if (!node->mNodeDescriptors.preprocess(*mRootNode->mpGraphMesh, (*mRootNode->mpGraphSegments)[node->mPartLevelID][node->mPartSegmentID])) return false;
	}

	int numAllNodes = (int)mAllNodes.size();

	vector<vector<TNode*>> graphLevels;
	if (!extractGraphLevels(graphLevels)) return false;
	int numLevels = (int)graphLevels.size();

	for (int level = 0; level < numLevels; level++) {

		auto &nodes = graphLevels[level];
		int numNodes = (int)nodes.size();

		// detect symmetry cliques (only within the same level)

		vector<int> cliqueMap(numNodes);
		vector<vector<int>> cliques(numNodes);
		for (int nodeID = 0; nodeID < numNodes; nodeID++) {
			cliqueMap[nodeID] = nodeID;
			cliques[nodeID].assign(1, nodeID);
		}

		for (int nodeID = 0; nodeID < numNodes; nodeID++) {
			cout << "\rDetecting symmetry " << (nodeID + 1) << " / " << numNodes << "         ";
			int nodeCliqueID = cliqueMap[nodeID];
			auto &nodeClique = cliques[nodeCliqueID];
			for (int otherID = nodeID + 1; otherID < numNodes; otherID++) {
				int otherCliqueID = cliqueMap[otherID];
				auto &otherClique = cliques[otherCliqueID];
				if (nodeCliqueID == otherCliqueID) continue;
				if (TNode::detectSymmetry(nodes[nodeID], nodes[otherID])) {
					nodeClique.insert(nodeClique.end(), otherClique.begin(), otherClique.end());
					for (int id : otherClique) cliqueMap[id] = nodeCliqueID;
					otherClique.clear();
				}
			}			
		}
		cout << endl;

		for (int nodeID = 0; nodeID < numNodes; nodeID++) {
			auto &clique = cliques[nodeID];
			int numNodesInClique = (int)clique.size();
			for (int id1 = 0; id1 < numNodesInClique - 1; id1++) {
				for (int id2 = id1 + 1; id2 < numNodesInClique; id2++) {
					nodes[clique[id1]]->mSymmetry.push_back(nodes[clique[id2]]);
					nodes[clique[id2]]->mSymmetry.push_back(nodes[clique[id1]]);
				}
			}

			if (level == 0 && numNodesInClique > 1) {
				cout << "Symmetry clique: ";
				for (int cid = 0; cid < (int)clique.size(); cid++) {
					cout << nodes[clique[cid]]->mUID << " ";
				}
				cout << endl;
				/*
				if (true) {
					PlyExporter pe;
					for (int cid = 0; cid < (int)clique.size(); cid++) {
						auto &mesh = nodes[clique[cid]]->mNodeDescriptors.mNodeMesh;
						if (!pe.addMesh(&mesh.indices, &mesh.positions, &mesh.normals)) return false;
					}
					if (!pe.output("symmetry.ply")) return false;
					system("pause");
				}
				*/
			}
		}
	}

	for (int nodeID = 0; nodeID < numAllNodes; nodeID++) {

		cout << "\rDetecting other context " << (nodeID + 1) << " / " << numAllNodes << "         ";
		
		// mark valid nodes (all other nodes except ancestors/descendants)
		vector<bool> validFlags(numAllNodes, true);
		validFlags[nodeID] = false;
		for (int otherID = nodeID + 1; otherID < numAllNodes; otherID++) {
			int parentID = mAllNodes[otherID]->mParent->mUID;
			if (parentID >= 0 && !validFlags[parentID]) validFlags[otherID] = false;
		}

		// detect co-centric parts
		for (int otherID = nodeID + 1; otherID < numAllNodes; otherID++) {
			if (!validFlags[otherID]) continue;
			if (TNode::detectCoCentricity(mAllNodes[nodeID], mAllNodes[otherID])) {
				mAllNodes[nodeID]->mCocentric.push_back(mAllNodes[otherID]);
				mAllNodes[otherID]->mCocentric.push_back(mAllNodes[nodeID]);
				/*
				if (true) {
					cout << "\nCo-centric: " << nodeID << " -- " << otherID << endl;
					auto &mesh1 = mAllNodes[nodeID]->mNodeDescriptors.mNodeMesh;
					auto &mesh2 = mAllNodes[otherID]->mNodeDescriptors.mNodeMesh;
					PlyExporter pe;
					if (!pe.addMesh(&mesh1.positions, &mesh1.normals, &mesh1.indices, vec3i(255, 0, 0))) return false;
					if (!pe.addMesh(&mesh2.positions, &mesh2.normals, &mesh2.indices, vec3i(0, 255, 255))) return false;
					if (!pe.output("cocentric.ply")) return false;
					system("pause");
				}
				*/
			}
		}

		// detect adjacent parts
		for (int otherID = nodeID + 1; otherID < numAllNodes; otherID++) {
			if (!validFlags[otherID]) continue;
			if (TNode::detectAdjacency(mAllNodes[nodeID], mAllNodes[otherID])) {
				mAllNodes[nodeID]->mAdjacent.push_back(mAllNodes[otherID]);
				mAllNodes[otherID]->mAdjacent.push_back(mAllNodes[nodeID]);
				int adjFlag = 0;
				if (TNode::detectContact(mAllNodes[nodeID], mAllNodes[otherID])) {
					mAllNodes[nodeID]->mContact.push_back(mAllNodes[otherID]);
					mAllNodes[otherID]->mContact.push_back(mAllNodes[nodeID]);
					adjFlag += 1;
				}
				if (TNode::detectSupport(mAllNodes[nodeID], mAllNodes[otherID])) {
					mAllNodes[nodeID]->mSupport.push_back(mAllNodes[otherID]);
					mAllNodes[otherID]->mSupport.push_back(mAllNodes[nodeID]);
					adjFlag += 2;
				}
				/*
				if (true) {
					cout << "\nAdjacent: " << nodeID << " -- " << otherID << ": " << adjFlag << endl;
					auto &mesh1 = mAllNodes[nodeID]->mNodeDescriptors.mNodeMesh;
					auto &mesh2 = mAllNodes[otherID]->mNodeDescriptors.mNodeMesh;
					PlyExporter pe;
					if (!pe.addMesh(&mesh1.positions, &mesh1.normals, &mesh1.indices, vec3i(255, 0, 0))) return false;
					if (!pe.addMesh(&mesh2.positions, &mesh2.normals, &mesh2.indices, vec3i(0, 255, 255))) return false;
					if (!pe.output("adjacent.ply")) return false;
					system("pause");
				}
				*/
			}
		}
	}
	cout << endl;

	// HACK: update orientation for sphere-like "leaf" node

	for (TNode *node : graphLevels[0]) { // only process "part" nodes
		if (node->mNodeDescriptors.mPrimitive == 2) { // sphere-like
			if (!node->mAdjacent.empty()) {
				Eigen::Vector3d nodeOrien = node->mNodeDescriptors.mMajorOrientation;
				int numSlots = 0;
				for (TNode *neighbor : node->mAdjacent) {
					// find slot position
					Eigen::Vector3d slotCenter;
					Eigen::AlignedBox3d slotBB;
					if (true) {
						SKDTree tree;
						SKDTreeData treeData;
						if (!SampleUtil::buildKdTree(node->mNodeDescriptors.mSampleMatP, tree, treeData)) return false;
						Eigen::VectorXi nnIndices;
						if (!SampleUtil::findNearestNeighbors(tree, neighbor->mNodeDescriptors.mSampleMatP, nnIndices)) return false;
						Eigen::Matrix3Xd matNNPoint;
						if (!SampleUtil::sliceMatrices(node->mNodeDescriptors.mSampleMatP, nnIndices, matNNPoint)) return false;
						Eigen::RowVectorXd vecNNDist = (matNNPoint - neighbor->mNodeDescriptors.mSampleMatP).colwise().norm();

						double maxSlotDist = min(node->mNodeDescriptors.mSizeEpsilon,
							neighbor->mNodeDescriptors.mSizeEpsilon) * 1.5; // UNDONE: param slot max distance threshold
						int numSlotPoints = 0;
						slotCenter.setZero();
						for (int id = 0; id < (int)vecNNDist.size(); id++) {
							if (vecNNDist[id] < maxSlotDist) {
								Eigen::Vector3d slotPoint = neighbor->mNodeDescriptors.mSampleMatP.col(id);
								slotCenter += slotPoint;
								slotBB.extend(slotPoint);
								numSlotPoints++;
							}
						}
						if (numSlotPoints) slotCenter *= 1.0 / numSlotPoints;
					}
					Eigen::Vector3d newDir = node->mNodeDescriptors.mSampleMatP.rowwise().mean() - slotCenter;
					if (true) {
						auto bbmin = slotBB.min();
						auto bbmax = slotBB.max();
						double radius = slotBB.diagonal().norm();
						int i = 0;
					}
					if (newDir.squaredNorm() && slotBB.diagonal().norm() < 0.1) { // UNDONE: param maximum slot BB radius
						nodeOrien = newDir.normalized();
						numSlots++;
						if (numSlots > 1) break;
					}
				}
				if (numSlots == 1) {
					node->mNodeDescriptors.mMajorOrientation = nodeOrien;
					cout << "Updated major orientation for node " << node->mUID << ":" << nodeOrien.transpose() << endl;
				}
			}
		}
	}

	return true;
}

bool ContextPartGraph::saveGraphHierarchy(string fileName) {

	ofstream file(fileName, ios::binary);
	if (!file.is_open()) {
		cout << "Error: cannot write to file " << fileName << endl;
		return false;
	}

	int numNodes = (int)mAllNodes.size();
	file.write((char*)&numNodes, sizeof(numNodes));
	for (TNode *node : mAllNodes) {
		if (!node->saveHierarchy(file)) return false;
	}

	file.close();

	return true;
}

bool ContextPartGraph::loadGraphHierarchy(
	string fileName,
	TNodeGen &nodeGenerator,
	TTriangleMesh *mesh,
	vector<vector<vector<int>>> *segments)
{

	ifstream file(fileName, ios::binary);
	if (!file.is_open()) {
		cout << "Error: cannot open file " << fileName << endl;
		return false;
	}

	mRootNode = nodeGenerator.generateNode();
	mRootNode->clearNode();
	mRootNode->mpGraphMesh = mesh;
	mRootNode->mpGraphSegments = segments;
	mRootNode->mUID = -1;

	int numNodes;
	file.read((char*)&numNodes, sizeof(numNodes));
	mAllNodes.resize(numNodes);
	for (int nodeID = 0; nodeID < numNodes; nodeID++) {
		mAllNodes[nodeID] = nodeGenerator.generateNode();
	}
	for (int nodeID = 0; nodeID < numNodes; nodeID++) {
		TNode *node = mAllNodes[nodeID];
		if (!node->loadHierarchy(file, &mAllNodes, mRootNode)) return false;
		node->mUID = nodeID;
		node->mpGraphMesh = mesh;
		node->mpGraphSegments = segments;
	}

	file.close();

	return true;
}

bool ContextPartGraph::saveGraphDescriptor(string fileName) {

	ofstream file(fileName, ios::binary);
	if (!file.is_open()) {
		cout << "Error: cannot write to file " << fileName << endl;
		return false;	
	}

	int numNodes = (int)mAllNodes.size();
	file.write((char*)&numNodes, sizeof(numNodes));
	for (TNode *node : mAllNodes) {
		if (!node->saveDescriptor(file)) return false;
	}

	file.close();

	return true;
}

bool ContextPartGraph::loadGraphDescriptor(string fileName) {

	ifstream file(fileName, ios::binary);
	if (!file.is_open()) {
		cout << "Error: cannot open file " << fileName << endl;
		return false;
	}

	int numNodes;
	file.read((char*)&numNodes, sizeof(numNodes));
	if (numNodes != mAllNodes.size()) return error("incorrect node number");
	for (int nodeID = 0; nodeID < numNodes; nodeID++) {
		TNode *node = mAllNodes[nodeID];
		if (!node->loadDescriptor(file)) return false;
		node->mUID = nodeID;
	}

	file.close();

	return true;
}

bool ContextPartGraph::saveGraphContext(string fileName) {

	ofstream file(fileName, ios::binary);
	if (!file.is_open()) {
		cout << "Error: cannot write to file " << fileName << endl;
		return false;
	}

	int numNodes = (int)mAllNodes.size();
	file.write((char*)&numNodes, sizeof(numNodes));
	for (TNode *node : mAllNodes) {
		if (!node->saveContext(file)) return false;
	}

	file.close();

	return true;
}

bool ContextPartGraph::loadGraphContext(string fileName) {

	ifstream file(fileName, ios::binary);
	if (!file.is_open()) {
		cout << "Error: cannot open file " << fileName << endl;
		return false;
	}

	int numNodes;
	file.read((char*)&numNodes, sizeof(numNodes));
	if (numNodes != mAllNodes.size()) return error("incorrect node number");
	for (int nodeID = 0; nodeID < numNodes; nodeID++) {
		TNode *node = mAllNodes[nodeID];
		if (!node->loadContext(file, &mAllNodes)) return false;
		node->mUID = nodeID;
	}

	file.close();

	return true;
}

bool ContextPartGraph::extractGraphLevels(vector<vector<TNode*>> &graphLevels) {

	graphLevels.clear();
	vector<TNode*> currentLevel(1, mRootNode);
	while (true) {
		vector<TNode*> nextLevel(0);
		for (TNode *node : currentLevel) {
			for (TNode *child : node->mChildren) {
				nextLevel.push_back(child);
			}
		}
		if (nextLevel.empty()) break;
		graphLevels.push_back(nextLevel);
		currentLevel.swap(nextLevel);
	}

	return true;
}