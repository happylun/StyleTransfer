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

#include "ContextPartGraphTabuSearchCurve.h"

#include <iostream>
#include <fstream>
#include <set>
#include <unordered_set>

#include "Mesh/MeshUtil.h"

#include "Context/ContextPartGraph.h"
#include "Context/ContextPartGraphAssembleResult.h"

#include "Curve/CurveDeformVertex.h"
#include "Curve/CurveUtil.h"

#include "Segment/SegmentUtil.h"

#include "Utility/FileUtil.h"
#include "Utility/PlyExporter.h"

#include "Data/DataUtil.h"
#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

#define OUTPUT_PROGRESS

ContextPartGraphTabuSearchCurve::ContextPartGraphTabuSearchCurve() {
}

ContextPartGraphTabuSearchCurve::~ContextPartGraphTabuSearchCurve() {
}

bool ContextPartGraphTabuSearchCurve::loadGraph(TGraph *source, TGraph *target) {

	mpSourceGraph = source;
	mpTargetGraph = target;

	return true;
}

bool ContextPartGraphTabuSearchCurve::loadAssemble(TSolution *solution) {

	mpAssembleSolution = solution;

	return true;
}

bool ContextPartGraphTabuSearchCurve::loadCurve(
	vector<vector<vector<vec3>>> &srcCurves,
	vector<vector<vector<vec3>>> &tgtCurves,
	vector<vector<int>> &srcNodes,
	vector<vector<int>> &tgtNodes)
{
	mMatchingSourceCurves = srcCurves;
	mMatchingTargetCurves = tgtCurves;
	mMatchingSourceNodes = srcNodes;
	mMatchingTargetNodes = tgtNodes;

	return true;
}

bool ContextPartGraphTabuSearchCurve::loadScore(
	vector<double> &sourceContrib,
	vector<double> &targetContrib,
	vector<double> &sourceSaliency,
	vector<double> &targetSaliency)
{
	mSourceNodeContribution = sourceContrib;
	mTargetNodeContribution = targetContrib;
	mSourceNodeSaliency = sourceSaliency;
	mTargetNodeSaliency = targetSaliency;

	return true;
}

bool ContextPartGraphTabuSearchCurve::loadNames(string solutionFolder) {

	mSolutionFolder = solutionFolder;

	return true;
}

bool ContextPartGraphTabuSearchCurve::process() {

	if (!initialize()) return false;
	if (!runTabuSearch()) return false;

	return true;
}

bool ContextPartGraphTabuSearchCurve::initialize() {

#ifdef OUTPUT_PROGRESS
	cout << "Initializing tabu search" << endl;
#endif

	// merge groups

	if (true) {

		int numGroups = (int)mMatchingSourceCurves.size();

		vector<bool> visitedFlags(numGroups, false);
		for (int groupID = 0; groupID < numGroups; groupID++) {
			if (visitedFlags[groupID]) continue;
			visitedFlags[groupID] = true;

			vector<vector<vec3>> &srcCurves = mMatchingSourceCurves[groupID];
			vector<vector<vec3>> &tgtCurves = mMatchingTargetCurves[groupID];
			vector<int> &srcNodes = mMatchingSourceNodes[groupID];
			vector<int> &tgtNodes = mMatchingTargetNodes[groupID];
			unordered_set<int> srcNodesSet(srcNodes.begin(), srcNodes.end());
			unordered_set<int> tgtNodesSet(tgtNodes.begin(), tgtNodes.end());

			for (int otherGroupID = groupID + 1; otherGroupID < numGroups; otherGroupID++) {
				if (visitedFlags[otherGroupID]) continue;

				// check whether all nodes are in the set

				vector<vector<vec3>> &otherSrcCurves = mMatchingSourceCurves[otherGroupID];
				vector<vector<vec3>> &otherTgtCurves = mMatchingTargetCurves[otherGroupID];
				vector<int> &otherSrcNodes = mMatchingSourceNodes[otherGroupID];
				vector<int> &otherTgtNodes = mMatchingTargetNodes[otherGroupID];
				bool canMerged = true;
				for (int nodeID : otherSrcNodes) {
					if (srcNodesSet.find(nodeID) == srcNodesSet.end()) {
						canMerged = false;
						break;
					}
				}
				if (!canMerged) continue;
				for (int nodeID : otherTgtNodes) {
					if (tgtNodesSet.find(nodeID) == tgtNodesSet.end()) {
						canMerged = false;
						break;
					}
				}
				if (!canMerged) continue;

				// merge curve group

				srcCurves.insert(srcCurves.end(), otherSrcCurves.begin(), otherSrcCurves.end());
				tgtCurves.insert(tgtCurves.end(), otherTgtCurves.begin(), otherTgtCurves.end());
				otherSrcCurves.clear();
				otherTgtCurves.clear();
				otherSrcNodes.clear();
				otherTgtNodes.clear();

				visitedFlags[otherGroupID] = true;
			}
		}

		// remove empty groups

		vector<vector<vector<vec3>>> newSrcCurves(0);
		vector<vector<vector<vec3>>> newTgtCurves(0);
		vector<vector<int>> newSrcNodes(0);
		vector<vector<int>> newTgtNodes(0);

		for (int groupID = 0; groupID < numGroups; groupID++) {
			if (!mMatchingSourceCurves[groupID].empty()) {
				newSrcCurves.push_back(mMatchingSourceCurves[groupID]);
				newTgtCurves.push_back(mMatchingTargetCurves[groupID]);
				newSrcNodes.push_back(mMatchingSourceNodes[groupID]);
				newTgtNodes.push_back(mMatchingTargetNodes[groupID]);
				cout << "Group:" << endl;
				cout << " source curves: " << mMatchingSourceCurves[groupID].size() << " curves" << endl;
				cout << " source nodes:";
				for (int id : mMatchingSourceNodes[groupID]) cout << " " << id;
				cout << endl;
				cout << " target curves: " << mMatchingTargetCurves[groupID].size() << " curves" << endl;
				cout << " target nodes:";
				for (int id : mMatchingTargetNodes[groupID]) cout << " " << id;
				cout << endl;
			}
		}

		mMatchingSourceCurves.swap(newSrcCurves);
		mMatchingTargetCurves.swap(newTgtCurves);
		mMatchingSourceNodes.swap(newSrcNodes);
		mMatchingTargetNodes.swap(newTgtNodes);
	}

	mSolutionScore.clear();

	return true;
}

bool ContextPartGraphTabuSearchCurve::runTabuSearch() {

#ifdef OUTPUT_PROGRESS
	cout << "Running tabu search" << endl;
#endif

	TTriangleMesh exemplarMesh = *mpSourceGraph->mRootNode->mpGraphMesh;
	TTriangleMesh assembleMesh = *mpAssembleSolution->mGraph.mRootNode->mpGraphMesh;

	int numMatchGroups = (int)mMatchingSourceCurves.size();
	long long maxCombination = 1ll << numMatchGroups;

	for (long long combination = 1; combination < maxCombination; combination++) {
	//for (long long combination = maxCombination-1ll; combination < maxCombination; combination++) {
	//for (long long combination = 1; combination < maxCombination; combination = combination * 2 + 1) {

#ifdef OUTPUT_PROGRESS
		cout << "=============== Running tabu search " << combination << " / " << maxCombination << " ===============" << endl;
#endif

		// names

		string subSolutionFolder = mSolutionFolder + "solution-" + to_string(combination) + "/";
		if (!FileUtil::makedir(subSolutionFolder)) return false;

		string outResultName = subSolutionFolder + "result.ply";
		string outColorMeshName = subSolutionFolder + "colorMesh.ply";
		string outVisualDeformName = subSolutionFolder + "visual-deform.ply";
		string outVisualMatchName = subSolutionFolder + "visual-match.ply";
		string outStyleScoreName = subSolutionFolder + "style-score.txt";

		// get group set

		vector<int> groupSet;
		if (true) {
			long long currentBits = combination;
			for (int groupID = 0; groupID < numMatchGroups; groupID++) {
				if (currentBits % 2) groupSet.push_back(groupID);
				currentBits >>= 1;
			}
		}

		// gather involved nodes

		set<int> srcNodesSet;
		set<int> tgtNodesSet;
		for (int groupID : groupSet) {
			for (int nodeID : mMatchingSourceNodes[groupID]) srcNodesSet.insert(nodeID);
			for (int nodeID : mMatchingTargetNodes[groupID]) tgtNodesSet.insert(nodeID);
		}

		PlyExporter collector;
		vector<int> deformableFaceNode(assembleMesh.indices.size(), -1);

		// mark faces on added node

		if (true) {
			int numNodes = (int)mpAssembleSolution->mReplaceMapping.size();
			for (int nodeID = 0; nodeID < numNodes; nodeID++) {
				vec2i rm = mpAssembleSolution->mReplaceMapping[nodeID];
				if (rm[0] >= 0 && rm[1] < 0) {
					TNode *node = mpAssembleSolution->mGraph.mAllNodes[nodeID];
					for (int faceID : (*node->mpGraphSegments)[node->mPartLevelID][node->mPartSegmentID]) {
						deformableFaceNode[faceID] = -2;
					}
				}
			}
		}

		// handle each target node

		int numDeformations = 0;
		for (int tgtNodeID : tgtNodesSet) {

			int mappedNodeID = mpAssembleSolution->mNodeMapping[tgtNodeID];
			if (mappedNodeID < 0) continue;
			Eigen::Affine3d nodeXform = mpAssembleSolution->mNodeTransformation[tgtNodeID];

			// gather curves
			vector<vector<vector<vec3>>> deformSourceCurvesGroup;
			vector<vector<vector<vec3>>> deformTargetCurvesGroup;
			for (int groupID : groupSet) {
				vector<vector<vec3>> srcCurves(0), tgtCurves(0);
				int numPairs = (int)mMatchingSourceCurves[groupID].size();
				for (int pairID = 0; pairID < numPairs; pairID++) {
					if (mMatchingTargetNodes[groupID][pairID] == tgtNodeID) {
						srcCurves.push_back(mMatchingSourceCurves[groupID][pairID]);
						tgtCurves.push_back(mMatchingTargetCurves[groupID][pairID]);
					}
				}
				deformSourceCurvesGroup.push_back(srcCurves);
				deformTargetCurvesGroup.push_back(tgtCurves);
			}

			// transform curves
			for (auto &group : deformTargetCurvesGroup) {
				for (auto &curve : group) {
					for (vec3 &point : curve) {
						Eigen::Vector3d v(vec3d(point).data());
						v = nodeXform * v;
						point = vec3d(v[0], v[1], v[2]);
					}
				}
			}

			// gather deformable mesh			
			TTriangleMesh deformableMesh;
			if (true) {
				TNode *node = mpAssembleSolution->mGraph.mAllNodes[mappedNodeID];
				vector<int> &segment = (*node->mpGraphSegments)[node->mPartLevelID][node->mPartSegmentID];
				for (int faceID : segment) deformableFaceNode[faceID] = mappedNodeID;
				if (!MeshUtil::extractSubMesh(assembleMesh, segment, deformableMesh)) return false;
			}

			// do deformation (separately for each node)
			if (true) {

				numDeformations++;
				string deformationFolder = subSolutionFolder + "deform-" + to_string(numDeformations) + "/";
				if (!FileUtil::makedir(deformationFolder)) return false;

				TTriangleMesh deformingMesh;
				if (!subdivideMesh(deformableMesh, deformingMesh)) return false;

				CurveDeformVertex cdv;
				if (!cdv.loadMesh(deformingMesh)) return false;
				if (!cdv.loadCurves(deformSourceCurvesGroup, deformTargetCurvesGroup)) return false;
				if (!cdv.loadNames(deformationFolder)) return false;
				if (!cdv.process()) return false;
				if (!cdv.visualize()) return false;
				if (!cdv.output(deformableMesh)) return false;
			}

			// add to result
			vec3i deformColor = SegmentUtil::colorMapping(tgtNodeID);
			//if (!collector.addMesh(&deformableMesh.positions, &deformableMesh.normals, &deformableMesh.indices, deformColor)) return false;
			if (!collector.addMesh(&deformableMesh.positions, &deformableMesh.normals, &deformableMesh.indices, vec3i(255, 255, 255))) return false;
		}

		// extract base mesh

		TTriangleMesh baseMesh;
		if (true) {
			vector<int> baseFaces;
			vector<vec3i> baseColors(0);
			for (int faceID = 0; faceID < (int)assembleMesh.indices.size(); faceID++) {
				if (deformableFaceNode[faceID] < 0) {
					baseFaces.push_back(faceID);
					if (deformableFaceNode[faceID] < -1) {
						baseColors.push_back(vec3i(128, 255, 128));
					} else {
						baseColors.push_back(vec3i(255, 255, 255));
					}
				}
			}
			if (!MeshUtil::extractSubMesh(assembleMesh, baseFaces, baseMesh)) return false;
			//if (!collector.addMesh(&baseMesh.positions, &baseMesh.normals, &baseMesh.indices, vec3i(127, 127, 127))) return false;
			if (!collector.addMesh(&baseMesh.positions, &baseMesh.normals, &baseMesh.indices, &baseColors)) return false;
		}

		// organize whole mesh

		TTriangleMesh wholeMesh;
		if (true) {
			wholeMesh.positions = collector.mVertices;
			wholeMesh.normals = collector.mNormals;
			wholeMesh.indices = collector.mFaceIndices;
			wholeMesh.amount = (int)wholeMesh.positions.size();
		}

		if (!MeshUtil::saveMesh(outResultName, wholeMesh)) return false;
		if (!collector.output(outColorMeshName)) return false;

		// compute style score

		if (true) {
			double score;
			if (!computeSolutionScore(srcNodesSet, tgtNodesSet, score)) return false;
			mSolutionScore.push_back(score);
			if (!DataUtil::saveValueListASCII(outStyleScoreName, vector<double>(1, score))) return false;
		}

		// visualize involved curves

		if (true) {

			vec3i baseColor(127, 127, 127);
			vec3i nodeColor(255, 255, 255);
			//vec3i sourceColor(255, 0, 0);
			//vec3i targetColor(0, 255, 0);
			vec3i sourceColor(0, 255, 0);
			vec3i targetColor(0, 255, 0);
			double tubeRadius = StyleSynthesisConfig::mCurve_VisualizationTubeRadius;

			vector<vec3i> exemplarFaceColors(exemplarMesh.indices.size(), baseColor);
			for (int srcNodeID : srcNodesSet) {
				TNode *node = mpSourceGraph->mAllNodes[srcNodeID];
				vector<int> &segment = (*node->mpGraphSegments)[node->mPartLevelID][node->mPartSegmentID];
				for (int faceID : segment) exemplarFaceColors[faceID] = nodeColor;
			}

			vector<vec3i> assembleFaceColors(assembleMesh.indices.size(), baseColor);
			for (int faceID = 0; faceID < (int)assembleMesh.indices.size(); faceID++) {
				if (deformableFaceNode[faceID] >= 0) assembleFaceColors[faceID] = nodeColor;
			}

			vec3 bbEMin, bbEMax;
			vec3 bbAMin, bbAMax;
			if (!MeshUtil::computeAABB(exemplarMesh, bbEMin, bbEMax)) return false;
			if (!MeshUtil::computeAABB(assembleMesh, bbAMin, bbAMax)) return false;
			vec3 offset((bbEMax[0] - bbEMin[0] + bbAMax[0] - bbAMin[0])*0.6f, 0.0f, 0.0f);

			PlyExporter pe;
			PlyExporter peE, peD;

			for (int groupID : groupSet) {
				int numPairs = (int)mMatchingSourceCurves[groupID].size();
				for (int pairID = 0; pairID < numPairs; pairID++) {
					auto &srcCurve = mMatchingSourceCurves[groupID][pairID];
					auto &tgtCurve = mMatchingTargetCurves[groupID][pairID];
					TTriangleMesh srcTube, tgtTube;
					if (!CurveUtil::makeTube(srcCurve, srcTube, (float)tubeRadius)) return false;
					if (!CurveUtil::makeTube(tgtCurve, tgtTube, (float)tubeRadius)) return false;

					if (!pe.addMesh(&srcTube.positions, &srcTube.normals, &srcTube.indices, sourceColor)) return false;
					if (!pe.addMesh(&tgtTube.positions, &tgtTube.normals, &tgtTube.indices, targetColor, offset)) return false;

					if (!peE.addMesh(&srcTube.positions, &srcTube.normals, &srcTube.indices, sourceColor)) return false;
					if (!peD.addMesh(&tgtTube.positions, &tgtTube.normals, &tgtTube.indices, targetColor)) return false;
				}
			}
			if (!pe.addMesh(&exemplarMesh.positions, &exemplarMesh.normals, &exemplarMesh.indices, &exemplarFaceColors)) return false;
			if (!pe.addMesh(&assembleMesh.positions, &assembleMesh.normals, &assembleMesh.indices, &assembleFaceColors, offset)) return false;

			if (!peE.addMesh(&exemplarMesh.positions, &exemplarMesh.normals, &exemplarMesh.indices, vec3i(255,255,255))) return false;
			if (!peD.addMesh(&assembleMesh.positions, &assembleMesh.normals, &assembleMesh.indices, vec3i(255,255,255))) return false;

			if (!pe.output(outVisualMatchName)) return false;

			if (!peE.output(subSolutionFolder + "curve-exemplar.ply")) return false;
			if (!peD.output(subSolutionFolder + "curve-seed.ply")) return false;
		}

		//system("pause");
	}

	return true;
}

bool ContextPartGraphTabuSearchCurve::subdivideMesh(TTriangleMesh &inMesh, TTriangleMesh &outMesh) {

	vector<int> faceIndices;
	vector<int> vertexIndices;
	double subdivRadius = StyleSynthesisConfig::mDeform_MeshSubdivisionRadius;

	if (!MeshUtil::removeDegeneratedFaces(inMesh, outMesh, faceIndices)) return false;
	if (!MeshUtil::removeDuplicateVertices(outMesh, outMesh, vertexIndices)) return false;
	while (true) {
		int numFaces = (int)outMesh.indices.size();
		if (!MeshUtil::subdivideMeshMidPoint(outMesh, outMesh, faceIndices, subdivRadius)) return false;
		if ((int)faceIndices.size() <= numFaces) break; // can not subdivide any more
	}

	return true;
}

bool ContextPartGraphTabuSearchCurve::computeSolutionScore(set<int> &srcNodes, set<int> &tgtNodes, double &outScore) {

	vector<int> srcContrib(0), tgtContrib(0);
	for (vec2i &nodePair : mpAssembleSolution->mReplaceMapping) {
		if (nodePair[0] >= 0) srcContrib.push_back(nodePair[0]);
		if (nodePair[1] >= 0) tgtContrib.push_back(nodePair[1]);
	}

	set<int> srcSet(srcContrib.begin(), srcContrib.end());
	set<int> tgtSet(tgtContrib.begin(), tgtContrib.end());

	for (auto &node : mpAssembleSolution->mGraph.mRootNode->mChildren) {
		int nodeID = node->mUID;
		if (mpAssembleSolution->mNodeMapping[nodeID] < 0 && tgtSet.find(nodeID) == tgtSet.end()) {
			// insert removed node
			tgtSet.insert(nodeID);
		}
	}

	outScore = 2.0;
	for (int nodeID : srcSet) outScore -= mSourceNodeContribution[nodeID];
	for (int nodeID : tgtSet) outScore -= mTargetNodeContribution[nodeID];
	for (int nodeID : tgtNodes) outScore -= mTargetNodeContribution[nodeID];
	outScore *= 0.5;

	// re-normalize saliency
	double oldSaliency = 0;
	for (TNode *node : mpTargetGraph->mRootNode->mChildren) oldSaliency += mTargetNodeSaliency[node->mUID];
	double newSaliency = oldSaliency;
	for (int nodeID : srcContrib) newSaliency += mSourceNodeSaliency[nodeID];
	for (int nodeID : tgtContrib) newSaliency -= mTargetNodeSaliency[nodeID];
	double srcArea = 0;
	double tgtArea = 0;
	for (int nodeID : srcNodes) {
		TNode *node = mpSourceGraph->mAllNodes[nodeID];
		TTriangleMesh nodeMesh;
		vector<int> &nodeSegment = (*node->mpGraphSegments)[node->mPartLevelID][node->mPartSegmentID];
		if (!MeshUtil::extractSubMesh(*node->mpGraphMesh, nodeSegment, nodeMesh)) return false;
		double nodeArea;
		if (!MeshUtil::computeFaceArea(nodeMesh, nodeArea)) return false;
		srcArea += nodeArea;
	}
	for (int nodeID : tgtNodes) {
		TNode *node = mpTargetGraph->mAllNodes[nodeID];
		TTriangleMesh nodeMesh;
		vector<int> &nodeSegment = (*node->mpGraphSegments)[node->mPartLevelID][node->mPartSegmentID];
		if (!MeshUtil::extractSubMesh(*node->mpGraphMesh, nodeSegment, nodeMesh)) return false;
		double nodeArea;
		if (!MeshUtil::computeFaceArea(nodeMesh, nodeArea)) return false;
		tgtArea += nodeArea;
	}
	double areaFactor = srcArea ? tgtArea / srcArea : 1.0;
	for (int nodeID : srcNodes) newSaliency += mSourceNodeSaliency[nodeID] * areaFactor;
	for (int nodeID : tgtNodes) newSaliency -= mTargetNodeSaliency[nodeID];

	cout << "New saliency = " << newSaliency << endl;
	if(newSaliency > 0) outScore *= oldSaliency / newSaliency;
	outScore += 1e-6; // prevent numerical unstability...

	cout << "Style score = " << outScore << endl;

	return true;
}