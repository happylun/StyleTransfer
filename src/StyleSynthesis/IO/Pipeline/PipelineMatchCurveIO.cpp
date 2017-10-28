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

#include "PipelineMatchCurveIO.h"

#include <iostream>
#include <fstream>

#include "Mesh/MeshUtil.h"

#include "Segment/SegmentGroupApxCvx.h"

#include "Context/ContextPartGraph.h"
#include "Context/ContextPartGraphNodeGenerator.h"
#include "Context/ContextPartGraphAssembleResult.h"
#include "Context/ContextPartGraphMatch.h"
#include "Context/ContextPartGraphMatchCurve.h"
#include "Context/ContextPartGraphTabuSearchCurve.h"

#include "Curve/CurveUtil.h"

#include "Data/StyleSynthesisConfig.h"
#include "Data/DataUtil.h"

using namespace std;
using namespace StyleSynthesis;

bool PipelineMatchCurveIO::process() {

	if (false) {
		string sourceName = StyleSynthesisConfig::mData_CustomString1;
		string targetName = StyleSynthesisConfig::mData_CustomString2;
		string sourceShortName = sourceName.substr(sourceName.find_last_of("/\\") + 1);
		string targetShortName = targetName.substr(targetName.find_last_of("/\\") + 1);
		string pairName = sourceShortName + "--" + targetShortName;

		if (!runModelPairs(sourceName, targetName, pairName)) return false;
		return true;
	}

	if (true) {

		string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

		vector<string> sourceList(0), targetList(0), pairList(0);

		string pairListName = datasetPrefix + "assemble/pair-list.txt";
		ifstream pairListFile(pairListName);
		if (!pairListFile.is_open()) return error("cannot open pair list file");
		while (!pairListFile.eof()) {
			string line;
			getline(pairListFile, line);
			line = StringUtil::trim(line);
			if (pairListFile.fail() || line.length() == 0) break;
			stringstream ss(line);
			string sourceName, targetName, pairName;
			ss >> sourceName >> targetName >> pairName;
			if (ss.fail()) return error("incorrect pair list token '" + line + "'");
			sourceList.push_back(sourceName);
			targetList.push_back(targetName);
			pairList.push_back(pairName);
		}
		pairListFile.close();

		int numPairs = (int)pairList.size();
		for (int pairID = 0; pairID < numPairs; pairID++) {
			cout << "=========== Processing pair " << (pairID + 1) << " / " << numPairs << " : " << pairList[pairID] << " ===========" << endl;

			if (!runModelPairs(sourceList[pairID], targetList[pairID], pairList[pairID])) return false;
			if (!organizeResults(sourceList[pairID], targetList[pairID], pairList[pairID])) return false;
		}
	}

	return true;
}

bool PipelineMatchCurveIO::runModelPairs(string sourceName, string targetName, string pairName) {

	// names

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string assembleFolder = datasetPrefix + "assemble/" + pairName + "/";
	string graphSimilarityName = assembleFolder + "similarity.txt";

	string solutionListName = assembleFolder + "solution-list.txt";
	string sourceNodeContributionName = assembleFolder + "source-contribution.txt";
	string targetNodeContributionName = assembleFolder + "target-contribution.txt";
	string sourceNodeSaliencyName = assembleFolder + "source-saliency.txt";
	string targetNodeSaliencyName = assembleFolder + "target-saliency.txt";

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

	string curveViewPointName = datasetPrefix + "curve/viewpoint.txt";
	string sourceCurveFolder = datasetPrefix + "curve/" + sourceName + "/";
	string targetCurveFolder = datasetPrefix + "curve/" + targetName + "/";

	string weightsFolder = datasetPrefix + "weights/";

	string deformFolder = datasetPrefix + "deform/" + pairName + "/";
	if (!FileUtil::makedir(deformFolder)) return false;

	// data

	TTriangleMesh sourceMesh;
	TTriangleMesh targetMesh;
	vector<vector<vector<int>>> sourceSegments;
	vector<vector<vector<int>>> targetSegments;

	if (!MeshUtil::loadMesh(sourceMeshName, sourceMesh)) return false;
	if (!MeshUtil::loadMesh(targetMeshName, targetMesh)) return false;

	if (!SegmentGroupApxCvx::loadSegments(sourceSegmentName, sourceSegments)) return false;
	if (!SegmentGroupApxCvx::loadSegments(targetSegmentName, targetSegments)) return false;

	ContextPartGraphNodeGenerator graphNodeGenerator;

	ContextPartGraph sourceGraph;
	if (!sourceGraph.loadGraphHierarchy(sourceGraphHName, graphNodeGenerator, &sourceMesh, &sourceSegments)) return false;
	if (!sourceGraph.loadGraphDescriptor(sourceGraphDName)) return false;
	if (!sourceGraph.loadGraphContext(sourceGraphCName)) return false;

	ContextPartGraph targetGraph;
	if (!targetGraph.loadGraphHierarchy(targetGraphHName, graphNodeGenerator, &targetMesh, &targetSegments)) return false;
	if (!targetGraph.loadGraphDescriptor(targetGraphDName)) return false;
	if (!targetGraph.loadGraphContext(targetGraphCName)) return false;

	Eigen::MatrixXd graphSimilarityMatrix;
	if (!DataUtil::loadMatrixASCII(graphSimilarityName, graphSimilarityMatrix)) return false;

	map<int, double> solutionScore;
	if (true) {
		Eigen::MatrixXd solutionScoreMatrix;
		if (!DataUtil::loadMatrixASCII(solutionListName, solutionScoreMatrix)) return false;
		for (int row = 0; row < (int)solutionScoreMatrix.rows(); row++) {
			int solID = (int)(solutionScoreMatrix(row, 0));
			double solScore = solutionScoreMatrix(row, 1);
			solutionScore[solID] = solScore;
		}
	}

	vector<double> sourceNodeContribution;
	vector<double> targetNodeContribution;
	vector<double> sourceNodeSaliency;
	vector<double> targetNodeSaliency;
	if (!DataUtil::loadValueListASCII(sourceNodeContributionName, sourceNodeContribution)) return false;
	if (!DataUtil::loadValueListASCII(targetNodeContributionName, targetNodeContribution)) return false;
	if (!DataUtil::loadValueListASCII(sourceNodeSaliencyName, sourceNodeSaliency)) return false;
	if (!DataUtil::loadValueListASCII(targetNodeSaliencyName, targetNodeSaliency)) return false;

	vector<vec3> curveViewPoints;
	if (!CurveUtil::loadViewPoints(curveViewPointName, curveViewPoints)) return false;

	// process each assembly solution

	vector<string> allResultList(0);
	vector<double> allScoreList(0);

	int assembleID = 0;
	while (true) {

		// names

		string assembleSolutionFolder = assembleFolder + "solution-" + to_string(assembleID) + "/";

		string assembleSolutionResultName = assembleSolutionFolder + "mesh.ply";
		if (!FileUtil::existsfile(assembleSolutionResultName)) break;

		string deformSolutionFolder = deformFolder + "solution-" + to_string(assembleID) + "/";
		if (!FileUtil::makedir(deformSolutionFolder)) return false;

		string deformSourceCurvesName = deformSolutionFolder + "sourceCurves.txt";
		string deformTargetCurvesName = deformSolutionFolder + "targetCurves.txt";
		string deformSourceNodesName = deformSolutionFolder + "sourceNodes.txt";
		string deformTargetNodesName = deformSolutionFolder + "targetNodes.txt";
		string deformCurveMatchingName = deformSolutionFolder + "matchCurves.txt";
		string deformCurveVisualName = deformSolutionFolder + "matchCurves.ply";
		string deformResultName = deformSolutionFolder + "result.txt";

		if (FileUtil::existsfile(deformResultName)) return true; // early quit

		// data

		ContextPartGraphAssembleResult assembleSolution;
		if (!assembleSolution.loadData(assembleSolutionFolder, graphNodeGenerator)) return false;

		// algorithm

		vector<vector<vector<vec3>>> deformSourceCurves;
		vector<vector<vector<vec3>>> deformTargetCurves;
		vector<vector<int>> deformSourceNodes;
		vector<vector<int>> deformTargetNodes;
		vector<int> deformCurveMatching; // I don't really need this. Just use it for manual override

		cout << "Matching curves..." << endl;

		if (FileUtil::existsfile(deformCurveMatchingName)) {

			if (!CurveUtil::loadCurveGroups(deformSourceCurvesName, deformSourceCurves)) return false;
			if (!CurveUtil::loadCurveGroups(deformTargetCurvesName, deformTargetCurves)) return false;
			if (!DataUtil::loadCellArraysBinary(deformSourceNodesName, deformSourceNodes)) return false;
			if (!DataUtil::loadCellArraysBinary(deformTargetNodesName, deformTargetNodes)) return false;

		} else {

			ContextPartGraphMatchCurve cpgmc;
			if (!cpgmc.loadWeights(weightsFolder)) return false;
			if (!cpgmc.loadGraph(&sourceGraph, &targetGraph, &assembleSolution, graphSimilarityMatrix)) return false;
			if (!cpgmc.loadCurve(sourceCurveFolder, targetCurveFolder, curveViewPoints)) return false;
			if (!cpgmc.process()) return false;
			if (!cpgmc.visualize(deformCurveVisualName)) return false;
			if (!cpgmc.output(deformSourceCurves, deformTargetCurves, deformSourceNodes, deformTargetNodes)) return false;

			int numGroups = (int)deformSourceCurves.size();
			deformCurveMatching.resize(numGroups);
			for (int k = 0; k < numGroups; k++) deformCurveMatching[k] = k;

			if (!CurveUtil::saveCurveGroups(deformSourceCurvesName, deformSourceCurves)) return false;
			if (!CurveUtil::saveCurveGroups(deformTargetCurvesName, deformTargetCurves)) return false;
			if (!DataUtil::saveCellArraysBinary(deformSourceNodesName, deformSourceNodes)) return false;
			if (!DataUtil::saveCellArraysBinary(deformTargetNodesName, deformTargetNodes)) return false;
			if (!DataUtil::saveIndexListASCII(deformCurveMatchingName, deformCurveMatching)) return false;
		}

		//system("pause");
		if (!DataUtil::loadIndexListASCII(deformCurveMatchingName, deformCurveMatching)) return false;
		if (true) {

			// HACK: only use the best curve group
			deformCurveMatching.resize(1);

			// prune groups
			vector<vector<vector<vec3>>> newSourceCurves(0);
			vector<vector<vector<vec3>>> newTargetCurves(0);
			vector<vector<int>> newSourceNodes(0);
			vector<vector<int>> newTargetNodes(0);
			for (int groupID : deformCurveMatching) {
				newSourceCurves.push_back(deformSourceCurves[groupID]);
				newTargetCurves.push_back(deformTargetCurves[groupID]);
				newSourceNodes.push_back(deformSourceNodes[groupID]);
				newTargetNodes.push_back(deformTargetNodes[groupID]);
			}
			deformSourceCurves.swap(newSourceCurves);
			deformTargetCurves.swap(newTargetCurves);
			deformSourceNodes.swap(newSourceNodes);
			deformTargetNodes.swap(newTargetNodes);
		}

		cout << "Running tabu search for curve deformation..." << endl;

		vector<double> deformScoreList;

		ContextPartGraphTabuSearchCurve cpgtsc;
		if (!cpgtsc.loadGraph(&sourceGraph, &targetGraph)) return false;
		if (!cpgtsc.loadAssemble(&assembleSolution)) return false;
		if (!cpgtsc.loadCurve(deformSourceCurves, deformTargetCurves, deformSourceNodes, deformTargetNodes)) return false;
		if (!cpgtsc.loadScore(sourceNodeContribution, targetNodeContribution, sourceNodeSaliency, targetNodeSaliency)) return false;
		if (!cpgtsc.loadNames(deformSolutionFolder)) return false;
		if (!cpgtsc.process()) return false;
		if (!cpgtsc.getSolutionScore(deformScoreList)) return false;

		if (assembleID > 0) { // skip null solution
			allResultList.push_back(assembleSolutionResultName);
			allScoreList.push_back(solutionScore[assembleID]);
		}
		for (int resID = 0; resID < (int)deformScoreList.size(); resID++) {
			string name = deformSolutionFolder + "solution-" + to_string(resID + 1) + "/result.ply";
			//string name = deformSolutionFolder + "solution-" + to_string(resID + 1) + "/colorMesh.ply";
			double score = deformScoreList[resID];
			allResultList.push_back(name);
			allScoreList.push_back(score);
		}

		if (!DataUtil::saveIndexListASCII(deformResultName, vector<int>(1, (int)deformScoreList.size()))) return false;

		assembleID++;
	}

	// organize all results

	cout << "Organizing all results..." << endl;

	if (true) {

		string allResultFolderName = deformFolder + "solution-all/";
		if (!FileUtil::makedir(allResultFolderName)) return false;

		int numResults = (int)allResultList.size();
		vector<int> orderList(numResults);
		for (int k = 0; k < numResults; k++) orderList[k] = k;
		sort(orderList.begin(), orderList.end(),
			[&allScoreList](int lhs, int rhs) {return allScoreList[lhs] < allScoreList[rhs]; });

		ofstream solutionSortedListFile(deformFolder + "solution-sorted-list.txt");
		for (int orderID = 0; orderID < numResults; orderID++) {
			int resultID = orderList[orderID];
			string fromName = allResultList[resultID];
			string toName = allResultFolderName + "mesh-" + to_string(orderID) + ".ply";
			if (!FileUtil::copyfile(fromName, toName)) return false;

			solutionSortedListFile << allScoreList[resultID] << "\t\t" << fromName << endl;
		}
		solutionSortedListFile.close();
	}

	return true;
}

bool PipelineMatchCurveIO::organizeResults(string sourceName, string targetName, string pairName) {

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;
	string assembleFolder = datasetPrefix + "assemble/" + pairName + "/solution-search/";
	string deformFolder = datasetPrefix + "deform/" + pairName + "/solution-all/";

	string outputFolder = datasetPrefix + "output/" + pairName + "/";
	if (!FileUtil::makedir(outputFolder)) return false;

	string exemplarName = assembleFolder + "exemplar.ply";
	string candidateName = assembleFolder + "candidate.ply";
	string outExemplarName = outputFolder + "exemplar.ply";
	string outCandidateName = outputFolder + "candidate.ply";
	if (!FileUtil::copyfile(exemplarName, outExemplarName)) return false;
	if (!FileUtil::copyfile(candidateName, outCandidateName)) return false;

	int outputID = 0;
	while (true) {
		string resultName = deformFolder + "mesh-" + to_string(outputID) + ".ply";
		if (!FileUtil::existsfile(resultName)) break; // no more results to copy out
		string outResultName = outputFolder + "mesh-" + to_string(outputID) + ".ply";
		//if (FileUtil::existsfile(outResultName)) break; // already copied out
		if (!FileUtil::copyfile(resultName, outResultName)) return false;
		outputID++;
	}

	return true;
}

bool PipelineMatchCurveIO::error(string s) {
	cout << "Error: " << s << endl;
	return false;
}