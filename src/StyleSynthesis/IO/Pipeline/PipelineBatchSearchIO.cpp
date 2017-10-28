#include "PipelineBatchSearchIO.h"

#include <iostream>
#include <fstream>
#include <set>

#include "Mesh/MeshUtil.h"

#include "Segment/SegmentGroupApxCvx.h"
#include "Segment/SegmentUtil.h"

#include "Similarity/SimilarityMetric.h"
#include "Similarity/SimilarityDistance.h"

#include "Context/ContextPartGraph.h"
#include "Context/ContextPartGraphNodeGenerator.h"
#include "Context/ContextPartGraphMatch.h"
#include "Context/ContextPartGraphAssemble.h"
#include "Context/ContextPartGraphTabuSearchPart.h"

#include "PipelineMatchPartIO.h"

#include "Data/StyleSynthesisConfig.h"
#include "Data/DataUtil.h"

using namespace std;
using namespace StyleSynthesis;

bool PipelineBatchSearchIO::process() {

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

	if (!StyleSynthesisConfig::saveConfig("base.cfg")) return false;

	int numPairs = (int)pairList.size();
	for (int pairID = 0; pairID < numPairs; pairID++) {
		cout << "=========== Processing pair " << (pairID+1) << " / " << numPairs << " : " << pairList[pairID] << " ===========" << endl;

		if (!StyleSynthesisConfig::loadConfig("base.cfg")) return false;
		if (!runModelPairs(sourceList[pairID], targetList[pairID], pairList[pairID])) return false;
		if (!organizeResults(sourceList[pairID], targetList[pairID], pairList[pairID])) return false;
	}

	return true;
}

bool PipelineBatchSearchIO::runModelPairs(string sourceName, string targetName, string pairName) {

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

	string assembleFolder = datasetPrefix + "assemble/" + pairName + "/";
	if (!FileUtil::makedir(assembleFolder)) return false;

	string paramDataName = assembleFolder + "params.cfg";
	string matchDataName = assembleFolder + "match.txt";
	string sourceNodeContributionName = assembleFolder + "source-contribution.txt";
	string targetNodeContributionName = assembleFolder + "target-contribution.txt";
	string sourceNodeSaliencyName = assembleFolder + "source-saliency.txt";
	string targetNodeSaliencyName = assembleFolder + "target-saliency.txt";
	string sourceKeyGroupName = assembleFolder + "source-key-group.txt";
	string targetKeyGroupName = assembleFolder + "target-key-group.txt";

	string searchResultName = assembleFolder + "search-result.txt";

	if (FileUtil::existsfile(searchResultName)) return true; // early quit

	// override case-specific param data

	cout << "Loading parameters..." << endl;

	if (!FileUtil::existsfile(paramDataName)) {
		string goodName = datasetPrefix + "assemble-good/" + pairName + "/params.cfg";
		if (!FileUtil::copyfile(goodName, paramDataName)) return false;
	}
	if (FileUtil::existsfile(paramDataName)) {

		bool runLaga = StyleSynthesisConfig::mContext_UseLagaDescriptors;

		if (!StyleSynthesisConfig::loadConfig(paramDataName)) return false;

		if (!runLaga) {
			//StyleSynthesisConfig::mAssemble_DebugNoBranching = true;
		} else {
			StyleSynthesisConfig::mAssemble_DebugNoBranching = true;
			StyleSynthesisConfig::mContext_UseLagaDescriptors = true;
			StyleSynthesisConfig::mAssemble_AllowAddingParts = false;
			StyleSynthesisConfig::mAssemble_AllowRemovingParts = false;
			StyleSynthesisConfig::mAssemble_AllowPostAlignment = false;

			//StyleSynthesisConfig::mAssemble_PrimitiveAngleThreshold = 80.0;
			//StyleSynthesisConfig::mAssemble_InitialAlignmentMode = -1;
			//StyleSynthesisConfig::mAssemble_GlobalAlignmentMode = -1;
			//StyleSynthesisConfig::mAssemble_BestGuessAlignmentMode.values = { 2, 4, 2 };
		}
	} else return error("missing param data file");
	
	// graph data

	cout << "Loading graph data..." << endl;

	TTriangleMesh sourceMesh;
	TTriangleMesh targetMesh;
	vector<vector<vector<int>>> sourceSegments;
	vector<vector<vector<int>>> targetSegments;

	if (!MeshUtil::loadMesh(sourceMeshName, sourceMesh)) return false;
	if (!MeshUtil::loadMesh(targetMeshName, targetMesh)) return false;

	if (!SegmentGroupApxCvx::loadSegments(sourceSegmentName, sourceSegments)) return false;
	if (!SegmentGroupApxCvx::loadSegments(targetSegmentName, targetSegments)) return false;

	ContextPartGraphNodeGenerator cpgng;

	ContextPartGraph sourceGraph;
	if (!sourceGraph.loadGraphHierarchy(sourceGraphHName, cpgng, &sourceMesh, &sourceSegments)) return false;
	if (!sourceGraph.loadGraphDescriptor(sourceGraphDName)) return false;
	if (!sourceGraph.loadGraphContext(sourceGraphCName)) return false;

	ContextPartGraph targetGraph;
	if (!targetGraph.loadGraphHierarchy(targetGraphHName, cpgng, &targetMesh, &targetSegments)) return false;
	if (!targetGraph.loadGraphDescriptor(targetGraphDName)) return false;
	if (!targetGraph.loadGraphContext(targetGraphCName)) return false;

	// matching data

	cout << "Loading matching data..." << endl;

	double shapeStyleDistance;
	vector<double> sourceNodeContribution;
	vector<double> targetNodeContribution;
	vector<double> sourceNodeSaliency;
	vector<double> targetNodeSaliency;
	if (FileUtil::existsfile(sourceNodeContributionName) &&
		FileUtil::existsfile(targetNodeContributionName) &&
		FileUtil::existsfile(sourceNodeSaliencyName) &&
		FileUtil::existsfile(targetNodeSaliencyName))
	{
		if (!DataUtil::loadValueListASCII(sourceNodeContributionName, sourceNodeContribution)) return false;
		if (!DataUtil::loadValueListASCII(targetNodeContributionName, targetNodeContribution)) return false;
		if (!DataUtil::loadValueListASCII(sourceNodeSaliencyName, sourceNodeSaliency)) return false;
		if (!DataUtil::loadValueListASCII(targetNodeSaliencyName, targetNodeSaliency)) return false;
	} else {

		if (!PipelineMatchPartIO::computeNodeContribution(
			sourceName, targetName,
			&sourceGraph, &targetGraph,
			&sourceSegments, &targetSegments,
			shapeStyleDistance,
			sourceNodeContribution, targetNodeContribution,
			sourceNodeSaliency, targetNodeSaliency)) return false;
		if (!DataUtil::saveValueListASCII(sourceNodeContributionName, sourceNodeContribution)) return false;
		if (!DataUtil::saveValueListASCII(targetNodeContributionName, targetNodeContribution)) return false;
		if (!DataUtil::saveValueListASCII(sourceNodeSaliencyName, sourceNodeSaliency)) return false;
		if (!DataUtil::saveValueListASCII(targetNodeSaliencyName, targetNodeSaliency)) return false;
	}

	vector<vector<vec2i>> matchings;
	if (FileUtil::existsfile(matchDataName)) {
		if (!ContextPartGraphMatch::loadMatchings(matchDataName, matchings)) return false;
	} else return error("missing matching data file");

	vector<vector<int>> sourceKeyGroups, targetKeyGroups;
	if (FileUtil::existsfile(sourceKeyGroupName) && FileUtil::existsfile(targetKeyGroupName)) {
		if (!DataUtil::loadGroupListASCII(sourceKeyGroupName, sourceKeyGroups)) return false;
		if (!DataUtil::loadGroupListASCII(targetKeyGroupName, targetKeyGroups)) return false;
	} else return error("missing key group files");

	// main search algorithm

	cout << "Running tabu search..." << endl;

	if (true) {

		ContextPartGraphTabuSearchPart cpgtsp;
		if (!cpgtsp.loadGraph(&cpgng, &sourceGraph, &targetGraph)) return false;
		if (!cpgtsp.loadMatchings(matchings)) return false;
		if (!cpgtsp.loadContributions(
			sourceNodeContribution, targetNodeContribution,
			sourceNodeSaliency, targetNodeSaliency)) return false;
		if (!cpgtsp.loadKeyGroups(sourceKeyGroups, targetKeyGroups)) return false;
		if (!cpgtsp.loadNames(weightsFolder, assembleFolder)) return false;
		if (!cpgtsp.process()) return false;

		int numSolutions = cpgtsp.getNumSolutions();
		if (!DataUtil::saveIndexListASCII(searchResultName, vector<int>(1, numSolutions))) return false;
	}

	return true;
}

bool PipelineBatchSearchIO::organizeResults(string sourceName, string targetName, string pairName) {

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;
	string assembleFolder = datasetPrefix + "assemble/" + pairName + "/";

	string searchResultName = assembleFolder + "search-result.txt";
	if (!FileUtil::existsfile(searchResultName)) return error("run batch search first");

	string allOutFolder = datasetPrefix + "assemble-ply/" + pairName + "/";
	if (!FileUtil::makedir(allOutFolder)) return false;

	string outFolder = assembleFolder + "solution-search/";
	if (!FileUtil::makedir(outFolder)) return false;

	if (true) {
		string inExemplarName = datasetPrefix + "segment/" + sourceName + ".ply";
		string inCandidateName = datasetPrefix + "segment/" + targetName + ".ply";

		string allOutExemplarName = allOutFolder + "exemplar.ply";
		string allOutCandidateName = allOutFolder + "candidate.ply";
		if (FileUtil::existsfile(allOutExemplarName)) return true; // early quit
		if (!FileUtil::copyfile(inExemplarName, allOutExemplarName)) return false;
		if (!FileUtil::copyfile(inCandidateName, allOutCandidateName)) return false;

		string outExemplarName = outFolder + "exemplar.ply";
		string outCandidateName = outFolder + "candidate.ply";
		if (!FileUtil::copyfile(inExemplarName, outExemplarName)) return false;
		if (!FileUtil::copyfile(inCandidateName, outCandidateName)) return false;
	}

	if (true) {
		string solutionListName = assembleFolder + "solution-sorted-list.txt";
		ifstream solutionListFile(solutionListName);
		if (!solutionListFile.is_open()) return error("cannot open solution list file");
		int numSolutions = 0;
		while (!solutionListFile.eof()) {
			int solutionID;
			double solutionDistance;
			solutionListFile >> solutionID >> solutionDistance;
			if (solutionID == 0 || solutionListFile.fail()) break;

			string inResultName = assembleFolder + "solution-" + to_string(solutionID) + "/mesh.ply";

			string allOutResultName = allOutFolder + "mesh-" + to_string(numSolutions) + ".ply";
			if (!FileUtil::copyfile(inResultName, allOutResultName)) return false;

			string outResultName = outFolder + "mesh-" + to_string(numSolutions) + ".ply";
			if (!FileUtil::copyfile(inResultName, outResultName)) return false;

			numSolutions++;
		}
		solutionListFile.close();
	}

	return true;
}

bool PipelineBatchSearchIO::error(string s) {
	cout << "Error: " << s << endl;
	return false;
}