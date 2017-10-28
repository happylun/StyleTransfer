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

#include "PipelineSimilarityIO.h"

#include <iostream>
#include <fstream>

#include "Sample/SampleUtil.h"

#include "Feature/FeatureSaliency.h"

#include "Data/DataUtil.h"

#include "Similarity/SimilarityData.h"
#include "Similarity/SimilarityDistance.h"
#include "Similarity/SimilarityMetric.h"

#include "Data/StyleSynthesisConfig.h"

using namespace std;
using namespace StyleSynthesis;

bool PipelineSimilarityIO::process() {

	//if (!learningPhase()) return false;
	if (!evaluatingPhase()) return false;
	//if (!rankingPhase()) return false;

	return true;
}

bool PipelineSimilarityIO::learningPhase() {

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string tripletFileName = datasetPrefix + "study/triplets.txt";
	vector<vector<string>> allTriplets;

	ifstream tripletFile(tripletFileName);
	while (!tripletFile.eof()) {
		string line;
		getline(tripletFile, line);
		if (line.empty()) break;
		vector<string> tokens;
		StringUtil::split(line, ' ', tokens);
		if (tokens.size() != 3) {
			cout << "Error: invalid triplet " << line << endl;
			return false;
		}
		for (string &name : tokens) name = StringUtil::trim(name);
		allTriplets.push_back(tokens);
	}
	tripletFile.close();

	int numTriplets = (int)allTriplets.size();

	if (!SimilarityMetric::initMetric(datasetPrefix + "weights/init/")) return false;

	for (int tripletID = 0; tripletID < numTriplets; tripletID++) {
		auto &triplet = allTriplets[tripletID];
		cout << "=========== Processing " << triplet[0] << "--" << triplet[1] << "--" << triplet[2] << " ===========" << endl;

		pair<string, string> pairAB(triplet[0], triplet[1]);
		pair<string, string> pairAC(triplet[0], triplet[2]);
		if (pairAB.first > pairAB.second) swap(pairAB.first, pairAB.second);
		if (pairAC.first > pairAC.second) swap(pairAC.first, pairAC.second);

		cout << "------ A--B ------" << endl;
		if (!runModelPairs(pairAB.first, pairAB.second)) return false;
		cout << "------ A--C ------" << endl;
		if (!runModelPairs(pairAC.first, pairAC.second)) return false;
		cout << "------ A--B--C ------" << endl;
		if (!runModelTriplets(triplet[0], triplet[1], triplet[2], tripletID)) return false;
	}

	return true;
}

bool PipelineSimilarityIO::evaluatingPhase() {

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string meshListFileName = datasetPrefix + "mesh/mesh-list-all.txt";
	vector<string> meshNameList;
	ifstream meshListFile(meshListFileName);
	while (!meshListFile.eof()) {
		string line;
		getline(meshListFile, line);
		if (line.empty()) break;
		meshNameList.push_back(StringUtil::trim(line));
	}
	meshListFile.close();

	int numMesh = (int)meshNameList.size();

	if (!SimilarityMetric::initMetric(datasetPrefix + "weights/")) return false;

	for (int meshID = 0; meshID < numMesh; meshID++) {

		string meshName = meshNameList[meshID];
		cout << "Evaluating model " << meshName << endl;
		if (!evaluateSingleModel(meshName)) return error("single model error");

		//system("pause");
	}

	return true;
}

bool PipelineSimilarityIO::rankingPhase() {

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string meshListFileName = datasetPrefix + "mesh/mesh-list-all.txt";
	vector<string> meshNameList;
	ifstream meshListFile(meshListFileName);
	while (!meshListFile.eof()) {
		string line;
		getline(meshListFile, line);
		if (line.empty()) break;
		meshNameList.push_back(StringUtil::trim(line));
	}
	meshListFile.close();

	int numMesh = (int)meshNameList.size();

	vector<string> sourceNames;
	vector<int> sourceIDs;
	if (true) {
		StringUtil::split(StyleSynthesisConfig::mData_CustomString1, ' ', sourceNames); // UNDONE: custom string
		sourceIDs.resize(sourceNames.size());
		for (int sourceID = 0; sourceID < numMesh; sourceID++) {
			for (int id = 0; id < (int)sourceNames.size(); id++) {
				if (sourceNames[id] == meshNameList[sourceID]) sourceIDs[id] = sourceID;
			}
		}
	}

	Eigen::MatrixXd distanceMatrix(numMesh, numMesh);
	distanceMatrix.setZero();

	//for (int sourceID = 0; sourceID < numMesh; sourceID++) {
		//string sourceName = meshNameList[sourceID];
	for (int id = 0; id < (int)sourceNames.size(); id++) {
		int sourceID = sourceIDs[id];
		string sourceName = sourceNames[id];

		//if (sourceName != StyleSynthesisConfig::mData_CustomString2) continue;

		for (int targetID = 0; targetID < numMesh; targetID++) {
			if (sourceID == targetID) continue;
			string targetName = meshNameList[targetID];

			//if (targetName != StyleSynthesisConfig::mData_CustomString3) continue;

			pair<string, string> pairName(sourceName, targetName);
			if (sourceName > targetName) swap(pairName.first, pairName.second);

			cout << "========= Processing " << sourceName << "--" << targetName << " =========" << endl;

			if (!SimilarityMetric::initMetric(datasetPrefix + "weights/init/")) return false;
			if (!runModelPairs(pairName.first, pairName.second)) return false;

			cout << "========= Evaluating " << sourceName << "--" << targetName << " =========" << endl;

			if (!SimilarityMetric::initMetric(datasetPrefix + "weights/")) return false;
			if (!evaluateModelPairs(pairName.first, pairName.second, distanceMatrix(sourceID, targetID))) return false;

			//system("pause");
		}
	}

	string distanceMatrixName = datasetPrefix + "rank/distance-matrix.txt";
	if (!FileUtil::makedir(distanceMatrixName)) return false;
	ofstream distanceMatrixFile(distanceMatrixName);
	distanceMatrixFile << distanceMatrix << endl;
	distanceMatrixFile.close();

	return true;
}

bool PipelineSimilarityIO::runModelPairs(string sourceName, string targetName) {

	// names

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string sourceShortName = sourceName.substr(sourceName.find_last_of("/\\") + 1);
	string targetShortName = targetName.substr(targetName.find_last_of("/\\") + 1);

	string elementFolder = datasetPrefix + "element/" + sourceShortName + "--" + targetShortName + "/";
	if (!FileUtil::makedir(elementFolder)) return false;

	string elementVisualName = elementFolder + "vis-element.ply";

	string anyFileName = elementFolder + "data-distance.txt";
	if (FileUtil::existsfile(anyFileName)) return true; // early quit
	
	// data

	SimilarityData data;
	if (!data.loadData(sourceName, targetName)) return false;

	// algorithm

	SimilarityDistance sd(&data);
	if (!sd.process()) return false;
	if (!sd.output(elementFolder)) return false;
	if (!sd.visualize(elementVisualName)) return false;

	return true;
}

bool PipelineSimilarityIO::runModelTriplets(string meshNameA, string meshNameB, string meshNameC, int tripletID) {

	// names

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string meshAShortName = meshNameA.substr(meshNameA.find_last_of("/\\") + 1);
	string meshBShortName = meshNameB.substr(meshNameB.find_last_of("/\\") + 1);
	string meshCShortName = meshNameC.substr(meshNameC.find_last_of("/\\") + 1);

	bool swapAB = !(meshAShortName < meshBShortName);
	bool swapAC = !(meshAShortName < meshCShortName);
	string pairABName = swapAB ? (meshBShortName + "--" + meshAShortName) : (meshAShortName + "--" + meshBShortName);
	string pairACName = swapAC ? (meshCShortName + "--" + meshAShortName) : (meshAShortName + "--" + meshCShortName);

	string elementABFolder = datasetPrefix + "element/" + pairABName + "/";
	string elementACFolder = datasetPrefix + "element/" + pairACName + "/";

	string tripletFolder = datasetPrefix + "triplet/" + to_string(tripletID + 1) + "/"; // 1-indexed
	if (!FileUtil::makedir(tripletFolder)) return false;

	FileUtil::copyfile(elementABFolder + "data-distance.txt", tripletFolder + "dAB.txt");
	FileUtil::copyfile(elementACFolder + "data-distance.txt", tripletFolder + "dAC.txt");

	if (!swapAB) {
		FileUtil::copyfile(elementABFolder + "data-source-element.txt", tripletFolder + "msAB.txt");
		FileUtil::copyfile(elementABFolder + "data-target-element.txt", tripletFolder + "mtAB.txt");
		FileUtil::copyfile(elementABFolder + "data-source-unmatch.txt", tripletFolder + "usAB.txt");
		FileUtil::copyfile(elementABFolder + "data-target-unmatch.txt", tripletFolder + "utAB.txt");
	} else {
		FileUtil::copyfile(elementABFolder + "data-target-element.txt", tripletFolder + "msAB.txt");
		FileUtil::copyfile(elementABFolder + "data-source-element.txt", tripletFolder + "mtAB.txt");
		FileUtil::copyfile(elementABFolder + "data-target-unmatch.txt", tripletFolder + "usAB.txt");
		FileUtil::copyfile(elementABFolder + "data-source-unmatch.txt", tripletFolder + "utAB.txt");
	}

	if (!swapAC) {
		FileUtil::copyfile(elementACFolder + "data-source-element.txt", tripletFolder + "msAC.txt");
		FileUtil::copyfile(elementACFolder + "data-target-element.txt", tripletFolder + "mtAC.txt");
		FileUtil::copyfile(elementACFolder + "data-source-unmatch.txt", tripletFolder + "usAC.txt");
		FileUtil::copyfile(elementACFolder + "data-target-unmatch.txt", tripletFolder + "utAC.txt");
	} else {
		FileUtil::copyfile(elementACFolder + "data-target-element.txt", tripletFolder + "msAC.txt");
		FileUtil::copyfile(elementACFolder + "data-source-element.txt", tripletFolder + "mtAC.txt");
		FileUtil::copyfile(elementACFolder + "data-target-unmatch.txt", tripletFolder + "usAC.txt");
		FileUtil::copyfile(elementACFolder + "data-source-unmatch.txt", tripletFolder + "utAC.txt");
	}

	return true;
}

bool PipelineSimilarityIO::evaluateSingleModel(string meshName) {

	// names

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string sampleName = datasetPrefix + "sample/" + meshName + ".ply";

	string featureFolder = datasetPrefix + "feature/" + meshName + "/";

	string saliencyName = featureFolder + "saliency.txt";
	string saliencyVName = featureFolder + "saliency.ply";

	// data

	TSampleSet samples;
	Eigen::MatrixXd saliency;

	if (!SampleUtil::loadSample(sampleName, samples)) return false;
	if (!DataUtil::loadMatrixBinary(saliencyName, saliency)) return false;

	// algorithm

	if (!FeatureSaliency::visualize(saliencyVName, samples, saliency)) return false;

	return true;
}

bool PipelineSimilarityIO::evaluateModelPairs(string sourceName, string targetName, double &distance) {

	// names

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string sourceShortName = sourceName.substr(sourceName.find_last_of("/\\") + 1);
	string targetShortName = targetName.substr(targetName.find_last_of("/\\") + 1);

	string elementFolder = datasetPrefix + "element/" + sourceShortName + "--" + targetShortName + "/";
	string elementEvaluationName = elementFolder + "data-evaluation.txt";
	string elementVisualName = elementFolder + "vis-element.ply";

	// data

	SimilarityData data;
	if (!data.loadData(sourceName, targetName)) return false;

	SimilarityDistanceData distanceData;
	if (!distanceData.loadData(elementFolder)) return false;

	// algorithm

	if (!SimilarityMetric::evaluateSimilarity(data, distanceData, distance)) return false;
	if (!SimilarityDistance::visualizeMatchingElements(elementVisualName, data, distanceData)) return false;

	ofstream elementEvaluationFile(elementEvaluationName);
	elementEvaluationFile << distance << endl;
	elementEvaluationFile.close();

	return true;
}

bool PipelineSimilarityIO::error(string s) {
	cout << "Error: " << s << endl;
	return false;
}