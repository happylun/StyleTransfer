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

#include "FeatureAsset.h"

#include <iostream>
#include <fstream>

using namespace StyleSynthesis;

FeatureAsset::FeatureAsset() {
}

FeatureAsset::~FeatureAsset() {
}

bool FeatureAsset::loadAllFeatures(string folderName) {

	string featureCurvatureName = folderName + "curvature.txt";
	string featureSDFName = folderName + "SDF.txt";
	string featureLFDName = folderName + "LFD.txt";
	string featureSDName = folderName + "SD.txt";
	string featureCurveName = folderName + "curve.txt";

	string featurePartLFDName = folderName + "part-LFD.txt";
	string featurePartSDName = folderName + "part-SD.txt";

	if (!FeatureSampleCurvature::loadFeature(featureCurvatureName, mCurvature)) return false;
	if (!loadFeature(featureSDFName, mSDF)) return false;
	if (!loadFeature(featureLFDName, mLFD)) return false;
	if (!loadFeature(featureSDName, mSD)) return false;
	if (!mCurve.loadFeature(featureCurveName)) return false;

	if (!loadPartFeature(featurePartLFDName, mPartLFD)) return false;
	if (!loadPartFeature(featurePartSDName, mPartSD)) return false;

	return true;
}

bool FeatureAsset::saveFeature(string fileName, vector<double> &feature) {

	ofstream outFile(fileName, ios::binary);
	if (!outFile.is_open()) {
		cout << "Error: cannot write to file " << fileName << endl;
		return false;
	}

	int dimension = (int)feature.size();
	outFile.write((char*)(&dimension), sizeof(dimension));

	for (auto &value : feature) {
		outFile.write((char*)(&value), sizeof(double));
	}
	outFile.close();

	return true;
}


bool FeatureAsset::loadFeature(string fileName, vector<double> &feature) {

	ifstream inFile(fileName, ios::binary);
	if (!inFile.is_open()) {
		cout << "Error: cannot open file " << fileName << endl;
		return false;
	}

	int dimension;
	inFile.read((char*)(&dimension), sizeof(dimension));
	feature.resize(dimension);

	for (auto &value : feature) {
		inFile.read((char*)(&value), sizeof(double));
	}
	inFile.close();

	return true;
}


bool FeatureAsset::savePartFeature(string fileName, vector<vector<double>> &feature) {

	ofstream outFile(fileName, ios::binary);
	if (!outFile.is_open()) {
		cout << "Error: cannot write to file " << fileName << endl;
		return false;
	}

	int numParts = (int)feature.size();
	int dimension = numParts ? (int)feature[0].size() : 0;
	outFile.write((char*)(&numParts), sizeof(numParts));
	outFile.write((char*)(&dimension), sizeof(dimension));

	for (auto &partFeature : feature) {
		if ((int)partFeature.size() != dimension) {
			cout << "Error: incompatible feature dimension" << endl;
			return false;
		}
		for (auto &value : partFeature) {
			outFile.write((char*)(&value), sizeof(double));
		}
	}
	outFile.close();

	return true;
}


bool FeatureAsset::loadPartFeature(string fileName, vector<vector<double>> &feature) {

	ifstream inFile(fileName, ios::binary);
	if (!inFile.is_open()) {
		cout << "Error: cannot open file " << fileName << endl;
		return false;
	}

	int numParts;
	int dimension;
	inFile.read((char*)(&numParts), sizeof(numParts));
	inFile.read((char*)(&dimension), sizeof(dimension));

	feature.resize(numParts);
	for (auto &partFeature : feature) {
		partFeature.resize(dimension);
		for (auto &value : partFeature) {
			inFile.read((char*)(&value), sizeof(double));
		}
	}
	inFile.close();

	return true;
}
