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

#include "DataUtil.h"

#include <fstream>

using namespace StyleSynthesis;

bool DataUtil::saveVectorASCII(string fileName, Eigen::VectorXd &inVector) {

	ofstream file(fileName);
	if (!file.is_open()) {
		cout << "Error: cannot write to file " << fileName << endl;
		return false;
	}
	file << inVector << endl;
	file.close();

	return true;
}

bool DataUtil::loadVectorASCII(string fileName, Eigen::VectorXd &outVector) {

	ifstream file(fileName);
	if (!file.is_open()) {
		cout << "Error: cannot open file " << fileName << endl;
		return false;
	}

	vector<double> data;
	while (!file.eof()) {
		double value;
		file >> value;
		if (file.fail()) break;
		data.push_back(value);
	}
	file.close();

	outVector.resize(data.size());
	for (int pos = 0; pos < (int)data.size(); pos++) {
		outVector(pos) = data[pos];
	}

	return true;
}

bool DataUtil::saveRowVectorASCII(string fileName, Eigen::RowVectorXd &inVector) {

	ofstream file(fileName);
	if (!file.is_open()) {
		cout << "Error: cannot write to file " << fileName << endl;
		return false;
	}
	file << inVector << endl;
	file.close();

	return true;
}

bool DataUtil::loadRowVectorASCII(string fileName, Eigen::RowVectorXd &outVector) {

	ifstream file(fileName);
	if (!file.is_open()) {
		cout << "Error: cannot open file " << fileName << endl;
		return false;
	}

	vector<double> data;
	while (!file.eof()) {
		double value;
		file >> value;
		if (file.fail()) break;
		data.push_back(value);
	}
	file.close();

	outVector.resize(data.size());
	for (int pos = 0; pos < (int)data.size(); pos++) {
		outVector(pos) = data[pos];
	}

	return true;
}

bool DataUtil::saveMatrixASCII(string fileName, Eigen::MatrixXd &inMatrix) {

	ofstream file(fileName);
	if (!file.is_open()) {
		cout << "Error: cannot write to file " << fileName << endl;
		return false;
	}
	file << inMatrix << endl;
	file.close();

	return true;
}

bool DataUtil::loadMatrixASCII(string fileName, Eigen::MatrixXd &outMatrix) {

	ifstream file(fileName);
	if (!file.is_open()) {
		cout << "Error: cannot open file " << fileName << endl;
		return false;
	}
	int dim = -1;
	vector<vector<double>> data;
	while (!file.eof()) {
		string line;
		getline(file, line);
		if (!file.good()) break;
		stringstream ss(line);
		vector<double> row;
		while (!ss.eof()) {
			double value;
			ss >> value;
			if (ss.fail()) break;
			row.push_back(value);
		}
		if (row.empty()) break;
		if (dim == -1) dim = (int)row.size();
		if (dim != (int)row.size()) {
			cout << "Error: input matrix has inconsistent dimensions" << endl;
			return false;
		}
		data.push_back(row);
	}
	file.close();

	if (dim < 0) {
		outMatrix.resize(0, 0);
		return true;
	}

	outMatrix.resize(data.size(), dim);
	for (int row = 0; row < (int)data.size(); row++) {
		for (int col = 0; col < dim; col++) {
			outMatrix(row, col) = data[row][col];
		}
	}

	return true;
}


bool DataUtil::savePairListASCII(string fileName, vector<vec2i> &inPairList) {

	ofstream file(fileName);
	if (!file.is_open()) {
		cout << "Error: cannot write to file " << fileName << endl;
		return false;
	}

	int numPairs = (int)inPairList.size();
	file << numPairs << endl;
	for (vec2i inPair : inPairList) {
		file << inPair << endl;
	}

	file.close();

	return true;
}

bool DataUtil::loadPairListASCII(string fileName, vector<vec2i> &outPairList) {

	ifstream file(fileName);
	if (!file.is_open()) {
		cout << "Error: cannot open file " << fileName << endl;
		return false;
	}

	int numPairs;
	file >> numPairs;
	outPairList.resize(numPairs);
	for (vec2i &outPair : outPairList) {
		file >> outPair[0] >> outPair[1];
	}

	file.close();

	return true;
}

bool DataUtil::saveIndexListASCII(string fileName, vector<int> &inIndexList) {

	ofstream file(fileName);
	if (!file.is_open()) {
		cout << "Error: cannot write to file " << fileName << endl;
		return false;
	}

	for (int index : inIndexList) file << index << " ";

	file.close();

	return true;
}

bool DataUtil::loadIndexListASCII(string fileName, vector<int> &outIndexList) {

	ifstream file(fileName);
	if (!file.is_open()) {
		cout << "Error: cannot open file " << fileName << endl;
		return false;
	}

	outIndexList.clear();
	while (!file.eof()) {
		int index;
		file >> index;
		if (file.fail()) break;
		outIndexList.push_back(index);
	}

	file.close();

	return true;
}

bool DataUtil::saveValueListASCII(string fileName, vector<double> &inValueList) {

	ofstream file(fileName);
	if (!file.is_open()) {
		cout << "Error: cannot write to file " << fileName << endl;
		return false;
	}

	for (double value : inValueList) file << value << " ";

	file.close();

	return true;
}

bool DataUtil::loadValueListASCII(string fileName, vector<double> &outValueList) {

	ifstream file(fileName);
	if (!file.is_open()) {
		cout << "Error: cannot open file " << fileName << endl;
		return false;
	}

	outValueList.clear();
	while (!file.eof()) {
		double value;
		file >> value;
		if (file.fail()) break;
		outValueList.push_back(value);
	}

	file.close();

	return true;
}

bool DataUtil::saveGroupListASCII(string fileName, vector<vector<int>> &inGroupList) {

	ofstream file(fileName);
	if (!file.is_open()) {
		cout << "Error: cannot write to file " << fileName << endl;
		return false;
	}

	file << inGroupList.size() << endl;
	for (auto &group : inGroupList) {
		file << group.size();
		for (int index : group) file << " " << index;
		file << endl;
	}

	file.close();

	return true;
}

bool DataUtil::loadGroupListASCII(string fileName, vector<vector<int>> &outGroupList) {

	ifstream file(fileName);
	if (!file.is_open()) {
		cout << "Error: cannot open file " << fileName << endl;
		return false;
	}

	int numGroups;
	file >> numGroups;
	outGroupList.resize(numGroups);
	for (auto &group : outGroupList) {
		int numEntries;
		file >> numEntries;
		group.resize(numEntries);
		for (int &index : group) file >> index;
	}

	file.close();

	return true;
}

bool DataUtil::saveVectorBinary(string fileName, Eigen::VectorXd &inVector) {

	ofstream file(fileName, ios::binary);
	if (!file.is_open()) {
		cout << "Error: cannot write to file " << fileName << endl;
		return false;
	}
	int length = (int)inVector.size();
	file.write((const char*)(&length), sizeof(length));
	for (int l = 0; l < length; l++) {
		double value = inVector(l);
		file.write((const char*)(&value), sizeof(value));
	}
	file.close();

	return true;
}

bool DataUtil::loadVectorBinary(string fileName, Eigen::VectorXd &outVector) {

	ifstream file(fileName, ios::binary);
	if (!file.is_open()) {
		cout << "Error: cannot open file " << fileName << endl;
		return false;
	}
	int length;
	file.read((char*)(&length), sizeof(length));
	outVector.resize(length);
	for (int l = 0; l < length; l++) {
		double value;
		file.read((char*)(&value), sizeof(value));
		outVector(l) = value;
	}
	file.close();

	return true;
}

bool DataUtil::saveRowVectorBinary(string fileName, Eigen::RowVectorXd &inVector) {

	ofstream file(fileName, ios::binary);
	if (!file.is_open()) {
		cout << "Error: cannot write to file " << fileName << endl;
		return false;
	}
	int length = (int)inVector.size();
	file.write((const char*)(&length), sizeof(length));
	for (int l = 0; l < length; l++) {
		double value = inVector(l);
		file.write((const char*)(&value), sizeof(value));
	}
	file.close();

	return true;
}

bool DataUtil::loadRowVectorBinary(string fileName, Eigen::RowVectorXd &outVector) {

	ifstream file(fileName, ios::binary);
	if (!file.is_open()) {
		cout << "Error: cannot open file " << fileName << endl;
		return false;
	}
	int length;
	file.read((char*)(&length), sizeof(length));
	outVector.resize(length);
	for (int l = 0; l < length; l++) {
		double value;
		file.read((char*)(&value), sizeof(value));
		outVector(l) = value;
	}
	file.close();

	return true;
}

bool DataUtil::saveMatrixBinary(string fileName, Eigen::MatrixXd &inMatrix) {

	ofstream file(fileName, ios::binary);
	if (!file.is_open()) {
		cout << "Error: cannot write to file " << fileName << endl;
		return false;
	}
	int rows = (int)inMatrix.rows();
	int cols = (int)inMatrix.cols();
	file.write((const char*)(&rows), sizeof(rows));
	file.write((const char*)(&cols), sizeof(cols));
	for (int r = 0; r < rows; r++) {
		for (int c = 0; c < cols; c++) {
			double value = inMatrix(r, c);
			file.write((const char*)(&value), sizeof(value));
		}
	}
	file.close();

	return true;
}

bool DataUtil::loadMatrixBinary(string fileName, Eigen::MatrixXd &outMatrix) {

	ifstream file(fileName, ios::binary);
	if (!file.is_open()) {
		cout << "Error: cannot open file " << fileName << endl;
		return false;
	}
	int rows, cols;
	file.read((char*)(&rows), sizeof(rows));
	file.read((char*)(&cols), sizeof(cols));
	outMatrix.resize(rows, cols);
	for (int r = 0; r < rows; r++) {
		for (int c = 0; c < cols; c++) {
			double value;
			file.read((char*)(&value), sizeof(value));
			outMatrix(r, c) = value;
		}
	}
	file.close();

	return true;
}

bool DataUtil::saveCellArrayBinary(string fileName, vector<int> &inArray) {

	ofstream file(fileName, ios::binary);
	if (!file.is_open()) {
		cout << "Error: cannot write to file " << fileName << endl;
		return false;
	}

	int length = (int)inArray.size();
	file.write((const char*)(&length), sizeof(length));
	for (int k = 0; k < length; k++) {
		int value = inArray[k];
		file.write((const char*)(&value), sizeof(value));
	}

	file.close();

	return true;
}

bool DataUtil::loadCellArrayBinary(string fileName, vector<int> &outArray) {

	ifstream file(fileName, ios::binary);
	if (!file.is_open()) {
		cout << "Error: cannot open file " << fileName << endl;
		return false;
	}

	int length;
	file.read((char*)(&length), sizeof(length));
	outArray.resize(length);
	for (int k = 0; k < length; k++) {
		int value;
		file.read((char*)(&value), sizeof(value));
		outArray[k] = value;
	}

	file.close();

	return true;
}

bool DataUtil::saveCellArraysBinary(string fileName, vector<vector<int>> &inArrays) {

	ofstream file(fileName, ios::binary);
	if (!file.is_open()) {
		cout << "Error: cannot write to file " << fileName << endl;
		return false;
	}

	int rows = (int)inArrays.size();
	file.write((const char*)(&rows), sizeof(rows));
	for (int r = 0; r < rows; r++) {
		auto &row = inArrays[r];
		int length = (int)row.size();
		file.write((const char*)(&length), sizeof(length));
		for (int k = 0; k < length; k++) {
			int value = row[k];
			file.write((const char*)(&value), sizeof(value));
		}
	}

	file.close();

	return true;
}

bool DataUtil::loadCellArraysBinary(string fileName, vector<vector<int>> &outArrays) {

	ifstream file(fileName, ios::binary);
	if (!file.is_open()) {
		cout << "Error: cannot open file " << fileName << endl;
		return false;
	}

	int rows;
	file.read((char*)(&rows), sizeof(rows));
	outArrays.resize(rows);
	for (int r = 0; r < rows; r++) {
		auto &row = outArrays[r];
		int length;
		file.read((char*)(&length), sizeof(length));
		row.resize(length);
		for (int k = 0; k < length; k++) {
			int value;
			file.read((char*)(&value), sizeof(value));
			row[k] = value;
		}
	}

	file.close();

	return true;
}