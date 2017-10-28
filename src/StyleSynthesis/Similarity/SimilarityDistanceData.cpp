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

#include "SimilarityDistanceData.h"

#include "Data/DataUtil.h"

#include <fstream>

using namespace std;
using namespace StyleSynthesis;

SimilarityDistanceData::SimilarityDistanceData() {
}

SimilarityDistanceData::~SimilarityDistanceData() {
}

bool SimilarityDistanceData::saveData(string folderName) {

	// names

	string dataIndicesName = folderName + "data-indices.txt";
	string dataDistanceName = folderName + "data-distance.txt";
	string dataSourceElementName = folderName + "data-source-element.txt";
	string dataTargetElementName = folderName + "data-target-element.txt";
	string dataSourceUnmatchName = folderName + "data-source-unmatch.txt";
	string dataTargetUnmatchName = folderName + "data-target-unmatch.txt";

	// output

	if (!DataUtil::savePairListASCII(dataIndicesName, mElementIndices)) return false;
	if (!DataUtil::saveMatrixBinary(dataDistanceName, mElementDistance)) return false;
	if (!DataUtil::saveCellArraysBinary(dataSourceElementName, mElementSourcePoints)) return false;
	if (!DataUtil::saveCellArraysBinary(dataTargetElementName, mElementTargetPoints)) return false;
	if (!DataUtil::saveCellArrayBinary(dataSourceUnmatchName, mUnmatchSourcePoints)) return false;
	if (!DataUtil::saveCellArrayBinary(dataTargetUnmatchName, mUnmatchTargetPoints)) return false;

	return true;
}

bool SimilarityDistanceData::loadData(string folderName) {

	// names

	string dataIndicesName = folderName + "data-indices.txt";
	string dataDistanceName = folderName + "data-distance.txt";
	string dataSourceElementName = folderName + "data-source-element.txt";
	string dataTargetElementName = folderName + "data-target-element.txt";
	string dataSourceUnmatchName = folderName + "data-source-unmatch.txt";
	string dataTargetUnmatchName = folderName + "data-target-unmatch.txt";

	// input

	if (!DataUtil::loadPairListASCII(dataIndicesName, mElementIndices)) return false;
	if (!DataUtil::loadMatrixBinary(dataDistanceName, mElementDistance)) return false;
	if (!DataUtil::loadCellArraysBinary(dataSourceElementName, mElementSourcePoints)) return false;
	if (!DataUtil::loadCellArraysBinary(dataTargetElementName, mElementTargetPoints)) return false;
	if (!DataUtil::loadCellArrayBinary(dataSourceUnmatchName, mUnmatchSourcePoints)) return false;
	if (!DataUtil::loadCellArrayBinary(dataTargetUnmatchName, mUnmatchTargetPoints)) return false;

	return true;
}