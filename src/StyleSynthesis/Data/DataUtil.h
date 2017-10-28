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

#pragma once

#include <string>
#include <vector>

#include <Eigen/Dense>

#include "Library/CMLHelper.h"

using namespace std;

namespace StyleSynthesis {

	class DataUtil {

	private:

		// make it non-instantiable
		DataUtil() {}
		~DataUtil() {}

	public:

		static bool saveVectorASCII(string fileName, Eigen::VectorXd &inVector);
		static bool loadVectorASCII(string fileName, Eigen::VectorXd &outVector);

		static bool saveRowVectorASCII(string fileName, Eigen::RowVectorXd &inVector);
		static bool loadRowVectorASCII(string fileName, Eigen::RowVectorXd &outVector);

		static bool saveMatrixASCII(string fileName, Eigen::MatrixXd &inMatrix);
		static bool loadMatrixASCII(string fileName, Eigen::MatrixXd &outMatrix);

		static bool savePairListASCII(string fileName, vector<vec2i> &inPairList);
		static bool loadPairListASCII(string fileName, vector<vec2i> &outPairList);

		static bool saveIndexListASCII(string fileName, vector<int> &inIndexList);
		static bool loadIndexListASCII(string fileName, vector<int> &outIndexList);

		static bool saveValueListASCII(string fileName, vector<double> &inValueList);
		static bool loadValueListASCII(string fileName, vector<double> &outValueList);

		static bool saveGroupListASCII(string fileName, vector<vector<int>> &inGroupList);
		static bool loadGroupListASCII(string fileName, vector<vector<int>> &outGroupList);

		static bool saveVectorBinary(string fileName, Eigen::VectorXd &inVector);
		static bool loadVectorBinary(string fileName, Eigen::VectorXd &outVector);

		static bool saveRowVectorBinary(string fileName, Eigen::RowVectorXd &inVector);
		static bool loadRowVectorBinary(string fileName, Eigen::RowVectorXd &outVector);

		static bool saveMatrixBinary(string fileName, Eigen::MatrixXd &inMatrix);
		static bool loadMatrixBinary(string fileName, Eigen::MatrixXd &outMatrix);

		static bool saveCellArrayBinary(string fileName, vector<int> &inArray);
		static bool loadCellArrayBinary(string fileName, vector<int> &outArray);

		static bool saveCellArraysBinary(string fileName, vector<vector<int>> &inArrays);
		static bool loadCellArraysBinary(string fileName, vector<vector<int>> &outArrays);
	};
}