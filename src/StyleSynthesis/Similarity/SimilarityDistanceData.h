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

#include <vector>

#include <Eigen/Dense>

#include "Library/CMLHelper.h"

using namespace std;

namespace StyleSynthesis {

	class SimilarityDistanceData {

		friend class SimilarityDistance;
		friend class SimilarityMetric;

	public:

		SimilarityDistanceData();
		~SimilarityDistanceData();

	public:

		bool saveData(string folderName);
		bool loadData(string folderName);

	protected:

		vector<vec2i> mElementIndices; // (source segment ID, target segment ID) : # of element pairs
		vector<vector<int>> mElementSourcePoints; // point ID : # of points on source element : # of element pairs
		vector<vector<int>> mElementTargetPoints; // point ID : # of points on target element : # of element pairs
		vector<int> mUnmatchSourcePoints; // point ID : # of unmatched points on source shape
		vector<int> mUnmatchTargetPoints; // point ID : # of unmatched points on target shape

		Eigen::MatrixXd mElementDistance; // # of element pairs X distance dimension
	};

}