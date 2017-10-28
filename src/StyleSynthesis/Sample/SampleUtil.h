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
#include <set>

#include "Eigen/Eigen"

#include "Library/TheaKDTreeHelper.h"

#include "Data/StyleSynthesisTypes.h"

using namespace std;

namespace StyleSynthesis {

	class SampleUtil {

		struct TPlySample {
			vec3 p; // position
			vec3 n; // normal
			int f; // face ID
		};

	private:

		// make it non-instantiable
		SampleUtil() {}
		~SampleUtil() {}

	public:

		static bool buildKdTree(vector<vec3> &samples, SKDTree &tree, SKDTreeData &treeData);
		static bool buildKdTree(Eigen::Matrix3Xd &samples, SKDTree &tree, SKDTreeData &treeData);
		static bool findNearestNeighbors(SKDTree &tree, Eigen::Matrix3Xd &inPoints, Eigen::VectorXi &outIndices);
		static bool findNearestNeighbors(SKDTree &tree, Eigen::Matrix3Xd &inPoints, set<int> &outIndices, double maxDist = -1);
		static bool sliceMatrices(Eigen::Matrix3Xd &inMatrix, Eigen::VectorXi &inIndices, Eigen::Matrix3Xd &outMatrix);

		static bool computeNormalPCA(TPointSet &inPoints, vec3 &outCenter, vec3 &outNormal);
		static bool computeAABB(TPointSet &inPoints, vec3 &outBBMin, vec3 &outBBMax);
		static bool computeVolume(TSampleSet &inPoints, double &outVolume);
		static bool recomputeNormals(TSampleSet &samples);

		static bool buildMatrices(vector<vec3> &inPoints, Eigen::Matrix3Xd &outMat);
		static bool buildMatrices(TPointSet &inPoints, Eigen::Matrix3Xd &outMat);
		static bool buildMatrices(TPointSet &inPoints, Eigen::Matrix3Xd &outMatP, Eigen::Matrix3Xd &outMatN);
		static bool subSampleMatrices(Eigen::Matrix3Xd &inMat, Eigen::VectorXi &outIdx, int subsample);

		static bool saveSample(string fileName, TSampleSet &samples);
		static bool loadSample(string fileName, TSampleSet &samples);

		static double getPreciseRandomNumber(); // higher precision random number generator
		static int samplingCDF(vector<double> &cdf); // generate random sample based on CDF
	};
}