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

#include <Eigen/Eigen>

#include "Library/TheaKDTreeHelper.h"

using namespace std;

namespace StyleSynthesis {

	class MatchCurveICP {

	private:

		MatchCurveICP() {}
		~MatchCurveICP() {}

	public:

		// source: 3xN source point positions
		// target: 3xN target point positions
		// transformation: aligned transformation (translation and scaling only)
		static bool prealign(
			Eigen::Matrix3Xd &source,
			Eigen::Matrix3Xd &target,
			Eigen::Affine3d &transformation);

		// sourceP: 3xN source point positions
		// sourceN: 3xN source point normals
		// targetP: 3xN target point positions
		// targetN: 3xN target point normals
		// transformation: initial and output transformation
		// success: output success flag
		// tree: optional kd tree for target points
		static bool run(
			int iteration,
			Eigen::Matrix3Xd &sourceP,
			Eigen::Matrix3Xd &sourceN,
			Eigen::Matrix3Xd &targetP,
			Eigen::Matrix3Xd &targetN,
			Eigen::Affine3d &transformation,
			bool &success,
			SKDTree *tree = 0);

		// inSource: 3xN source points
		// inTarget: 3xN target points
		// outSource: output source inliner points
		// outTarget: output target inliner points
		// transformation: ICP result
		// tree: optional kd tree for target points
		static bool extractInliner(
			Eigen::Matrix3Xd &inSource,
			Eigen::Matrix3Xd &inTarget,
			Eigen::Matrix3Xd &outSource,
			Eigen::Matrix3Xd &outTarget,
			Eigen::Affine3d &transformation,
			SKDTree *tree = 0);

		// inSourcePoint1/2: two 3D points
		// inTargetPoint1/2: two 3D points
		// transformation: output transformation to align two source points with two target points
		static bool alignTwoPoints(
			Eigen::Vector3d inSourcePoint1, Eigen::Vector3d inSourcePoint2,
			Eigen::Vector3d inTargetPoint1, Eigen::Vector3d inTargetPoint2,
			Eigen::Affine3d &transformation);

		// inSource: 3xN source points
		// inTarget: 3xN target points
		// outCorrespondences: source ID : # of target points
		static bool findCorrespondences(
			Eigen::Matrix3Xd &inSource,
			Eigen::Matrix3Xd &inTarget,
			Eigen::VectorXi &outCorrespondences);

		// source: 3xNs source point sets
		// target: 3xNt target point sets
		// distance: average squared distance to closest point
		static bool error(
			Eigen::Matrix3Xd &source,
			Eigen::Matrix3Xd &target,
			double &outError);

		// filename: output .ply file
		// source: untransformed source point positions
		// target: target point positions
		// transformation: ICP result
		static bool visualize(
			string filename,
			Eigen::Matrix3Xd &source,
			Eigen::Matrix3Xd &target,
			Eigen::Affine3d &transformation);

	private:

		static bool buildKDTree(
			Eigen::Matrix3Xd &points,
			SKDTree &tree,
			SKDTreeData &data);

		static bool findNearestNeighbors(
			SKDTree &tree,
			Eigen::Matrix3Xd &inPoints,
			vector<int> &outIndices);

		static bool findMatchedNeighbors(
			Eigen::Matrix3Xd &inSourceP,
			Eigen::Matrix3Xd &inSourceN,
			Eigen::Matrix3Xd &inTargetP,
			Eigen::Matrix3Xd &inTargetN,
			vector<int> &outIndices);

		static bool sliceMatrices(
			Eigen::Matrix3Xd &inMatrix,
			vector<int> &inIndices,
			Eigen::Matrix3Xd &outMatrix);

		static bool extractTransformation(
			Eigen::Matrix3Xd &source,
			Eigen::Matrix3Xd &target,
			Eigen::Affine3d &transformation);
	};
}