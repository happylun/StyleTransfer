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

	class MatchPrimitiveICP {

	private:

		MatchPrimitiveICP() {}
		~MatchPrimitiveICP() {}

	public:

		// mode: alignment mode (primitive specific)
		// primitive: 0 - stick; 1 - plane; 2 - sphere
		// source: 3xN source point positions
		// target: 3xN target point positions
		// orientation: major orientation
		// transformation: output transformation
		static bool run(
			int iteration,
			int mode,
			int primitive,
			Eigen::Matrix3Xd &source,
			Eigen::Matrix3Xd &target,
			Eigen::Vector3d &orientation,
			Eigen::Affine3d &transformation);

		// mode: alignment mode (primitive specific)
		// primitive: 0 - stick; 1 - plane; 2 - sphere
		// source: 3xN source point positions
		// target: 3xN target point positions
		// orientation: major orientation
		// transformation: output transformation
		static bool fit(
			int mode,
			int primitive,
			Eigen::Matrix3Xd &source,
			Eigen::Matrix3Xd &target,
			Eigen::Vector3d &orientation,
			Eigen::Affine3d &transformation);

		// source: 3xNs source point sets
		// target: 3xNt target point sets
		// error: average squared distance to closest point
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

		static bool regularizeAxis(
			Eigen::Vector3d &axis);

	private:

		static bool preAlign(
			Eigen::Matrix3Xd &source,
			Eigen::Matrix3Xd &target,
			Eigen::Vector3d &orientation,
			Eigen::Affine3d &preXform,
			Eigen::Affine3d &postXform);

		static bool buildKDTree(
			Eigen::Matrix3Xd &points,
			SKDTree &tree,
			SKDTreeData &data);

		static bool findNearestNeighbors(
			SKDTree &tree,
			Eigen::Matrix3Xd &inPoints,
			vector<int> &outIndices);

		static bool findMatchedNeighbors(
			Eigen::Matrix3Xd &inSource,
			Eigen::Matrix3Xd &inTarget,
			vector<int> &outIndices);

		static bool sliceMatrices(
			Eigen::Matrix3Xd &inMatrix,
			vector<int> &inIndices,
			Eigen::Matrix3Xd &outMatrix);

		static bool extractTransformation(
			int mode,
			int primitive,
			Eigen::Matrix3Xd &source,
			Eigen::Matrix3Xd &target,
			Eigen::Affine3d &transformation);

		static bool extractFittingMass(
			int mode,
			int primitive,
			Eigen::Matrix3Xd &source,
			Eigen::Matrix3Xd &target,
			Eigen::Affine3d &transformation);

		static bool extractFittingBoundingBox(
			int mode,
			int primitive,
			Eigen::Matrix3Xd &source,
			Eigen::Matrix3Xd &target,
			Eigen::Affine3d &transformation);

		static bool extractPlanarRotation(
			Eigen::Matrix3Xd &source,
			Eigen::Matrix3Xd &target,
			Eigen::Matrix3d &rotation);

		static bool extractFreeRotation(
			Eigen::Matrix3Xd &source,
			Eigen::Matrix3Xd &target,
			Eigen::Matrix3d &rotation);
	};
}