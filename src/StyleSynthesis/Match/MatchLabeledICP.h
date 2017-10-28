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
#include <string>

#include <Eigen/Eigen>

#include "Library/TheaKDTreeHelper.h"

using namespace std;

namespace StyleSynthesis {

	class MatchLabeledICP {

	private:

		MatchLabeledICP() {}
		~MatchLabeledICP() {}

	public:

		// mode: alignment mode (primitive specific)
		// primitive: 0 - stick; 1 - plane; 2 - sphere
		// sourcePoints: 3xN source point positions
		// targetPoints: 3xN target point positions
		// sourceLabels: Nx1 source point labels
		// targetLabels: Nx1 target point labels
		// orientation: major orientation
		// transformation: output transformation
		// labelWeights: optional Nx1 label weights
		static bool run(
			int iteration,
			int mode,
			int primitive,
			Eigen::Matrix3Xd &sourcePoints,
			Eigen::Matrix3Xd &targetPoints,
			Eigen::VectorXi &sourceLabels,
			Eigen::VectorXi &targetLabels,
			Eigen::Vector3d &orientation,
			Eigen::Affine3d &transformation,
			Eigen::VectorXd *labelWeights = 0);

		// sourcePoints: 3xN source point positions
		// targetPoints: 3xN target point positions
		// sourceLabels: Nx1 source point labels
		// targetLabels: Nx1 target point labels
		// outLabelError: squared distance for closest points for each label
		static bool error(
			Eigen::Matrix3Xd &sourcePoints,
			Eigen::Matrix3Xd &targetPoints,
			Eigen::VectorXi &sourceLabels,
			Eigen::VectorXi &targetLabels,
			Eigen::VectorXd &outLabelError);

		// filename: output .ply file
		// source: untransformed source points
		// target: untransformed target points
		// transformation: ICP result
		static bool visualize(
			string filename,
			Eigen::Matrix3Xd &sourcePoints,
			Eigen::Matrix3Xd &targetPoints,
			Eigen::VectorXi &sourceLabels,
			Eigen::VectorXi &targetLabels,
			Eigen::Affine3d &transformation);

		static bool visualizeLink(
			string filename,
			Eigen::Matrix3Xd &sourcePoints,
			Eigen::Matrix3Xd &targetPoints,
			Eigen::VectorXi &sourceLabels,
			Eigen::VectorXi &targetLabels,
			Eigen::Affine3d &transformation);

		static bool buildMatrices(
			vector<Eigen::Vector3d> &inPoints,
			vector<int> &inLabels,
			Eigen::Matrix3Xd &outMat,
			Eigen::VectorXi &outVec);

	private:

		static bool preAlign(
			Eigen::Matrix3Xd &source,
			Eigen::Matrix3Xd &target,
			Eigen::Vector3d &orientation,
			Eigen::Affine3d &preXform,
			Eigen::Affine3d &postXform);

		static bool fit(
			Eigen::Matrix3Xd &sourcePoints,
			Eigen::Matrix3Xd &targetPoints,
			Eigen::VectorXi &sourceLabels,
			Eigen::VectorXi &targetLabels,
			Eigen::Affine3d &transformation);

		static bool buildKDTrees(
			Eigen::Matrix3Xd &points,
			Eigen::VectorXi &labels,
			vector<SKDTree*> &trees,
			vector<SKDTreeData> &treesData);

		static bool findMatchedNeighbors(
			Eigen::Matrix3Xd &inoutSourcePoints,
			Eigen::VectorXi &inoutSourceLabels,
			Eigen::Matrix3Xd &inoutTargetPoints,
			Eigen::VectorXi &inoutTargetLabels,
			Eigen::VectorXd &outPointWeights,
			Eigen::VectorXd *labelWeights = 0);

		static bool findNearestNeighbors(
			vector<SKDTree*> &trees,
			Eigen::Matrix3Xd &inPoints,
			Eigen::VectorXi &inLabels,
			vector<int> &outIndices);

		static bool findNearestByLabel(
			Eigen::Matrix3Xd &inPoints,
			Eigen::Matrix3Xd &slicedPoints,
			Eigen::VectorXi &inLabels,
			vector<int> &outIndices);

		static bool findFarthestByLabel(
			Eigen::Matrix3Xd &inPoints,
			Eigen::Matrix3Xd &slicedPoints,
			Eigen::VectorXi &inLabels,
			vector<int> &outIndices);

		static bool sliceMatrices(
			vector<int> &inIndices,
			Eigen::Matrix3Xd &inoutPoints,
			Eigen::VectorXi &inoutLabels);

		static bool extractTransformation(
			int mode,
			int primitive,
			Eigen::Matrix3Xd &source,
			Eigen::Matrix3Xd &target,
			Eigen::VectorXd &weights,
			Eigen::Affine3d &transformation);

		static bool extractPlanarRotation(
			Eigen::Matrix3Xd &source,
			Eigen::Matrix3Xd &target,
			Eigen::VectorXd &weights,
			Eigen::Matrix3d &rotation);

		static bool extractFreeRotation(
			Eigen::Matrix3Xd &source,
			Eigen::Matrix3Xd &target,
			Eigen::VectorXd &weights,
			Eigen::Matrix3d &rotation);
	};
}