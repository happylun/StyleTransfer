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

#include "MatchUniformICP.h"

#include <vector>
#include <set>

#include "Utility/PlyExporter.h"

#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

bool MatchUniformICP::runRegularShape(
	int iteration,
	Eigen::Matrix3Xd &sourceP,
	Eigen::Matrix3Xd &sourceN,
	Eigen::Matrix3Xd &targetP,
	Eigen::Matrix3Xd &targetN,
	Eigen::Affine3d &transformation)
{
	// align two shapes with regular initial transformations
	// will initialize transformations using 4 horizontal rotation (w/ or w/o horizontal flipping) and vertical flipping

	// initial transformations
	vector<Eigen::Affine3d> transformations;
	for (int rotID = 0; rotID < 4; rotID++) {
		Eigen::AngleAxisd rotation(cml::rad(90.0*rotID), Eigen::Vector3d(0, 1, 0));
		Eigen::Matrix3d hflip;
		hflip << -1, 0, 0, 0, 1, 0, 0, 0, 1;
		Eigen::Affine3d xform;
		// only rotation
		xform.setIdentity();
		xform.prerotate(rotation);
		transformations.push_back(xform);
		// horizontal flipping + rotation
		xform.setIdentity();
		xform.linear() = rotation * hflip;
		transformations.push_back(xform);
	}
	if (true) {
		// vertical flipping
		Eigen::Matrix3d vflip;
		vflip << 1, 0, 0, 0, -1, 0, 0, 0, 1;
		Eigen::Affine3d xform;
		xform.setIdentity();
		xform.linear() = vflip;
		transformations.push_back(xform);
	}

	transformation.setIdentity();
	double minError = DBL_MAX;

	for (int transID = 0; transID < (int)transformations.size(); transID++) {

		Eigen::Affine3d xform = transformations[transID];
		if (!run(iteration, sourceP, sourceN, targetP, targetN, xform)) return false;

		double currentError;
		Eigen::Matrix3Xd sourceXP = xform * sourceP;
		if (!MatchUniformICP::error(targetP, sourceXP, currentError)) return false; // NOTE: swap source and target

		if (currentError < minError) {
			minError = currentError;
			transformation = xform;
		}
	}

	return true;
}

bool MatchUniformICP::run(
	int iteration,
	Eigen::Matrix3Xd &sourceP,
	Eigen::Matrix3Xd &sourceN,
	Eigen::Matrix3Xd &targetP,
	Eigen::Matrix3Xd &targetN,
	Eigen::Affine3d &transformation)
{

	// initialization
	SKDTree tree;
	SKDTreeData treeData;
	if (!buildKDTree(targetP, tree, treeData)) return false;
	Eigen::Affine3d initialTransformation = transformation;

	// ICP iteration
	for (int iterID = 0; iterID<iteration; iterID++) {

		// find nearest neighbors
		Eigen::Matrix3d rotation = transformation.rotation();
		Eigen::Matrix3Xd matXSP = transformation * sourceP;
		Eigen::Matrix3Xd matXSN = rotation * sourceN;
		Eigen::Matrix3Xd matTP, matTN;
		vector<int> slices;
		if (!findNearestNeighbors(tree, matXSP, slices)) return false;
		if (!sliceMatrices(targetP, slices, matTP)) return false;
		if (!sliceMatrices(targetN, slices, matTN)) return false;
		if (!findMatchedNeighbors(matXSP, matXSN, matTP, matTN, slices)) return false;
		if (slices.empty()) break; // no matched points
		if (!sliceMatrices(matXSP, slices, matXSP)) return false;
		if (!sliceMatrices(matXSN, slices, matXSN)) return false;
		if (!sliceMatrices(matTP, slices, matTP)) return false;
		if (!sliceMatrices(matTN, slices, matTN)) return false;

		// align matched points
		Eigen::Affine3d newTransformation;
		if (!extractTransformation(matXSP, matTP, newTransformation)) return false;
		transformation = newTransformation * transformation;
	}

	if (!transformation.matrix().allFinite()) {
		transformation = initialTransformation;
	}

	return true;
}

bool MatchUniformICP::error(
	Eigen::Matrix3Xd &source,
	Eigen::Matrix3Xd &target,
	double &outError)
{
	SKDTree tree;
	SKDTreeData treeData;
	if (!buildKDTree(target, tree, treeData)) return false;

	Eigen::Matrix3Xd matched;
	vector<int> neighbors;
	if (!findNearestNeighbors(tree, source, neighbors)) return false;
	if (!sliceMatrices(target, neighbors, matched)) return false;
	outError = (matched - source).squaredNorm() / source.cols();

	return true;
}

bool MatchUniformICP::buildKDTree(
	Eigen::Matrix3Xd &points,
	SKDTree &tree,
	SKDTreeData &data)
{

	data.resize(points.cols());
	for (int i = 0; i<(int)points.cols(); i++) {
		data[i] = SKDT::NamedPoint((float)points(0, i), (float)points(1, i), (float)points(2, i), (size_t)i);
	}
	tree.init(data.begin(), data.end());

	return true;
}

bool MatchUniformICP::findNearestNeighbors(
	SKDTree &tree,
	Eigen::Matrix3Xd &inPoints,
	vector<int> &outIndices)
{

	outIndices.resize(inPoints.cols());
#pragma omp parallel for
	for (int j = 0; j < inPoints.cols(); j++) {
		SKDT::NamedPoint queryPoint((float)inPoints(0, j), (float)inPoints(1, j), (float)inPoints(2, j));
		Thea::BoundedSortedArray<SKDTree::Neighbor> queryResult(1);
		tree.kClosestElements<Thea::MetricL2>(queryPoint, queryResult);
		if (queryResult.isEmpty()) outIndices[j] = 0; // too far away, whatever...
		else outIndices[j] = (int)tree.getElements()[queryResult[0].getIndex()].id;
	}

	return true;
}

bool MatchUniformICP::findMatchedNeighbors(
	Eigen::Matrix3Xd &inSourceP,
	Eigen::Matrix3Xd &inSourceN,
	Eigen::Matrix3Xd &inTargetP,
	Eigen::Matrix3Xd &inTargetN,
	vector<int> &outIndices)
{
	if (inSourceP.cols() == 0) return true; // empty

	const double rejectDistanceThreshold = 5.0;

	Eigen::ArrayXd vecD = (inSourceP - inTargetP).colwise().norm().array();
	Eigen::ArrayXd vecN = (inSourceN.transpose() * inTargetN).diagonal().array();
	double maxDist;
	if (true) {
		// use r times median length as clamping distance (r = 5?)
		vector<double> vDist(vecD.data(), vecD.data() + vecD.size());
		nth_element(vDist.begin(), vDist.begin() + vecD.size() / 2, vDist.end());
		maxDist = vDist[vecD.size() / 2] * rejectDistanceThreshold;
	}

	auto filter = vecD < maxDist && vecN > 0.0;
	outIndices.clear();
	outIndices.reserve((int)filter.count());
	for (int j = 0; j < filter.size(); j++) { // don't parallelize
		if (filter(j)) {
			outIndices.push_back(j);
		}
	}

	return true;
}

bool MatchUniformICP::sliceMatrices(
	Eigen::Matrix3Xd &inMatrix,
	vector<int> &inIndices,
	Eigen::Matrix3Xd &outMatrix)
{
	Eigen::Matrix3Xd tmpMatrix;
	tmpMatrix.resize(inMatrix.rows(), inIndices.size());
	for (int j = 0; j < (int)inIndices.size(); j++) {
		tmpMatrix.col(j) = inMatrix.col(inIndices[j]);
	}
	outMatrix.swap(tmpMatrix);

	return true;
}

bool MatchUniformICP::extractTransformation(
	Eigen::Matrix3Xd &source,
	Eigen::Matrix3Xd &target,
	Eigen::Affine3d &transformation)
{
	// classical ICP solution

	Eigen::Vector3d vecSCenter = source.rowwise().mean();
	Eigen::Vector3d vecTCenter = target.rowwise().mean();
	Eigen::Matrix3Xd transSource = source.colwise() - vecSCenter;
	Eigen::Matrix3Xd transTarget = target.colwise() - vecTCenter;
	double scaSScale = transSource.norm();
	double scaTScale = transTarget.norm();
	double scaScale = scaTScale / scaSScale;
	transSource *= scaScale;
	Eigen::JacobiSVD<Eigen::Matrix3d> svd(transTarget*transSource.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix3d matRotate = svd.matrixU() * svd.matrixV().transpose();
	Eigen::Vector3d vecTranslate = vecTCenter - matRotate * vecSCenter * scaScale;
	transformation.setIdentity();
	transformation.prescale(scaScale);
	transformation.prerotate(matRotate);
	transformation.pretranslate(vecTranslate);
	
	return true;
}

bool MatchUniformICP::visualize(
	string filename,
	Eigen::Matrix3Xd &source,
	Eigen::Matrix3Xd &target,
	Eigen::Affine3d &transformation)
{

	vector<vec3> sourcePoints(source.cols());
	vector<vec3> targetPoints(target.cols());

	for (int j = 0; j < source.cols(); j++) {
		sourcePoints[j] = vec3d(source(0, j), source(1, j), source(2, j));
	}
	for (int j = 0; j < target.cols(); j++) {
		targetPoints[j] = vec3d(target(0, j), target(1, j), target(2, j));
	}

	matrix4d mat(
		transformation(0, 0), transformation(0, 1), transformation(0, 2), transformation(0, 3),
		transformation(1, 0), transformation(1, 1), transformation(1, 2), transformation(1, 3),
		transformation(2, 0), transformation(2, 1), transformation(2, 2), transformation(2, 3),
		0.0, 0.0, 0.0, 1.0);

	PlyExporter pe;
	if (!pe.addPoint(&sourcePoints, 0, mat, vec3i(255, 0, 0))) return false;
	if (!pe.addPoint(&targetPoints, 0, cml::identity_4x4(), vec3i(0, 255, 0))) return false;
	if (!pe.output(filename)) return false;

	return true;
}