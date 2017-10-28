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

#include "MatchCurveICP.h"

#include <vector>
#include <set>

#include "Sample/SampleUtil.h"

#include "Utility/PlyExporter.h"

#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

bool MatchCurveICP::prealign(
	Eigen::Matrix3Xd &source,
	Eigen::Matrix3Xd &target,
	Eigen::Affine3d &transformation)
{
	Eigen::Vector3d srcCenter = source.rowwise().mean();
	Eigen::Vector3d tgtCenter = target.rowwise().mean();

	int srcDim, tgtDim;
	Eigen::AlignedBox3d srcBB(source.rowwise().minCoeff(), source.rowwise().maxCoeff());
	Eigen::AlignedBox3d tgtBB(target.rowwise().minCoeff(), target.rowwise().maxCoeff());
	double scaSScale = srcBB.sizes().maxCoeff(&srcDim);
	double scaTScale = tgtBB.sizes().maxCoeff(&tgtDim);
	double scale = scaTScale / scaSScale;
	//srcCenter[srcDim] = (srcCenter[srcDim] + srcBB.center()[srcDim]) / 2;
	//tgtCenter[tgtDim] = (tgtCenter[tgtDim] + tgtBB.center()[tgtDim]) / 2;
	srcCenter[srcDim] = srcBB.center()[srcDim];
	tgtCenter[tgtDim] = tgtBB.center()[tgtDim];

	// put everything together
	transformation.setIdentity();
	transformation.pretranslate(-srcCenter);
	transformation.prescale(scale);
	transformation.pretranslate(tgtCenter);

	return true;
}

bool MatchCurveICP::run(
	int iteration,
	Eigen::Matrix3Xd &sourceP,
	Eigen::Matrix3Xd &sourceN,
	Eigen::Matrix3Xd &targetP,
	Eigen::Matrix3Xd &targetN,
	Eigen::Affine3d &transformation,
	bool &success,
	SKDTree *tree)
{

	// initialization
	SKDTree tmpTree;
	SKDTreeData tmpTreeData;
	if (!tree) {
		if (!buildKDTree(targetP, tmpTree, tmpTreeData)) return false;
		tree = &tmpTree;
	}
	Eigen::Affine3d initialTransformation = transformation;

	// ICP iteration
	success = true;
	for (int iterID = 0; iterID<iteration; iterID++) {

		// find nearest neighbors
		Eigen::Matrix3d rotation = transformation.rotation();
		Eigen::Matrix3Xd matXSP = transformation * sourceP;
		Eigen::Matrix3Xd matXSN = rotation * sourceN;
		Eigen::Matrix3Xd matTP, matTN;

		vector<int> slices;
		if (!findNearestNeighbors(*tree, matXSP, slices)) return false;
		if (!sliceMatrices(targetP, slices, matTP)) return false;
		if (!sliceMatrices(targetN, slices, matTN)) return false;
		/* // HACK: no filtering
		if (!findMatchedNeighbors(matXSP, matXSN, matTP, matTN, slices)) return false;
		if (slices.size() < 20) { // UNDONE: param minimum matched point
			success = false;
			break;
		}
		if (!sliceMatrices(matXSP, slices, matXSP)) return false;
		if (!sliceMatrices(matXSN, slices, matXSN)) return false;
		if (!sliceMatrices(matTP, slices, matTP)) return false;
		if (!sliceMatrices(matTN, slices, matTN)) return false;
		*/
		// align matched points
		Eigen::Affine3d newTransformation;
		if (!extractTransformation(matXSP, matTP, newTransformation)) return false;
		transformation = newTransformation * transformation;

		/*
		// visualize each iteration
		if (!visualize("iter-match.ply", matXSP, matTP, newTransformation)) return false;
		if (!visualize("iter-all.ply", sourceP, targetP, transformation)) return false;
		system("pause");
		//*/
	}

	if (!transformation.matrix().allFinite()) {
		transformation = initialTransformation;
	}

	return true;
}

bool MatchCurveICP::extractInliner(
	Eigen::Matrix3Xd &inSource,
	Eigen::Matrix3Xd &inTarget,
	Eigen::Matrix3Xd &outSource,
	Eigen::Matrix3Xd &outTarget,
	Eigen::Affine3d &transformation,
	SKDTree *tree)
{
	// initialization

	SKDTree tmpTree;
	SKDTreeData tmpTreeData;
	if (!tree) {
		if (!buildKDTree(inTarget, tmpTree, tmpTreeData)) return false;
		tree = &tmpTree;
	}

	outSource.resize(3, 0);
	outTarget.resize(3, 0);

	// find inliner slice
	
	vector<int> slices;
	Eigen::Matrix3Xd nnSource = transformation * inSource;
	Eigen::Matrix3Xd nnTarget; // have the same size with nnSource
	if (true) {
		if (!findNearestNeighbors(*tree, nnSource, slices)) return false;
		if (!sliceMatrices(inTarget, slices, nnTarget)) return false;

		// fake normals - just used to pass the matching criteria
		Eigen::Matrix3Xd tmpSrcNormal = Eigen::Matrix3Xd::Zero(3, nnSource.cols());
		Eigen::Matrix3Xd tmpTgtNormal = Eigen::Matrix3Xd::Zero(3, nnTarget.cols());
		tmpSrcNormal.row(1).setOnes();
		tmpTgtNormal.row(1).setOnes();

		if (!findMatchedNeighbors(nnSource, tmpSrcNormal, nnTarget, tmpTgtNormal, slices)) return false;
		if (slices.empty()) return true; // no matched points
	}

	// get longest section from slice

	vector<int> longestSection(0);
	if (true) {
		vector<int> currentSection(1, slices[0]);
		for (int sliceID = 1; sliceID < (int)slices.size(); sliceID++) {
			if (slices[sliceID] == currentSection.back() + 1) {
				currentSection.push_back(slices[sliceID]);
			} else {
				if (currentSection.size() > longestSection.size()) longestSection.swap(currentSection);
				currentSection.assign(1, slices[sliceID]);
			}
		}
		if (currentSection.size() > longestSection.size()) longestSection.swap(currentSection);
	}

	if (!sliceMatrices(inSource, longestSection, outSource)) return false;
	if (!sliceMatrices(nnTarget, longestSection, outTarget)) return false;

	return true;
}

bool MatchCurveICP::alignTwoPoints(
	Eigen::Vector3d inSourcePoint1, Eigen::Vector3d inSourcePoint2,
	Eigen::Vector3d inTargetPoint1, Eigen::Vector3d inTargetPoint2,
	Eigen::Affine3d &transformation)
{
	double dist11 = (inSourcePoint1 - inTargetPoint1).norm();
	double dist22 = (inSourcePoint2 - inTargetPoint2).norm();
	if (dist22 < dist11) {
		swap(inSourcePoint1, inSourcePoint2);
		swap(inTargetPoint1, inTargetPoint2);
	}

	Eigen::Vector3d vecS = inSourcePoint2 - inSourcePoint1;
	Eigen::Vector3d vecT = inTargetPoint2 - inTargetPoint1;
	Eigen::Quaterniond rotation = Eigen::Quaterniond::FromTwoVectors(vecS, vecT);
	double scale = vecT.norm() / vecS.norm();

	transformation.setIdentity();
	transformation.pretranslate(-inSourcePoint1);
	transformation.prescale(scale);
	transformation.prerotate(rotation);
	transformation.pretranslate(inTargetPoint1);

	return true;
}

bool MatchCurveICP::findCorrespondences(
	Eigen::Matrix3Xd &inSource,
	Eigen::Matrix3Xd &inTarget,
	Eigen::VectorXi &outCorrespondences)
{

	int numSteps = 10; // UNDONE: param iterative alignment steps

	int numSourcePoints = (int)inSource.cols();
	int numTargetPoints = (int)inTarget.cols();

	// build KD tree

	SKDTree tree;
	SKDTreeData treeData;
	if (!SampleUtil::buildKdTree(inSource, tree, treeData)) return false;

	// progressive alignment
	
	Eigen::Matrix3Xd movingPoints = inTarget;
	for (int step = 0; step < numSteps; step++) {
		Eigen::Matrix3Xd matchedPoints;
		if (!SampleUtil::findNearestNeighbors(tree, movingPoints, outCorrespondences)) return false;

		Eigen::Matrix3Xd sourceOffset(3, numSourcePoints);
		sourceOffset.setZero();
		vector<bool> sourceUsedFlags(numSourcePoints, false);
		for (int tgtID = 0; tgtID < numTargetPoints; tgtID++) {
			int srcID = outCorrespondences[tgtID];
			Eigen::Vector3d newOffset = inSource.col(srcID) - movingPoints.col(tgtID);
			if (!sourceUsedFlags[srcID] || newOffset.norm() < sourceOffset.col(srcID).norm()) {
				sourceOffset.col(srcID) = newOffset;
				sourceUsedFlags[srcID] = true;
			}
		}
		sourceOffset /= (double)(numSteps - step);

		for (int tgtID = 0; tgtID < numTargetPoints; tgtID++) {
			int srcID = outCorrespondences[tgtID];
			movingPoints.col(tgtID) += sourceOffset.col(srcID);
		}
	}

	// HACK: fill gaps

	if (true) {

		// get match set for source points

		vector<vector<int>> sourceMatchSet(numSourcePoints, vector<int>(0));
		for (int tgtID = 0; tgtID < numTargetPoints; tgtID++) {
			int srcID = outCorrespondences[tgtID];
			sourceMatchSet[srcID].push_back(tgtID);
		}

		// find source match set pivot (nearest point)

		vector<int> sourceMatchPivot(numSourcePoints, -1);
		for (int srcID = 0; srcID < numSourcePoints; srcID++) {
			int numMatch = (int)sourceMatchSet[srcID].size();
			double minDist = DBL_MAX;
			for (int matchID = 0; matchID < numMatch; matchID++) {
				int tgtID = sourceMatchSet[srcID][matchID];
				double dist = (inSource.col(srcID) - inTarget.col(tgtID)).norm();
				if (dist < minDist) {
					minDist = dist;
					sourceMatchPivot[srcID] = matchID;
				}
			}
		}

		// process each gap

		vector<int> targetPool(0);
		vector<int> sourceProcessedFlags(numSourcePoints, false);
		int lastProcessedPoint = -1;
		for (int tgtID = 0; tgtID <= numTargetPoints; tgtID++) {
			int srcID;
			int pivot;
			if (tgtID < numTargetPoints) {
				srcID = outCorrespondences[tgtID];
				if (sourceProcessedFlags[srcID]) continue;
				if (sourceMatchSet[srcID].empty()) continue;

				// add candidate to target pool
				vector<int> &matchSet = sourceMatchSet[srcID];
				pivot = sourceMatchPivot[srcID];
				for (int k = 0; k < pivot; k++) targetPool.push_back(matchSet[k]);
			} else {
				srcID = numSourcePoints;
			}

			// build source pool
			vector<int> sourcePool;
			for (int id = lastProcessedPoint + 1; id < srcID; id++) sourcePool.push_back(id);

			// match gap
			if (!sourcePool.empty() && !targetPool.empty()) {
				int srcSize = (int)sourcePool.size();
				int tgtSize = (int)targetPool.size();
				int srcPtr = 0;
				int tgtPtr = 0;
				while (srcPtr < srcSize && tgtPtr < tgtSize) {
					outCorrespondences[targetPool[tgtPtr]] = sourcePool[srcPtr];
					double srcPrg = (srcPtr + 1) / (double)srcSize;
					double tgtPrg = (tgtPtr + 1) / (double)tgtSize;
					if (srcPrg < tgtPrg) srcPtr++;
					else tgtPtr++;
				}
				if (srcPtr + 1 < srcSize || tgtPtr + 1 < tgtSize) {
					cout << "Error: failed to fill gap" << endl;
					return false;
				}
			}

			if (tgtID < numTargetPoints) {
				// update target pool
				targetPool.clear();
				vector<int> &matchSet = sourceMatchSet[srcID];
				for (int k = pivot + 1; k < (int)matchSet.size(); k++) targetPool.push_back(matchSet[k]);

				sourceProcessedFlags[srcID] = true;
				lastProcessedPoint = max(lastProcessedPoint, srcID);
			}
		}
	}


	return true;
}

bool MatchCurveICP::error(
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

bool MatchCurveICP::buildKDTree(
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

bool MatchCurveICP::findNearestNeighbors(
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

bool MatchCurveICP::findMatchedNeighbors(
	Eigen::Matrix3Xd &inSourceP,
	Eigen::Matrix3Xd &inSourceN,
	Eigen::Matrix3Xd &inTargetP,
	Eigen::Matrix3Xd &inTargetN,
	vector<int> &outIndices)
{
	if (inSourceP.cols() == 0) return true; // empty

	const double rejectDistanceThreshold = 5.0;

	Eigen::ArrayXd vecD = (inSourceP - inTargetP).colwise().norm().array();
	Eigen::ArrayXd vecN = (inSourceN.transpose() * inTargetN).diagonal().array().abs(); // only care about unsigned angle
	double maxDist;
	if (true) {
		// use r times median length as clamping distance (r = 5?)
		vector<double> vDist(vecD.data(), vecD.data() + vecD.size());
		nth_element(vDist.begin(), vDist.begin() + vecD.size() / 2, vDist.end());
		maxDist = vDist[vecD.size() / 2] * rejectDistanceThreshold;
	}

	auto filter = vecD < maxDist && vecN > 0.5; // UNDONE: matching threshold for |cosA|
	outIndices.clear();
	outIndices.reserve((int)filter.count());
	for (int j = 0; j < filter.size(); j++) { // don't parallelize
		if (filter(j)) {
			outIndices.push_back(j);
		}
	}

	return true;
}

bool MatchCurveICP::sliceMatrices(
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

bool MatchCurveICP::extractTransformation(
	Eigen::Matrix3Xd &source,
	Eigen::Matrix3Xd &target,
	Eigen::Affine3d &transformation)
{
	// classical ICP solution

	// center
	Eigen::Vector3d vecSCenter = source.rowwise().mean();
	Eigen::Vector3d vecTCenter = target.rowwise().mean();
	Eigen::Matrix3Xd transSource = source.colwise() - vecSCenter;
	Eigen::Matrix3Xd transTarget = target.colwise() - vecTCenter;

	// scale
	double scaScale = 1.0;
	/* // HACK: no scale during ICP
	if (true) {
		// uniform scaling
		//double scaSScale = transSource.squaredNorm();
		//double scaTScale = transTarget.squaredNorm();
		//double scaSScale = (transSource*transSource.transpose()).diagonal().maxCoeff();
		//double scaTScale = (transTarget*transTarget.transpose()).diagonal().maxCoeff();
		double scaSScale = (transSource.rowwise().maxCoeff() - transSource.rowwise().minCoeff()).maxCoeff();
		double scaTScale = (transTarget.rowwise().maxCoeff() - transTarget.rowwise().minCoeff()).maxCoeff();
		scaScale = scaTScale / scaSScale;
	}
	*/

	// rotate
	Eigen::Matrix3d matRotate;
	matRotate.setIdentity();
	if (true) {
		// horizontal rotation only
		Eigen::Matrix2Xd xzSource, xzTarget;
		xzSource.resize(2, transSource.cols());
		xzTarget.resize(2, transTarget.cols());
		xzSource << transSource.row(0), transSource.row(2);
		xzTarget << transTarget.row(0), transTarget.row(2);
		Eigen::JacobiSVD<Eigen::Matrix2d> svd(xzTarget*xzSource.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::Matrix2d matRotateXZ = svd.matrixU() * svd.matrixV().transpose();
		matRotate(0, 0) = matRotateXZ(0, 0);
		matRotate(0, 2) = matRotateXZ(0, 1);
		matRotate(2, 0) = matRotateXZ(1, 0);
		matRotate(2, 2) = matRotateXZ(1, 1);
	}

	// translate
	Eigen::Vector3d vecTranslate = vecTCenter - scaScale * matRotate * vecSCenter;

	// put everything together
	transformation.setIdentity();
	transformation.prescale(scaScale);
	transformation.prerotate(matRotate);
	transformation.pretranslate(vecTranslate);
	
	return true;
}

bool MatchCurveICP::visualize(
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