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

#include "MatchLabeledICP.h"

#include <vector>
#include <set>

#include "Segment/SegmentUtil.h"

#include "Utility/PlyExporter.h"

#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

bool MatchLabeledICP::run(
	int iteration,
	int mode,
	int primitive,
	Eigen::Matrix3Xd &sourcePoints,
	Eigen::Matrix3Xd &targetPoints,
	Eigen::VectorXi &sourceLabels,
	Eigen::VectorXi &targetLabels,
	Eigen::Vector3d &orientation,
	Eigen::Affine3d &transformation,
	Eigen::VectorXd *labelWeights)
{

	transformation.setIdentity();
	if (!sourcePoints.size() || !targetPoints.size()) return true;

	// pre-alignment
	Eigen::Matrix3Xd preSource = sourcePoints;
	Eigen::Matrix3Xd preTarget = targetPoints;
	Eigen::Affine3d preXform, postXform;
	if (!preAlign(preSource, preTarget, orientation, preXform, postXform)) return false;
	if (!fit(preSource, preTarget, sourceLabels, targetLabels, transformation)) return false; // HACK: skip...

	bool debugging = false;
	// ICP iteration
	double lastError = DBL_MAX;
	for (int iterID = 0; iterID<iteration; iterID++) {

		// find matched neighbors
		Eigen::Matrix3Xd matSP = transformation * preSource;
		Eigen::Matrix3Xd matTP = preTarget;
		Eigen::VectorXi vecSL = sourceLabels;
		Eigen::VectorXi vecTL = targetLabels;
		Eigen::VectorXd matW;
		if(debugging) if (!visualize("ICP-match-before.ply", preSource, preTarget, vecSL, vecTL, transformation)) return false;
		if (!findMatchedNeighbors(matSP, vecSL, matTP, vecTL, matW, labelWeights)) return false;
		if (matSP.size() == 0 || matTP.size() == 0) break; // no matched point pairs found

		// align matched points
		Eigen::Affine3d newTransformation = Eigen::Affine3d::Identity();
		if (debugging) if (!visualize("ICP-iteration-before.ply", matSP, matTP, vecSL, vecTL, newTransformation)) return false;
		if (debugging) if (!visualizeLink("ICP-link-before.ply", matSP, matTP, vecSL, vecTL, newTransformation)) return false;
		if (!extractTransformation(mode, primitive, matSP, matTP, matW, newTransformation)) return false;
		transformation = newTransformation * transformation;
		if (debugging) if (!visualize("ICP-iteration.ply", matSP, matTP, vecSL, vecTL, newTransformation)) return false;
		if (debugging) if (!visualizeLink("ICP-link.ply", matSP, matTP, vecSL, vecTL, newTransformation)) return false;

		// check convergence
		// NOTE: don't check it like this -- will easily lead to local minimum
		matSP = newTransformation * matSP;
		Eigen::VectorXd pointError = (matTP - matSP).colwise().squaredNorm();
		double currentError = pointError.dot(matW);
		//cout << "Current error: " << currentError << endl;
		//if (currentError > lastError*0.99) break; // UNDONE: param ICP convergence threshold
		lastError = currentError;

		if (debugging) system("pause");
	}

	if (!transformation.matrix().allFinite()) {
		transformation.setIdentity();
	}

	transformation = postXform * transformation * preXform;

	return true;
}

bool MatchLabeledICP::error(
	Eigen::Matrix3Xd &sourcePoints,
	Eigen::Matrix3Xd &targetPoints,
	Eigen::VectorXi &sourceLabels,
	Eigen::VectorXi &targetLabels,
	Eigen::VectorXd &outLabelError)
{

	// build KD trees
	vector<SKDTree*> targetTrees;
	vector<SKDTreeData> targetTreesData;
	if (!buildKDTrees(targetPoints, targetLabels, targetTrees, targetTreesData)) return false;

	// get nearest neighbors
	vector<int> tgtNNIndices;
	if (!findNearestNeighbors(targetTrees, sourcePoints, sourceLabels, tgtNNIndices)) return false;
	Eigen::Matrix3Xd slicedTargetPoints = targetPoints;
	Eigen::VectorXi slicedTargetLabels = targetLabels;
	if (!sliceMatrices(tgtNNIndices, slicedTargetPoints, slicedTargetLabels)) return false;

	if (false) {
		// use average square distance

		Eigen::VectorXd pointSquareError = (sourcePoints - slicedTargetPoints).colwise().squaredNorm();

		int numPoints = (int)sourcePoints.cols();
		int numLabels = sourceLabels.maxCoeff() + 1;
		Eigen::VectorXd labelSumSquareError(numLabels);
		Eigen::VectorXi labelPointCount(numLabels);
		labelSumSquareError.setZero();
		labelPointCount.setZero();
		for (int pointID = 0; pointID < numPoints; pointID++) {
			int labelID = sourceLabels[pointID];
			labelSumSquareError[labelID] += pointSquareError[pointID];
			labelPointCount[labelID]++;
		}
		labelPointCount = labelPointCount.cwiseMax(1);
		outLabelError = labelSumSquareError.array() / labelPointCount.cast<double>().array();
	}

	if (true) {
		// use square distance of nearest points

		vector<int> srcLabelIndices;
		if (!findNearestByLabel(sourcePoints, slicedTargetPoints, sourceLabels, srcLabelIndices)) return false;
		int numLabels = (int)srcLabelIndices.size();

		outLabelError.resize(numLabels);
		for (int labelID = 0; labelID < numLabels; labelID++) {
			int sourcePointID = srcLabelIndices[labelID];
			int targetPointID = tgtNNIndices[sourcePointID];
			outLabelError[labelID] = (sourcePoints.col(sourcePointID) - targetPoints.col(targetPointID)).squaredNorm();
		}
	}

	if (false) {
		// use square distance of farthest points

		vector<int> srcLabelIndices;
		if (!findFarthestByLabel(sourcePoints, slicedTargetPoints, sourceLabels, srcLabelIndices)) return false;
		int numLabels = (int)srcLabelIndices.size();

		outLabelError.resize(numLabels);
		for (int labelID = 0; labelID < numLabels; labelID++) {
			int sourcePointID = srcLabelIndices[labelID];
			int targetPointID = tgtNNIndices[sourcePointID];
			outLabelError[labelID] = (sourcePoints.col(sourcePointID) - targetPoints.col(targetPointID)).squaredNorm();
		}
	}

	return true;
}

bool MatchLabeledICP::preAlign(
	Eigen::Matrix3Xd &source,
	Eigen::Matrix3Xd &target,
	Eigen::Vector3d &orientation,
	Eigen::Affine3d &preXform,
	Eigen::Affine3d &postXform)
{

	// regularize major axis
	Eigen::Vector3d majorAxis = orientation;
	if (true) {
		int majorDim;
		majorAxis.cwiseAbs().maxCoeff(&majorDim);
		if (majorAxis[majorDim] < 0) majorAxis = -majorAxis;
		double angle = cml::deg(cml::acos_safe(majorAxis[majorDim]));
		if (angle < StyleSynthesisConfig::mAssemble_OrientationAngleThreshold) {
			majorAxis = Eigen::Vector3d::Unit(majorDim);
		}
	}

	// rotate to align major axis to "up"
	Eigen::AngleAxisd rotation = Eigen::AngleAxisd::Identity();
	Eigen::Vector3d normal = majorAxis.cross(Eigen::Vector3d::UnitY());
	if (normal.squaredNorm()) {
		Eigen::Vector3d axisX = normal.normalized();
		Eigen::Vector3d axisY = majorAxis;
		Eigen::Vector3d axisZ = axisX.cross(axisY).normalized();
		Eigen::Matrix3d matCS;
		matCS << axisX, axisY, axisZ;
		rotation = matCS.transpose();
	}

	// apply pre-alignment to points
	source = rotation.matrix() * source;
	target = rotation.matrix() * target;

	// output transformations
	preXform.setIdentity();
	preXform.prerotate(rotation);
	postXform.setIdentity();
	postXform.prerotate(rotation.inverse());

	return true;
}

bool MatchLabeledICP::fit(
	Eigen::Matrix3Xd &sourcePoints,
	Eigen::Matrix3Xd &targetPoints,
	Eigen::VectorXi &sourceLabels,
	Eigen::VectorXi &targetLabels,
	Eigen::Affine3d &transformation)
{
	// align center (label weighted)

	int numSourcePoints = (int)sourcePoints.cols();
	int numTargetPoints = (int)targetPoints.cols();
	int numSourceLabels = sourceLabels.maxCoeff() + 1;
	int numTargetLabels = targetLabels.maxCoeff() + 1;

	Eigen::Matrix3Xd srcLabelCenter(3, numSourceLabels);
	Eigen::Matrix3Xd tgtLabelCenter(3, numTargetLabels);
	Eigen::VectorXi srcLabelPointCount(numSourceLabels);
	Eigen::VectorXi tgtLabelPointCount(numTargetLabels);
	srcLabelCenter.setZero();
	tgtLabelCenter.setZero();
	srcLabelPointCount.setZero();
	tgtLabelPointCount.setZero();

	for (int pointID = 0; pointID < numSourcePoints; pointID++) {
		int labelID = sourceLabels[pointID];
		srcLabelCenter.col(labelID) += sourcePoints.col(pointID);
		srcLabelPointCount[labelID]++;
	}
	for (int pointID = 0; pointID < numTargetPoints; pointID++) {
		int labelID = targetLabels[pointID];
		tgtLabelCenter.col(labelID) += targetPoints.col(pointID);
		tgtLabelPointCount[labelID]++;
	}

	Eigen::Vector3d srcCenter = Eigen::Vector3d::Zero();
	Eigen::Vector3d tgtCenter = Eigen::Vector3d::Zero();
	int numSourceValidLabels = 0;
	int numTargetValidLabels = 0;
	for (int labelID = 0; labelID < numSourceLabels; labelID++) {
		if (srcLabelPointCount[labelID]) {
			srcCenter += srcLabelCenter.col(labelID) / srcLabelPointCount[labelID];
			numSourceValidLabels++;
		}
	}
	for (int labelID = 0; labelID < numTargetLabels; labelID++) {
		if (tgtLabelPointCount[labelID]) {
			tgtCenter += tgtLabelCenter.col(labelID) / tgtLabelPointCount[labelID];
			numTargetValidLabels++;
		}
	}
	if (numSourceValidLabels) srcCenter *= 1.0 / numSourceValidLabels;
	if (numTargetValidLabels) tgtCenter *= 1.0 / numTargetValidLabels;

	transformation.setIdentity();
	transformation.pretranslate(tgtCenter - srcCenter);

	return true;
}

bool MatchLabeledICP::buildKDTrees(
	Eigen::Matrix3Xd &points,
	Eigen::VectorXi &labels,
	vector<SKDTree*> &trees,
	vector<SKDTreeData> &treesData)
{	
	int numLabels = labels.maxCoeff() + 1;
	treesData.resize(numLabels);
	for (int labelID = 0; labelID < numLabels; labelID++) {
		treesData[labelID].clear();
	}
	
	int numPoints = (int)points.cols();
	for (int pointID = 0; pointID < numPoints; pointID++) {
		int labelID = labels[pointID];
		SKDT::NamedPoint point((float)points(0, pointID), (float)points(1, pointID), (float)points(2, pointID), (size_t)pointID);
		treesData[labelID].push_back(point);
	}
	
	trees.resize(numLabels);
	for (int labelID = 0; labelID < numLabels; labelID++) {
		trees[labelID] = new SKDTree(treesData[labelID].begin(), treesData[labelID].end());
	}
	
	return true;
}

bool MatchLabeledICP::findMatchedNeighbors(
	Eigen::Matrix3Xd &inoutSourcePoints,
	Eigen::VectorXi &inoutSourceLabels,
	Eigen::Matrix3Xd &inoutTargetPoints,
	Eigen::VectorXi &inoutTargetLabels,
	Eigen::VectorXd &outWeights,
	Eigen::VectorXd *labelWeights)
{

	// build KD trees
	vector<SKDTree*> targetTrees;
	vector<SKDTreeData> targetTreesData;
	if (!buildKDTrees(inoutTargetPoints, inoutTargetLabels, targetTrees, targetTreesData)) return false;

	// get nearest neighbors
	vector<int> tgtNNIndices;
	if (!findNearestNeighbors(targetTrees, inoutSourcePoints, inoutSourceLabels, tgtNNIndices)) return false;
	Eigen::Matrix3Xd slicedTargetPoints = inoutTargetPoints;
	Eigen::VectorXi slicedTargetLabels = inoutTargetLabels;
	if (!sliceMatrices(tgtNNIndices, slicedTargetPoints, slicedTargetLabels)) return false;

	if (true) {

		// use weighted ICP scheme (each label sums up to same weights)

		int numLabels = inoutSourceLabels.maxCoeff() + 1;

		// prune outlier points

		vector<int> prunedIndices;
		Eigen::Matrix3Xd prunedSourcePoints = inoutSourcePoints;
		Eigen::Matrix3Xd prunedTargetPoints = slicedTargetPoints;
		Eigen::VectorXi prunedSourceLabels = inoutSourceLabels;
		Eigen::VectorXi prunedTargetLabels = slicedTargetLabels;

		if (StyleSynthesisConfig::mAssemble_SlotsMatchingAggressiveFiltering) {

			// keep only closest for each target point

			Eigen::RowVectorXd pointDists = (inoutSourcePoints - slicedTargetPoints).colwise().norm();

			int numTargetNNPoints;
			for (int nnid : tgtNNIndices) numTargetNNPoints = max(numTargetNNPoints, nnid);
			numTargetNNPoints++;

			vector<double> targetNNPointDistance(numTargetNNPoints, DBL_MAX);
			for (int id = 0; id < (int)tgtNNIndices.size(); id++) {
				int nnid = tgtNNIndices[id];
				double dist = pointDists[id];
				if (dist < targetNNPointDistance[nnid]) {
					targetNNPointDistance[nnid] = dist;
				}
			}

			prunedIndices.clear();
			for (int id = 0; id < (int)tgtNNIndices.size(); id++) {
				int nnid = tgtNNIndices[id];
				double dist = pointDists[id];
				if (fabs(targetNNPointDistance[nnid] - dist) < 1e-6) {
					prunedIndices.push_back(id);
				}
			}
			//cout << "first pruning: " << prunedSourcePoints.cols() << " => " << prunedIndices.size() << endl;

			if (!sliceMatrices(prunedIndices, prunedSourcePoints, prunedSourceLabels)) return false;
			if (!sliceMatrices(prunedIndices, prunedTargetPoints, prunedTargetLabels)) return false;
		}

		if (true) {

			// prune by distance percentile

			// compute pruning distance
			double multiplier = StyleSynthesisConfig::mAssemble_SlotPointsPruningMultiplier;
			double percentile = StyleSynthesisConfig::mAssemble_SlotPointsPruningPercentile;
			vector<double> labelPruningDist(numLabels);
			vector<vector<double>> labelDists(numLabels);
			Eigen::RowVectorXd pointDists = (prunedSourcePoints - prunedTargetPoints).colwise().norm();
			for (int pointID = 0; pointID < (int)prunedSourceLabels.size(); pointID++) {
				int label = prunedSourceLabels[pointID];
				labelDists[label].push_back(pointDists[pointID]);
			}
			for (int labelID = 0; labelID < numLabels; labelID++) {
				vector<double> dists = labelDists[labelID];
				if ((int)dists.size() >= 1) {
					int n = cml::clamp((int)(dists.size() * percentile), 0, (int)dists.size() - 1);
					nth_element(dists.begin(), dists.begin() + n, dists.end());
					labelPruningDist[labelID] = dists[n] * multiplier;
				} else {
					labelPruningDist[labelID] = 0; // no points for this label
				}
			}

			// slice inliner points
			prunedIndices.clear();
			for (int pointID = 0; pointID < (int)prunedSourceLabels.size(); pointID++) {
				int label = prunedSourceLabels[pointID];
				if (pointDists[pointID] <= labelPruningDist[label]) {
					prunedIndices.push_back(pointID);
				}
			}
			//cout << "second pruning: " << prunedSourcePoints.cols() << " => " << prunedIndices.size() << endl;

			if (!sliceMatrices(prunedIndices, prunedSourcePoints, prunedSourceLabels)) return false;
			if (!sliceMatrices(prunedIndices, prunedTargetPoints, prunedTargetLabels)) return false;
		}

		// count number of points for each label
		
		Eigen::VectorXi pointCount(numLabels);
		pointCount.setZero();
		for (int id = 0; id < (int)prunedSourceLabels.size(); id++) {
			int labelID = prunedSourceLabels[id];
			pointCount[labelID] += 1;
		}
		Eigen::VectorXd unnormalizedLabelWeights = pointCount.cast<double>().array().max(1.0).inverse();
		if (labelWeights) {
			if ((int)labelWeights->size() != numLabels) {
				cout << "Error: incorrect number of labels" << endl;
				return false;
			}
			unnormalizedLabelWeights = unnormalizedLabelWeights.array() * labelWeights->array();
		}

		// compute weights
		Eigen::VectorXd weights((int)prunedSourceLabels.size());
		if (true) {
			weights.setOnes();

			// compute point weights (unestablished approach)
			//Eigen::VectorXd pointDist = prunedSourcePoints.colwise().norm();
			//double maxDist = pointDist.maxCoeff();
			//double minDist = pointDist.minCoeff();
			//Eigen::VectorXd pointWeights = ((pointDist.array() - maxDist) / (maxDist - minDist)).exp();
			//weights = pointWeights;

			for (int id = 0; id < (int)prunedSourceLabels.size(); id++) {
				int labelID = prunedSourceLabels[id];
				weights[id] *= unnormalizedLabelWeights[labelID];
			}
		}
		weights /= weights.sum(); // normalization

		// output stuffs
		inoutSourcePoints.swap(prunedSourcePoints);
		inoutSourceLabels.swap(prunedSourceLabels);
		inoutTargetPoints.swap(prunedTargetPoints);
		inoutTargetLabels.swap(prunedTargetLabels);
		outWeights.swap(weights);
	}

	if (false) {

		// classical ICP on closest points for each label

		vector<int> labelIndices;
		if (!findNearestByLabel(inoutSourcePoints, slicedTargetPoints, inoutSourceLabels, labelIndices)) return false;
		if (!sliceMatrices(labelIndices, inoutSourcePoints, inoutSourceLabels)) return false;
		if (!sliceMatrices(labelIndices, inoutTargetPoints, inoutTargetLabels)) return false;
		outWeights = Eigen::VectorXd::Ones((int)labelIndices.size());
		outWeights /= outWeights.sum();
	}

	// don't forget to delete things
	for (auto tree : targetTrees) delete tree;
	
	return true;
}

bool MatchLabeledICP::findNearestNeighbors(
	vector<SKDTree*> &trees,
	Eigen::Matrix3Xd &inPoints,
	Eigen::VectorXi &inLabels,
	vector<int> &outIndices)
{
	outIndices.resize(inPoints.cols());
#pragma omp parallel for
	for (int j = 0; j < inPoints.cols(); j++) {
		SKDT::NamedPoint queryPoint((float)inPoints(0, j), (float)inPoints(1, j), (float)inPoints(2, j));
		Thea::BoundedSortedArray<SKDTree::Neighbor> queryResult(1);
		int labelID = inLabels[j];
		trees[labelID]->kClosestElements<Thea::MetricL2>(queryPoint, queryResult);
		if (queryResult.isEmpty()) outIndices[j] = 0; // too far away, whatever...
		else outIndices[j] = (int)trees[labelID]->getElements()[queryResult[0].getIndex()].id;
	}

	return true;
}

bool MatchLabeledICP::sliceMatrices(
	vector<int> &inIndices,
	Eigen::Matrix3Xd &inoutPoints,
	Eigen::VectorXi &inoutLabels)
{	
	Eigen::Matrix3Xd tmpPoints(inoutPoints.rows(), inIndices.size());
	Eigen::VectorXi tmpLabels(inIndices.size());
	for (int j = 0; j < (int)inIndices.size(); j++) {
		tmpPoints.col(j) = inoutPoints.col(inIndices[j]);
		tmpLabels[j] = inoutLabels[inIndices[j]];
	}
	inoutPoints.swap(tmpPoints);
	inoutLabels.swap(tmpLabels);
	
	return true;
}

static bool sliceLabels(
	Eigen::VectorXi &inLabels,
	vector<int> &inIndices,
	Eigen::VectorXi &outLabels)
{

}

bool MatchLabeledICP::findNearestByLabel(
	Eigen::Matrix3Xd &inPoints,
	Eigen::Matrix3Xd &slicedPoints,
	Eigen::VectorXi &inLabels,
	vector<int> &outIndices)
{
	int numLabels = (int)inLabels.maxCoeff() + 1;
	vector<int> labelMinIndices(numLabels, -1);
	vector<double> labelMinDistances(numLabels, DBL_MAX);

	Eigen::RowVectorXd vecDist = (inPoints - slicedPoints).colwise().norm();
	for (int pointID = 0; pointID < vecDist.size(); pointID++) {
		int labelID = inLabels[pointID];
		double dist = vecDist[pointID];
		if (dist < labelMinDistances[labelID]) {
			labelMinDistances[labelID] = dist;
			labelMinIndices[labelID] = pointID;
		}
	}

	outIndices = labelMinIndices;

	return true;
}

bool MatchLabeledICP::findFarthestByLabel(
	Eigen::Matrix3Xd &inPoints,
	Eigen::Matrix3Xd &slicedPoints,
	Eigen::VectorXi &inLabels,
	vector<int> &outIndices)
{
	int numLabels = (int)inLabels.maxCoeff() + 1;
	vector<int> labelMaxIndices(numLabels, -1);
	vector<double> labelMaxDistances(numLabels, 0);

	Eigen::RowVectorXd vecDist = (inPoints - slicedPoints).colwise().norm();
	for (int pointID = 0; pointID < vecDist.size(); pointID++) {
		int labelID = inLabels[pointID];
		double dist = vecDist[pointID];
		if (dist > labelMaxDistances[labelID]) {
			labelMaxDistances[labelID] = dist;
			labelMaxIndices[labelID] = pointID;
		}
	}

	outIndices = labelMaxIndices;

	return true;
}

bool MatchLabeledICP::extractTransformation(
	int mode,
	int primitive,
	Eigen::Matrix3Xd &source,
	Eigen::Matrix3Xd &target,
	Eigen::VectorXd &weights,
	Eigen::Affine3d &transformation)
{
	// weighted ICP

	Eigen::Vector3d srcCenter = source * weights;
	Eigen::Vector3d tgtCenter = target * weights;
	Eigen::Matrix3Xd matS = source.colwise() - srcCenter;
	Eigen::Matrix3Xd matT = target.colwise() - tgtCenter;

	Eigen::Vector3d scale = Eigen::Vector3d::Ones();
	Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();

	Eigen::Vector3d vecSScale = (matS * weights.asDiagonal() * matS.transpose()).diagonal().cwiseSqrt();
	Eigen::Vector3d vecTScale = (matT * weights.asDiagonal() * matT.transpose()).diagonal().cwiseSqrt();

	if (primitive == 0) {

		/////////////// stick ///////////////

		if (mode >= 1 && mode != 4 && mode != 6) {
			// scaling along major axis
			if (vecSScale[1]) scale[1] = vecTScale[1] / vecSScale[1];
		}
		if (mode >= 2 && mode != 4 && mode != 6) {
			// uniform scaling on cross section
			double scaleT = sqrt(cml::sqr(vecTScale[0]) + cml::sqr(vecTScale[2]));
			double scaleS = sqrt(cml::sqr(vecSScale[0]) + cml::sqr(vecSScale[2]));
			if(scaleS) scale[0] = scale[2] = scaleT / scaleS;
		}
		if (mode == 4 || mode == 6) {
			// uniform scaling
			double scaleT = vecTScale.norm();
			double scaleS = vecSScale.norm();
			if (scaleS) scale.setConstant(scaleT / scaleS);
		}
		for (int dim = 0; dim < 3; dim++) { // prevent weird scaling
			if (fabs(scale[dim]) > 5.0 || fabs(scale[dim]) < 0.2) scale[dim] = 1.0;
		}
		matS = scale.asDiagonal() * matS;

		if (mode == 3) {
			// rotation along major axis
			if (!extractPlanarRotation(matS, matT, weights, rotation)) return false;
		}
		if (mode >= 4 && mode != 6) {
			// free rotation
			if (!extractFreeRotation(matS, matT, weights, rotation)) return false;
		}

	} else if (primitive == 1) {

		/////////////// plane ///////////////

		if (mode == 1 || mode == 4 || mode > 6) {
			// scaling along X axis
			if (vecSScale[0]) scale[0] = vecTScale[0] / vecSScale[0];
		}
		if (mode == 2 || mode == 4 || mode > 6) {
			// scaling along Z axis
			if (vecSScale[2]) scale[2] = vecTScale[2] / vecSScale[2];
		}
		if (mode == 3) {
			// uniform scaling on cross section
			double scaleT = sqrt(cml::sqr(vecTScale[0]) + cml::sqr(vecTScale[2]));
			double scaleS = sqrt(cml::sqr(vecSScale[0]) + cml::sqr(vecSScale[2]));
			if (scaleS) scale[0] = scale[2] = scaleT / scaleS;
		}
		if (mode >= 7) {
			// scaling along Y axis
			if (vecSScale[1]) scale[1] = vecTScale[1] / vecSScale[1];
		}
		if (mode == 5 || mode == 6) {
			// uniform scaling
			double scaleT = vecTScale.norm();
			double scaleS = vecSScale.norm();
			if (scaleS) scale.setConstant(scaleT / scaleS);
		}
		for (int dim = 0; dim < 3; dim++) { // prevent weird scaling
			if (fabs(scale[dim]) > 5.0 || fabs(scale[dim]) < 0.2) scale[dim] = 1.0;
		}
		matS = scale.asDiagonal() * matS;

		if (mode == 5) {
			// rotation along major axis
			if (!extractPlanarRotation(matS, matT, weights, rotation)) return false;
		}
		if (mode >= 6) {
			// free rotation
			if (!extractFreeRotation(matS, matT, weights, rotation)) return false;
		}

	} else if (primitive == 2) {

		/////////////// sphere ///////////////

		if (mode == 1 || mode == 3) {
			// uniform scaling
			double scaleT = vecTScale.norm();
			double scaleS = vecSScale.norm();
			if (scaleS) scale.setConstant(scaleT / scaleS);
		}
		if (mode == 2 || mode >= 4) {
			// non-uniform scaling
			for (int dim = 0; dim < 3; dim++) {
				if (vecSScale[dim]) scale[dim] = vecTScale[dim] / vecSScale[dim];
			}
		}
		for (int dim = 0; dim < 3; dim++) { // prevent weird scaling
			if (fabs(scale[dim]) > 5.0 || fabs(scale[dim]) < 0.2) scale[dim] = 1.0;
		}
		matS = scale.asDiagonal() * matS;

		if (mode >= 3) {
			// free rotation
			if (!extractFreeRotation(matS, matT, weights, rotation)) return false;
		}
	} else {
		cout << "Error: unrecognized primitive " << primitive << endl;
		return false;
	}

	Eigen::Vector3d translation = tgtCenter - rotation * scale.asDiagonal() * srcCenter;
	//cout << "Xform: " << mode << ", " << primitive << endl;
	//cout << "scale: " << scale.transpose() << endl;
	//cout << "translate: " << translation.transpose() << endl;
	//cout << "rotation:" << endl;
	//cout << rotation << endl;

	transformation.setIdentity();
	transformation.prescale(scale);
	transformation.prerotate(rotation);
	transformation.pretranslate(translation);
	
	return true;
}


bool MatchLabeledICP::extractPlanarRotation(
	Eigen::Matrix3Xd &source,
	Eigen::Matrix3Xd &target,
	Eigen::VectorXd &weights,
	Eigen::Matrix3d &rotation)
{

	Eigen::Matrix2Xd projSource(2, source.cols());
	Eigen::Matrix2Xd projTarget(2, target.cols());
	projSource << source.row(0), source.row(2);
	projTarget << target.row(0), target.row(2);

	Eigen::JacobiSVD<Eigen::Matrix2d> svd(projTarget*weights.asDiagonal()*projSource.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix2d matU = svd.matrixU();
	Eigen::Matrix2d matV = svd.matrixV();
	Eigen::Matrix2d rotation2D = matU * matV.transpose();
	if (rotation2D.determinant() < 0) {
		rotation2D = matU * Eigen::Vector2d(1, -1).asDiagonal() * matV.transpose();
	}
	rotation.setIdentity();
	rotation(0, 0) = rotation2D(0, 0);
	rotation(0, 2) = rotation2D(0, 1);
	rotation(2, 0) = rotation2D(1, 0);
	rotation(2, 2) = rotation2D(1, 1);

	return true;
}

bool MatchLabeledICP::extractFreeRotation(
	Eigen::Matrix3Xd &source,
	Eigen::Matrix3Xd &target,
	Eigen::VectorXd &weights,
	Eigen::Matrix3d &rotation)
{
	Eigen::JacobiSVD<Eigen::Matrix3d> svd(target*weights.asDiagonal()*source.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix3d matU = svd.matrixU();
	Eigen::Matrix3d matV = svd.matrixV();
	rotation = matU * matV.transpose();
	if (rotation.determinant() < 0) { // HACK: no reflection
		rotation = matU * Eigen::Vector3d(1, 1, -1).asDiagonal() * matV.transpose();
	}

	return true;
}

bool MatchLabeledICP::visualize(
	string filename,
	Eigen::Matrix3Xd &sourcePoints,
	Eigen::Matrix3Xd &targetPoints,
	Eigen::VectorXi &sourceLabels,
	Eigen::VectorXi &targetLabels,
	Eigen::Affine3d &transformation)
{

	int numSourcePoints = (int)sourcePoints.cols();
	int numTargetPoints = (int)targetPoints.cols();

	Eigen::Matrix3Xd transformedPoints = transformation * sourcePoints;
	
	vector<vec3> srcVec(numSourcePoints);
	vector<vec3> tgtVec(numTargetPoints);
	vector<vec3i> srcColor(numSourcePoints);
	vector<vec3i> tgtColor(numTargetPoints);

	for (int pointID = 0; pointID < transformedPoints.cols(); pointID++) {
		for (int dim = 0; dim < 3; dim++) {
			srcVec[pointID][dim] = (float)transformedPoints(dim, pointID);
		}
		srcColor[pointID] = SegmentUtil::colorMapping(sourceLabels[pointID]) / 2; // use darker color for source points
	}
	for (int pointID = 0; pointID < targetPoints.cols(); pointID++) {
		for (int dim = 0; dim < 3; dim++) {
			tgtVec[pointID][dim] = (float)targetPoints(dim, pointID);
		}
		tgtColor[pointID] = SegmentUtil::colorMapping(targetLabels[pointID]);
	}

	PlyExporter pe;
	if (!pe.addPoint(&tgtVec, 0, &tgtColor)) return false;
	if (!pe.addPoint(&srcVec, 0, &srcColor)) return false;
	if (!pe.output(filename)) return false;
	
	return true;
}

bool MatchLabeledICP::visualizeLink(
	string filename,
	Eigen::Matrix3Xd &sourcePoints,
	Eigen::Matrix3Xd &targetPoints,
	Eigen::VectorXi &sourceLabels,
	Eigen::VectorXi &targetLabels,
	Eigen::Affine3d &transformation)
{

	int numSourcePoints = (int)sourcePoints.cols();
	int numTargetPoints = (int)targetPoints.cols();
	if (numSourcePoints != numTargetPoints) {
		cout << "Error: incorrect source/target point amount" << endl;
		return false;
	}

	Eigen::Matrix3Xd transformedPoints = transformation * sourcePoints;

	vector<vec3> lines;

	for (int pointID = 0; pointID < numSourcePoints; pointID++) {
		vec3 srcPoint, tgtPoint;
		for (int dim = 0; dim < 3; dim++) {
			srcPoint[dim] = (float)transformedPoints(dim, pointID);
			tgtPoint[dim] = (float)targetPoints(dim, pointID);
		}
		lines.push_back(srcPoint);
		lines.push_back(tgtPoint);

		vec3 normal = tgtPoint - srcPoint;
		vec3 axis1 = cml::normalize(cml::cross(fabs(normal[0]) > 0.5f ? vec3(0.0f, 1.0f, 0.0f) : vec3(1.0f, 0.0f, 0.0f), normal));
		vec3 axis2 = cml::normalize(cml::cross(normal, axis1));
		float len = normal.length() * 0.05f;
		lines.push_back(srcPoint - axis1*len);
		lines.push_back(srcPoint + axis1*len);
		lines.push_back(srcPoint - axis2*len);
		lines.push_back(srcPoint + axis2*len);
	}

	PlyExporter pe;
	if (!pe.addLine(&lines)) return false;
	if (!pe.output(filename)) return false;

	return true;
}

bool MatchLabeledICP::buildMatrices(
	vector<Eigen::Vector3d> &inPoints,
	vector<int> &inLabels,
	Eigen::Matrix3Xd &outMat,
	Eigen::VectorXi &outVec)
{

	outMat.resize(3, inPoints.size());
	for (int id = 0; id < (int)inPoints.size(); id++) {
		outMat.col(id) = inPoints[id];
	}
	outVec = Eigen::Map<Eigen::VectorXi>(inLabels.data(), inLabels.size());

	return true;
}