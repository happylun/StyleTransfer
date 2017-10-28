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

#include "MatchPrimitiveICP.h"

#include <vector>
#include <set>

#include "Utility/PlyExporter.h"

#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

bool MatchPrimitiveICP::run(
	int iteration,
	int mode,
	int primitive,
	Eigen::Matrix3Xd &source,
	Eigen::Matrix3Xd &target,
	Eigen::Vector3d &orientation,
	Eigen::Affine3d &transformation)
{

	// pre-alignment
	Eigen::Matrix3Xd preSource = source;
	Eigen::Matrix3Xd preTarget = target;
	Eigen::Affine3d preXform, postXform;
	if (!preAlign(preSource, preTarget, orientation, preXform, postXform)) return false;

	// initialization
	SKDTree tree;
	SKDTreeData treeData;
	if (!buildKDTree(preTarget, tree, treeData)) return false;

	Eigen::Affine3d initialTransformation = Eigen::Affine3d::Identity();
	if (!extractFittingBoundingBox(mode, primitive, preSource, preTarget, initialTransformation)) return false;
	//if (!extractFittingMass(mode, primitive, preSource, preTarget, initialTransformation)) return false;
	transformation = initialTransformation;
	//if (!visualize("ICP-init.ply", preSource, preTarget, initialTransformation)) return false;

	// ICP iteration
	double lastError = DBL_MAX;
	for (int iterID = 0; iterID<iteration; iterID++) {

		// find nearest neighbors
		Eigen::Matrix3Xd matS = transformation * preSource;
		Eigen::Matrix3Xd matT;
		vector<int> slices;
		if (!findNearestNeighbors(tree, matS, slices)) return false;
		if (!sliceMatrices(preTarget, slices, matT)) return false;
		if (!findMatchedNeighbors(matS, matT, slices)) return false;
		if (slices.empty()) break; // no matched points
		if (!sliceMatrices(matS, slices, matS)) return false;
		if (!sliceMatrices(matT, slices, matT)) return false;

		// align matched points
		Eigen::Affine3d newTransformation;
		if (!extractTransformation(mode, primitive, matS, matT, newTransformation)) return false;
		//if (!visualize("ICP-xform.ply", matS, matT, newTransformation)) return false;
		transformation = newTransformation * transformation;
		//if (!visualize("ICP-iteration.ply", preSource, preTarget, transformation)) return false;

		// check convergence
		matS = transformation * matS;
		double currentError = (matT - matS).squaredNorm() / matS.cols();
		if (currentError > lastError*0.99) break; // UNDONE: param ICP convergence threshold
		lastError = currentError;
	}

	if (!transformation.matrix().allFinite()) {
		transformation = initialTransformation;
	}

	transformation = postXform * transformation * preXform;

	return true;
}

bool MatchPrimitiveICP::fit(
	int mode,
	int primitive,
	Eigen::Matrix3Xd &source,
	Eigen::Matrix3Xd &target,
	Eigen::Vector3d &orientation,
	Eigen::Affine3d &transformation)
{

	// pre-alignment
	Eigen::Matrix3Xd preSource = source;
	Eigen::Matrix3Xd preTarget = target;
	Eigen::Affine3d preXform, postXform;
	if (!preAlign(preSource, preTarget, orientation, preXform, postXform)) return false;

	// get fitting transformation
	if (!extractFittingBoundingBox(mode, primitive, preSource, preTarget, transformation)) return false;
	//if (!extractFittingMass(mode, primitive, preSource, preTarget, transformation)) return false;

	transformation = postXform * transformation * preXform;

	return true;
}

bool MatchPrimitiveICP::error(
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

bool MatchPrimitiveICP::preAlign(
	Eigen::Matrix3Xd &source,
	Eigen::Matrix3Xd &target,
	Eigen::Vector3d &orientation,
	Eigen::Affine3d &preXform,
	Eigen::Affine3d &postXform)
{

	// move center to origin
	Eigen::Vector3d srcCenter = source.rowwise().mean();
	Eigen::Vector3d tgtCenter = target.rowwise().mean();

	// regularize major axis
	Eigen::Vector3d majorAxis = orientation;
	if (!regularizeAxis(majorAxis)) return false;

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
	source.colwise() -= srcCenter;
	target.colwise() -= tgtCenter;
	source = rotation.matrix() * source;
	target = rotation.matrix() * target;

	// output transformations
	preXform.setIdentity();
	preXform.pretranslate(-srcCenter);
	preXform.prerotate(rotation);
	postXform.setIdentity();
	postXform.prerotate(rotation.inverse());
	postXform.pretranslate(tgtCenter);

	return true;
}

bool MatchPrimitiveICP::regularizeAxis(
	Eigen::Vector3d &axis)
{
	int majorDim;
	axis.cwiseAbs().maxCoeff(&majorDim);
	double angle = cml::deg(cml::acos_safe(fabs(axis[majorDim])));
	if (angle < StyleSynthesisConfig::mAssemble_OrientationAngleThreshold) {
		axis = Eigen::Vector3d::Unit(majorDim) * (axis[majorDim] > 0 ? 1 : -1);
	}

	if (StyleSynthesisConfig::mContext_HandleRotationalSymmetry && fabs(axis[1]) < 1.0) {
		// make axis parallel to ground plane
		axis[1] = 0;
		axis.normalize();
	}

	return true;
}

bool MatchPrimitiveICP::buildKDTree(
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

bool MatchPrimitiveICP::findNearestNeighbors(
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

bool MatchPrimitiveICP::findMatchedNeighbors(
	Eigen::Matrix3Xd &inSource,
	Eigen::Matrix3Xd &inTarget,
	vector<int> &outIndices)
{
	outIndices.clear();
	if (inSource.cols() == 0) return true; // empty

	const double rejectDistanceThreshold = 2.0;

	Eigen::ArrayXd vecD = (inSource - inTarget).colwise().norm().array();
	double maxDist;
	if (true) {
		// use r times a percentile length as clamping distance (r = 2, a = 90%?)
		vector<double> vDist(vecD.data(), vecD.data() + vecD.size());
		int n = (int)(vecD.size() * 0.1);
		nth_element(vDist.begin(), vDist.begin() + n, vDist.end());
		maxDist = vDist[n] * rejectDistanceThreshold;
	}

	auto filter = vecD < maxDist;
	outIndices.clear();
	outIndices.reserve((int)filter.count());
	for (int j = 0; j < filter.size(); j++) { // don't parallelize
		if (filter(j)) outIndices.push_back(j);
	}

	return true;
}

bool MatchPrimitiveICP::sliceMatrices(
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

bool MatchPrimitiveICP::extractTransformation(
	int mode,
	int primitive,
	Eigen::Matrix3Xd &source,
	Eigen::Matrix3Xd &target,
	Eigen::Affine3d &transformation)
{
	// classic ICP

	Eigen::AlignedBox3d srcBB(source.rowwise().minCoeff(), source.rowwise().maxCoeff());
	Eigen::AlignedBox3d tgtBB(target.rowwise().minCoeff(), target.rowwise().maxCoeff());
	Eigen::Vector3d srcBBSize = srcBB.sizes();
	Eigen::Vector3d tgtBBSize = tgtBB.sizes();
	Eigen::Vector3d srcCenter = srcBB.center();
	Eigen::Vector3d tgtCenter = tgtBB.center();

	Eigen::Matrix3Xd matS = source.colwise() - srcCenter;
	Eigen::Matrix3Xd matT = target.colwise() - tgtCenter;
	Eigen::Matrix3Xd matMS = source.colwise() - source.rowwise().mean();
	Eigen::Matrix3Xd matMT = target.colwise() - target.rowwise().mean();

	Eigen::Vector3d scale = Eigen::Vector3d::Ones();
	Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();

	bool useMassScale = false; // use mass scale or BB scale

	if (primitive == 0) {

		/////////////// stick ///////////////

		if (mode >= 1 && mode != 4 && mode != 6) {
			// scaling along major axis
			double scaleT, scaleS;
			if (useMassScale) {
				scaleT = matMT.row(1).norm();
				scaleS = matMS.row(1).norm();
			} else {
				scaleT = tgtBBSize[1];
				scaleS = srcBBSize[1];
			}
			if(scaleS) scale[1] = scaleT / scaleS;
		}
		if (mode >= 2 && mode != 4 && mode != 6) {
			// uniform scaling on cross section
			double scaleT = sqrt(matMT.row(0).squaredNorm() + matMT.row(2).squaredNorm());
			double scaleS = sqrt(matMS.row(0).squaredNorm() + matMS.row(2).squaredNorm());
			//double scaleT = sqrt(cml::sqr(tgtBBSize[0]) + cml::sqr(tgtBBSize[2]));
			//double scaleS = sqrt(cml::sqr(srcBBSize[0]) + cml::sqr(srcBBSize[2]));
			if(scaleS) scale[0] = scale[2] = scaleT / scaleS;
		}
		if (mode == 4 || mode == 6) {
			// uniform scaling
			double scaleT = matMT.norm();
			double scaleS = matMS.norm();
			if (scaleS) scale.setConstant(scaleT / scaleS);
		}
		for (int dim = 0; dim < 3; dim++) { // prevent weird scaling
			if (fabs(scale[dim]) > 5.0) scale[dim] = 1.0;
		}
		matS = scale.asDiagonal() * matS;

		if (mode == 3) {
			// rotation along major axis
			if (!extractPlanarRotation(matS, matT, rotation)) return false;
		}
		if (mode >= 4 && mode != 6) {
			// free rotation
			if (!extractFreeRotation(matS, matT, rotation)) return false;
		}

	} else if (primitive == 1) {

		/////////////// plane ///////////////

		if (mode == 1 || mode == 4 || mode > 6) {
			// scaling along X axis
			double scaleT, scaleS;
			if (useMassScale) {
				scaleT = matMT.row(0).norm();
				scaleS = matMS.row(0).norm();
			} else {
				scaleT = tgtBBSize[0];
				scaleS = srcBBSize[0];
			}
			if (scaleS) scale[0] = scaleT / scaleS;
		}
		if (mode == 2 || mode == 4 || mode > 6) {
			// scaling along Z axis
			double scaleT, scaleS;
			if (useMassScale) {
				scaleT = matMT.row(2).norm();
				scaleS = matMS.row(2).norm();
			} else {
				scaleT = tgtBBSize[2];
				scaleS = srcBBSize[2];
			}
			if (scaleS) scale[2] = scaleT / scaleS;
		}
		if (mode == 3) {
			// uniform scaling on cross section
			double scaleT = sqrt(matMT.row(0).squaredNorm() + matMT.row(2).squaredNorm());
			double scaleS = sqrt(matMS.row(0).squaredNorm() + matMS.row(2).squaredNorm());
			if(scaleS) scale[0] = scale[2] = scaleT / scaleS;
		}
		if (mode >= 7) {
			// scaling along Y axis
			double scaleT, scaleS;
			if (useMassScale) {
				scaleT = matMT.row(1).norm();
				scaleS = matMS.row(1).norm();
			} else {
				scaleT = tgtBBSize[1];
				scaleS = srcBBSize[1];
			}
			if (scaleS) scale[1] = scaleT / scaleS;
		}
		if (mode == 5 || mode == 6) {
			// uniform scaling
			double scaleT = matMT.norm();
			double scaleS = matMS.norm();
			if (scaleS) scale.setConstant(scaleT / scaleS);
		}
		for (int dim = 0; dim < 3; dim++) { // prevent weird scaling
			if (fabs(scale[dim]) > 5.0) scale[dim] = 1.0;
		}
		matS = scale.asDiagonal() * matS;

		if (mode == 5) {
			// rotation along major axis
			if (!extractPlanarRotation(matS, matT, rotation)) return false;
		}
		if (mode >= 6) {
			// free rotation
			if (!extractFreeRotation(matS, matT, rotation)) return false;
		}

	} else if(primitive == 2) {

		/////////////// sphere ///////////////

		if (mode == 1 || mode == 3) {
			// uniform scaling
			double scaleT = matMT.norm();
			double scaleS = matMS.norm();
			if(scaleS) scale.setConstant(scaleT / scaleS);
		}
		if (mode == 2 || mode >= 4) {
			// non-uniform scaling
			scale.setOnes();
			for (int dim = 0; dim < 3; dim++) {
				double scaleT, scaleS;
				if (useMassScale) {
					scaleT = matMT.row(dim).norm();
					scaleS = matMS.row(dim).norm();
				} else {
					scaleT = tgtBBSize[dim];
					scaleS = srcBBSize[dim];
				}
				if (scaleS) scale(dim) = scaleT / scaleS;
			}
		}
		for (int dim = 0; dim < 3; dim++) { // prevent weird scaling
			if (fabs(scale[dim]) > 5.0) scale[dim] = 1.0;
		}
		matS = scale.asDiagonal() * matS;

		if (mode >= 3) {
			// free rotation
			if (!extractFreeRotation(matS, matT, rotation)) return false;
		}
	} else {
		cout << "Error: unrecognized primitive " << primitive << endl;
		return false;
	}

	Eigen::Vector3d translation = tgtCenter - rotation * scale.asDiagonal() * srcCenter;
	transformation.setIdentity();
	transformation.prescale(scale);
	transformation.prerotate(rotation);
	transformation.pretranslate(translation);
	
	return true;
}

bool MatchPrimitiveICP::extractFittingMass(
	int mode,
	int primitive,
	Eigen::Matrix3Xd &source,
	Eigen::Matrix3Xd &target,
	Eigen::Affine3d &transformation)
{
	// fit mass distribution

	Eigen::Vector3d srcCenter = source.rowwise().mean();
	Eigen::Vector3d tgtCenter = target.rowwise().mean();
	Eigen::Matrix3Xd matS = source.colwise() - srcCenter;
	Eigen::Matrix3Xd matT = target.colwise() - tgtCenter;
	int numS = (int)matS.cols();
	int numT = (int)matT.cols();

	Eigen::Vector3d scale = Eigen::Vector3d::Ones();

	if (primitive == 0) {

		/////////////// stick ///////////////

		if (mode >= 1 && mode != 4 && mode != 6) {
			// scaling along major axis
			double scaleT = matT.row(1).norm() / numT;
			double scaleS = matS.row(1).norm() / numS;
			if (scaleS) scale[1] = scaleT / scaleS;
		}
		if (mode >= 2 && mode != 4 && mode != 6) {
			// uniform scaling on cross section
			double scaleT = sqrt(matT.row(0).squaredNorm() + matT.row(2).squaredNorm()) / numT;
			double scaleS = sqrt(matS.row(0).squaredNorm() + matS.row(2).squaredNorm()) / numS;
			if (scaleS) scale[0] = scale[2] = scaleT / scaleS;
		}
		if (mode == 4 || mode == 6) {
			// uniform scaling
			double scaleT = matT.norm() / numT;
			double scaleS = matS.norm() / numS;
			if (scaleS) scale.setConstant(scaleT / scaleS);
		}

	} else if (primitive == 1) {

		/////////////// plane ///////////////

		if (mode == 1 || mode == 4 || mode > 6) {
			// scaling along X axis
			double scaleT = matT.row(0).norm() / numT;
			double scaleS = matS.row(0).norm() / numS;
			if (scaleS) scale[0] = scaleT / scaleS;
		}
		if (mode == 2 || mode == 4 || mode > 6) { 
			// scaling along Z axis
			double scaleT = matT.row(2).norm() / numT;
			double scaleS = matS.row(2).norm() / numS;
			if (scaleS) scale[2] = scaleT / scaleS;
		}
		if (mode == 3) {
			// uniform scaling on cross section
			double scaleT = sqrt(matT.row(0).squaredNorm() + matT.row(2).squaredNorm()) / numT;
			double scaleS = sqrt(matS.row(0).squaredNorm() + matS.row(2).squaredNorm()) / numS;
			if (scaleS) scale[0] = scale[2] = scaleT / scaleS;
		}
		if (mode == 5 || mode == 6) {
			// uniform scaling
			double scaleT = matT.norm() / numT;
			double scaleS = matS.norm() / numS;
			if (scaleS) scale.setConstant(scaleT / scaleS);
		}

	} else if (primitive == 2) {

		/////////////// sphere ///////////////

		if (mode == 1 || mode == 3) {
			// uniform scaling
			double scaleT = matT.norm() / numT;
			double scaleS = matS.norm() / numS;
			if (scaleS) scale.setConstant(scaleT / scaleS);
		}
		if (mode == 2 || mode >= 4) {
			// non-uniform scaling
			scale.setOnes();
			for (int dim = 0; dim < 3; dim++) {
				double scaleT = matT.row(dim).norm();
				double scaleS = matS.row(dim).norm();
				if (scaleS) scale(dim) = scaleT / scaleS;
			}
		}
		
	} else {
		cout << "Error: unrecognized primitive " << primitive << endl;
		return false;
	}

	transformation.setIdentity();
	transformation.pretranslate(-srcCenter);
	transformation.prescale(scale);
	transformation.pretranslate(tgtCenter);

	return true;
}

bool MatchPrimitiveICP::extractFittingBoundingBox(
	int mode,
	int primitive,
	Eigen::Matrix3Xd &source,
	Eigen::Matrix3Xd &target,
	Eigen::Affine3d &transformation)
{
	// fit bounding box

	Eigen::AlignedBox3d srcBB(source.rowwise().minCoeff(), source.rowwise().maxCoeff());
	Eigen::AlignedBox3d tgtBB(target.rowwise().minCoeff(), target.rowwise().maxCoeff());

	Eigen::Vector3d srcCenter = srcBB.center();
	Eigen::Vector3d tgtCenter = tgtBB.center();
	Eigen::Vector3d srcBBSize = srcBB.sizes();
	Eigen::Vector3d tgtBBSize = tgtBB.sizes();

	Eigen::Vector3d scale = Eigen::Vector3d::Ones();

	if (primitive == 0) {

		/////////////// stick ///////////////

		if (mode >= 1 && mode != 4 && mode != 6) {
			// scaling along major axis
			double scaleT = tgtBBSize[1];
			double scaleS = srcBBSize[1];
			if (scaleS) scale[1] = scaleT / scaleS;
		}
		if (mode >= 2 && mode != 4 && mode != 6) {
			// uniform scaling on cross section
			double scaleT = sqrt(cml::sqr(tgtBBSize[0]) + cml::sqr(tgtBBSize[2]));
			double scaleS = sqrt(cml::sqr(srcBBSize[0]) + cml::sqr(srcBBSize[2]));
			if (scaleS) scale[0] = scale[2] = scaleT / scaleS;
		}
		if (mode == 4 || mode == 6) {
			// uniform scaling
			double scaleT = tgtBB.diagonal().norm();
			double scaleS = srcBB.diagonal().norm();
			if (scaleS) scale.setConstant(scaleT / scaleS);
		}

	} else if (primitive == 1) {

		/////////////// plane ///////////////

		if (mode == 1 || mode == 4 || mode > 6) {
			// scaling along X axis
			double scaleT = tgtBBSize[0];
			double scaleS = srcBBSize[0];
			if (scaleS) scale[0] = scaleT / scaleS;
		}
		if (mode == 2 || mode == 4 || mode > 6) {
			// scaling along Z axis
			double scaleT = tgtBBSize[2];
			double scaleS = srcBBSize[2];
			if (scaleS) scale[2] = scaleT / scaleS;
		}
		if (mode == 3) {
			// uniform scaling on cross section
			double scaleT = sqrt(cml::sqr(tgtBBSize[0]) + cml::sqr(tgtBBSize[2]));
			double scaleS = sqrt(cml::sqr(srcBBSize[0]) + cml::sqr(srcBBSize[2]));
			if (scaleS) scale[0] = scale[2] = scaleT / scaleS;
		}
		if (mode == 5 || mode == 6) {
			// uniform scaling
			double scaleT = tgtBB.diagonal().norm();
			double scaleS = srcBB.diagonal().norm();
			if (scaleS) scale.setConstant(scaleT / scaleS);
		}

	} else if (primitive == 2) {

		/////////////// sphere ///////////////

		if (mode == 1 || mode == 3) {
			// uniform scaling
			double scaleT = tgtBB.diagonal().norm();
			double scaleS = srcBB.diagonal().norm();
			if (scaleS) scale.setConstant(scaleT / scaleS);
		}
		if (mode == 2 || mode >= 4) {
			// non-uniform scaling
			scale.setOnes();
			for (int dim = 0; dim < 3; dim++) {
				double scaleT = tgtBBSize[dim];
				double scaleS = srcBBSize[dim];
				if (scaleS) scale(dim) = scaleT / scaleS;
			}
		}

	} else {
		cout << "Error: unrecognized primitive " << primitive << endl;
		return false;
	}

	transformation.setIdentity();
	transformation.pretranslate(-srcCenter);
	transformation.prescale(scale);
	transformation.pretranslate(tgtCenter);

	return true;
}

bool MatchPrimitiveICP::extractPlanarRotation(
	Eigen::Matrix3Xd &source,
	Eigen::Matrix3Xd &target,
	Eigen::Matrix3d &rotation)
{

	Eigen::Matrix2Xd projSource(2, source.cols());
	Eigen::Matrix2Xd projTarget(2, target.cols());
	projSource << source.row(0), source.row(2);
	projTarget << target.row(0), target.row(2);

	Eigen::JacobiSVD<Eigen::Matrix2d> svd(projTarget*projSource.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix2d rotation2D = svd.matrixU() * svd.matrixV().transpose();
	rotation.setIdentity();
	rotation(0, 0) = rotation2D(0, 0);
	rotation(0, 2) = rotation2D(0, 1);
	rotation(2, 0) = rotation2D(1, 0);
	rotation(2, 2) = rotation2D(1, 1);

	return true;
}

bool MatchPrimitiveICP::extractFreeRotation(
	Eigen::Matrix3Xd &source,
	Eigen::Matrix3Xd &target,
	Eigen::Matrix3d &rotation)
{
	Eigen::JacobiSVD<Eigen::Matrix3d> svd(target*source.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
	rotation = svd.matrixU() * svd.matrixV().transpose();

	return true;
}

bool MatchPrimitiveICP::visualize(
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