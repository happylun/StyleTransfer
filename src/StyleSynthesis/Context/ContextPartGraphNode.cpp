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

#include "ContextPartGraphNode.h"

#include "Mesh/MeshUtil.h"
#include "Sample/SampleUtil.h"
#include "Feature/FeatureUtil.h"

#include "Sample/SampleSimplePoissonDisk.h"

#include "Match/MatchRigidICP.h"
#include "Match/MatchSimpleICP.h"

#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

ContextPartGraphNode::ContextPartGraphNode() {
	clearNode();
}

ContextPartGraphNode::~ContextPartGraphNode() {
}

bool ContextPartGraphNode::clearNode() {

	mParent = 0;
	mChildren.clear();
	mSymmetry.clear();
	mCocentric.clear();
	mAdjacent.clear();
	mContact.clear();
	mSupport.clear();

	mUID = -1;

	return true;
}

bool ContextPartGraphNode::computeDescriptor(int mode) {

	if (!mNodeDescriptors.preprocess(*mpGraphMesh, (*mpGraphSegments)[mPartLevelID][mPartSegmentID])) return false;
	if (!mNodeDescriptors.compute(*mpGraphMesh, mode)) return false;

	return true;
}

bool ContextPartGraphNode::copyMetaData(ContextPartGraphNode *otherNode) {

	if (!clearNode()) return false;

	mpGraphMesh = otherNode->mpGraphMesh;
	mpGraphSegments = otherNode->mpGraphSegments;
	mPartLevelID = otherNode->mPartLevelID;
	mPartSegmentID = otherNode->mPartSegmentID;

	return true;
}

bool ContextPartGraphNode::copyNode(ContextPartGraphNode *otherNode) {

	if (!copyMetaData(otherNode)) return false;

	mNodeDescriptors = otherNode->mNodeDescriptors;

	return true;
}

bool ContextPartGraphNode::transformNode(ContextPartGraphNode *otherNode, Eigen::Affine3d &transformation) {
	
	if (!copyMetaData(otherNode)) return false;

	// get other node mesh

	if (!otherNode->mNodeDescriptors.preprocess(*otherNode->mpGraphMesh,
		(*otherNode->mpGraphSegments)[otherNode->mPartLevelID][otherNode->mPartSegmentID])) return false;
	TTriangleMesh &otherNodeMesh = otherNode->mNodeDescriptors.mNodeMesh;

	// transform mesh

	TTriangleMesh &thisNodeMesh = mNodeDescriptors.mNodeMesh;
	thisNodeMesh = otherNodeMesh;
	Eigen::Matrix3d rotation = transformation.rotation();
	for (vec3 &position : thisNodeMesh.positions) {
		Eigen::Vector3d vec(vec3d(position).data());
		vec = transformation * vec;
		position = vec3d(vec[0], vec[1], vec[2]);
	}
	for (vec3 &normal : thisNodeMesh.normals) {
		Eigen::Vector3d vec(vec3d(normal).data());
		vec = rotation * vec;
		normal = vec3d(vec[0], vec[1], vec[2]);
	}

	// update node descriptor

	if (!mNodeDescriptors.preprocess()) return false;
	if (!mNodeDescriptors.compute(*mpGraphMesh)) return false;

	return true;
}

bool ContextPartGraphNode::detectSymmetry(ContextPartGraphNode *node1, ContextPartGraphNode *node2) {

	if (StyleSynthesisConfig::mContext_SymmetryCheckingByVertex) {
		return detectSymmetryRigid(node1, node2) || detectSymmetryRigid(node2, node1);
	} else {
		return detectSymmetryNonUniform(node1, node2) || detectSymmetryNonUniform(node2, node1);
	}
	return true;
}

bool ContextPartGraphNode::detectSymmetryRigid(ContextPartGraphNode *node1, ContextPartGraphNode *node2) {

	double errorFactor = StyleSynthesisConfig::mContext_SymmetryCheckingThreshold;

	auto &desc1 = node1->mNodeDescriptors;
	auto &desc2 = node2->mNodeDescriptors;

	// quick pruning by height
	if (StyleSynthesisConfig::mContext_SymmetryCheckingPruneByHeight) {
		double eps = max(desc1.mSizeEpsilon, desc2.mSizeEpsilon) * 1.0;
		double diffHeight = fabs(desc1.mCenterHeight - desc2.mCenterHeight);
		if (diffHeight > eps) return false;
	}

	// quick pruning by PCA variance
	if (true) {
		Eigen::Vector3d &v1 = desc1.mPrincipalVariance;
		Eigen::Vector3d &v2 = desc2.mPrincipalVariance;
		Eigen::Vector3d relVar1 = v1.cwiseQuotient(v2);
		Eigen::Vector3d relVar2 = v2.cwiseQuotient(v1);
		double error = relVar1.cwiseMax(relVar2).sum() - 3.0;
		if (error > 0.5) return false; // UNDONE: param variance difference threshold for pruning
	}

	// compare mesh vertices
	Eigen::Matrix3Xd matSP, matSN, matTP, matTN;
	if (!SampleUtil::buildMatrices(node1->mNodeDescriptors.mNodeMesh, matSP, matSN)) return false;
	if (!SampleUtil::buildMatrices(node2->mNodeDescriptors.mNodeMesh, matTP, matTN)) return false;
	
	Eigen::Affine3d xform;
	if (!MatchRigidICP::runRegularShape(10, matSP, matSN, matTP, matTN, xform)) return false;
	//if (!MatchRigidICP::visualize("regularICP.ply", matSP, matTP, xform)) return false;

	double errorST;
	Eigen::Matrix3Xd matXSP = xform * matSP;
	if (!MatchRigidICP::error(matXSP, matTP, errorST)) return false;

	double errorTS;
	Eigen::Matrix3Xd matXTP = xform.inverse() * matTP;
	if (!MatchRigidICP::error(matXTP, matSP, errorTS)) return false;

	if (errorST < cml::sqr(desc2.mSizeEpsilon * errorFactor) &&
		errorTS < cml::sqr(desc1.mSizeEpsilon * errorFactor)) return true;

	return false;
}

bool ContextPartGraphNode::detectSymmetryNonUniform(ContextPartGraphNode *node1, ContextPartGraphNode *node2) {

	double errorFactor = StyleSynthesisConfig::mContext_SymmetryCheckingThreshold;

	auto &desc1 = node1->mNodeDescriptors;
	auto &desc2 = node2->mNodeDescriptors;

	// quick pruning by height
	if (StyleSynthesisConfig::mContext_SymmetryCheckingPruneByHeight) {
		double eps = max(desc1.mSizeEpsilon, desc2.mSizeEpsilon) * 1.0;
		double diffHeight = fabs(desc1.mCenterHeight - desc2.mCenterHeight);
		if (diffHeight > eps) return false;
	}

	// compare sample points
	Eigen::Matrix3Xd &matSP = desc1.mSampleMatP;
	Eigen::Matrix3Xd &matSN = desc1.mSampleMatN;
	Eigen::Matrix3Xd &matTP = desc2.mSampleMatP;
	Eigen::Matrix3Xd &matTN = desc2.mSampleMatN;

	//// compare mesh vertices
	//Eigen::Matrix3Xd matSP, matSN, matTP, matTN;
	//if (!SampleUtil::buildMatrices(node1->mNodeDescriptors.mNodeMesh, matSP, matSN)) return false;
	//if (!SampleUtil::buildMatrices(node2->mNodeDescriptors.mNodeMesh, matTP, matTN)) return false;

	Eigen::Affine3d xform;
	if (!MatchSimpleICP::runRegularShape(10, matSP, matSN, matTP, matTN, xform)) return false;
	//if (!MatchRigidICP::visualize("regularICP.ply", matSP, matTP, xform)) return false;

	// check scaling factor

	if (true) {
		Eigen::Matrix3d linearMat = xform.linear();
		Eigen::JacobiSVD<Eigen::Matrix3d> svd(linearMat);
		Eigen::Vector3d scaleFactor = svd.singularValues();
		if (scaleFactor.maxCoeff() > 5.0 || scaleFactor.minCoeff() < 0.2) { // weird scaling
			return false;
		}
	}

	// check point-to-point distance

	double errorST;
	Eigen::Matrix3Xd matXSP = xform * matSP;
	if (!MatchSimpleICP::error(matXSP, matTP, errorST)) return false;

	double errorTS;
	Eigen::Matrix3Xd matXTP = xform.inverse() * matTP;
	if (!MatchSimpleICP::error(matXTP, matSP, errorTS)) return false;

	if (errorST < cml::sqr(desc2.mSizeEpsilon * errorFactor) &&
		errorTS < cml::sqr(desc1.mSizeEpsilon * errorFactor)) return true;

	return false;
}

bool ContextPartGraphNode::detectCoCentricity(ContextPartGraphNode *node1, ContextPartGraphNode *node2) {

	double errorFactor = 2.0; // UNDONE: param "co-centric" error factor

	auto &desc1 = node1->mNodeDescriptors;
	auto &desc2 = node2->mNodeDescriptors;

	Eigen::Vector3d center1 = desc1.mBoundingBox.center();
	Eigen::Vector3d center2 = desc2.mBoundingBox.center();
	Eigen::Vector3d axisDir = center2 - center1;
	if (axisDir.norm() > 0) {
		axisDir.normalize();
	} else {
		// use axis with largest variance
		int index;
		desc1.mAxisAlignedVariance.maxCoeff(&index, (int*)0);
		axisDir = Eigen::Vector3d(0, 0, 0);
		axisDir[index] = 1.0;
	}

	Eigen::AngleAxisd rotation(cml::rad(180.0), axisDir);

	Eigen::Affine3d trans1;
	trans1.setIdentity();
	trans1.pretranslate(-center1);
	trans1.prerotate(rotation);
	trans1.pretranslate(center1);

	double error1;
	Eigen::Matrix3Xd matXP1 = trans1 * desc1.mSampleMatP;
	if (!MatchRigidICP::error(matXP1, desc1.mSampleMatP, error1)) return false;
	if (error1 > cml::sqr(desc1.mSizeEpsilon * errorFactor)) return false;

	Eigen::Affine3d trans2;
	trans2.setIdentity();
	trans2.pretranslate(-center2);
	trans2.prerotate(rotation);
	trans2.pretranslate(center2);

	double error2;
	Eigen::Matrix3Xd matXP2 = trans2 * desc2.mSampleMatP;
	if (!MatchRigidICP::error(matXP2, desc2.mSampleMatP, error2)) return false;
	if (error2 > cml::sqr(desc2.mSizeEpsilon * errorFactor)) return false;

	return true;
}

bool ContextPartGraphNode::detectAdjacency(ContextPartGraphNode *node1, ContextPartGraphNode *node2) {

	auto &desc1 = node1->mNodeDescriptors;
	auto &desc2 = node2->mNodeDescriptors;

	Eigen::AlignedBox3d bb1 = desc1.mBoundingBox;
	Eigen::AlignedBox3d bb2 = desc2.mBoundingBox;
	double offset1 = desc1.mSizeEpsilon * 1; // expand BB a bit
	double offset2 = desc2.mSizeEpsilon * 1;
	double szFactor = 0.1;
	Eigen::Vector3d minExt1 = bb1.min() - (bb1.sizes()*szFactor).cwiseMin(offset1);
	Eigen::Vector3d maxExt1 = bb1.max() + (bb1.sizes()*szFactor).cwiseMin(offset1);
	Eigen::Vector3d minExt2 = bb2.min() - (bb2.sizes()*szFactor).cwiseMin(offset2);
	Eigen::Vector3d maxExt2 = bb2.max() + (bb2.sizes()*szFactor).cwiseMin(offset2);
	bb1.extend(minExt1);
	bb1.extend(maxExt1);
	bb2.extend(minExt2);
	bb2.extend(maxExt2);

	if (bb1.intersects(bb2)) return true;

	return false;
}

bool ContextPartGraphNode::detectContact(ContextPartGraphNode *node1, ContextPartGraphNode *node2) {

	// assumes adjacent

	auto &desc1 = node1->mNodeDescriptors;
	auto &desc2 = node2->mNodeDescriptors;

	double threshold = desc1.mSizeEpsilon + desc2.mSizeEpsilon;
	Eigen::AlignedBox3d &bb1 = desc1.mBoundingBox;
	Eigen::AlignedBox3d &bb2 = desc2.mBoundingBox;
	if (fabs(bb1.max()[0] - bb2.min()[0]) < threshold ||
		fabs(bb1.min()[0] - bb2.max()[0]) < threshold ||
		fabs(bb1.max()[2] - bb2.min()[2]) < threshold ||
		fabs(bb1.min()[2] - bb2.max()[2]) < threshold) return true;

	return false;
}

bool ContextPartGraphNode::detectSupport(ContextPartGraphNode *node1, ContextPartGraphNode *node2) {

	// assumes adjacent

	auto &desc1 = node1->mNodeDescriptors;
	auto &desc2 = node2->mNodeDescriptors;

	double threshold = desc1.mSizeEpsilon + desc2.mSizeEpsilon;
	Eigen::AlignedBox3d &bb1 = desc1.mBoundingBox;
	Eigen::AlignedBox3d &bb2 = desc2.mBoundingBox;
	if (fabs(bb1.max()[1] - bb2.min()[1]) < threshold ||
		fabs(bb1.min()[1] - bb2.max()[1]) < threshold) return true;

	return false;
}

bool ContextPartGraphNode::computeNodeDistance(
	ContextPartGraphNode *node1, ContextPartGraphNode *node2,
	vector<double> &distance)
{
	auto &desc1 = node1->mNodeDescriptors;
	auto &desc2 = node2->mNodeDescriptors;

	if (StyleSynthesisConfig::mContext_UseLagaDescriptors) {
		// special handling: override with Laga's descriptors

#pragma omp critical
		{
			bool fail = false;
			if (!desc1.preprocess(*node1->mpGraphMesh, (*node1->mpGraphSegments)[node1->mPartLevelID][node1->mPartSegmentID])) fail = true;
			if (desc1.mShapeDistributions.empty() && !desc1.computeShapeDistributions(desc1.mNodeSamples, desc1.mShapeDistributions)) fail = true;

			if (!desc2.preprocess(*node2->mpGraphMesh, (*node2->mpGraphSegments)[node2->mPartLevelID][node2->mPartSegmentID])) fail = true;
			if (desc2.mShapeDistributions.empty() && !desc2.computeShapeDistributions(desc2.mNodeSamples, desc2.mShapeDistributions)) fail = true;
			if (fail) {
				cout << "Error: fail to preprocess node descriptors for Laga's features" << endl;
				system("pause");
			}
		}

		// shape descriptor
		double distD = FeatureUtil::computeEMD(desc1.mShapeDistributions, desc2.mShapeDistributions);
		// component size
		double distS = fabs(desc1.mBoundingBox.diagonal().norm() - desc2.mBoundingBox.diagonal().norm());
		// component aspect
		double distA = (desc1.mPrincipalVariance - desc2.mPrincipalVariance).norm();

		distance.clear();
		distance.push_back(distD);
		distance.push_back(distS);
		distance.push_back(distA);

		return true;
	}

	double distCenter = (desc1.mMassCenter - desc2.mMassCenter).norm();
	double distCenterHeight = fabs(desc1.mCenterHeight - desc2.mCenterHeight);
	double distLowHeight = fabs(desc1.mLowHeight - desc2.mLowHeight);
	double distHighHeight = fabs(desc1.mHighHeight - desc2.mHighHeight);
	double distRadius = fabs(desc1.mRadialDistance - desc2.mRadialDistance);

	Eigen::Vector3d bbSize1 = desc1.mBoundingBox.sizes().cwiseMax(desc1.mSizeEpsilon);
	Eigen::Vector3d bbSize2 = desc2.mBoundingBox.sizes().cwiseMax(desc2.mSizeEpsilon);
	Eigen::Vector3d quoSize = bbSize1.cwiseQuotient(bbSize2);
	Eigen::Vector3d quoVariance = desc1.mAxisAlignedVariance.cwiseQuotient(desc2.mAxisAlignedVariance);
	Eigen::Vector3d distSize = quoSize.cwiseMax(quoSize.cwiseInverse()) - Eigen::Vector3d::Ones();
	Eigen::Vector3d distVariance = quoVariance.cwiseMax(quoVariance.cwiseInverse()) - Eigen::Vector3d::Ones();

	double distOrientation = 1.0 - min(1.0, fabs(desc1.mMajorOrientation.dot(desc2.mMajorOrientation)));

	double distMassDist = FeatureUtil::computeL1D(desc1.mMassDistribution, desc2.mMassDistribution);

	distance.clear();
	distance.push_back(distCenter);
	distance.push_back(distCenterHeight);
	distance.push_back(distLowHeight);
	distance.push_back(distHighHeight);
	distance.push_back(distRadius);
	for (int dim = 0; dim < 3; dim++) distance.push_back(distSize[dim]);
	for (int dim = 0; dim < 3; dim++) distance.push_back(distVariance[dim]);
	distance.push_back(distOrientation);
	distance.push_back(distMassDist);

	return true;
}

bool ContextPartGraphNode::computeEdgeDistance(
	ContextPartGraphNode *node11, ContextPartGraphNode *node12,
	ContextPartGraphNode *node21, ContextPartGraphNode *node22,
	vector<double> &distance)
{
	auto &desc11 = node11->mNodeDescriptors;
	auto &desc12 = node12->mNodeDescriptors;
	auto &desc21 = node21->mNodeDescriptors;
	auto &desc22 = node22->mNodeDescriptors;

	if (StyleSynthesisConfig::mContext_UseLagaDescriptors) {
		// special handling: override with Laga's descriptors

		distance.clear();
		distance.push_back(0);

		return true;
	}

	Eigen::Vector3d relativeCenter1 = desc12.mMassCenter - desc11.mMassCenter;
	Eigen::Vector3d relativeCenter2 = desc22.mMassCenter - desc21.mMassCenter;
	double distCenter = (relativeCenter1 - relativeCenter2).norm();

	double relativeCenterHeight1 = desc12.mCenterHeight - desc11.mCenterHeight;
	double relativeCenterHeight2 = desc22.mCenterHeight - desc21.mCenterHeight;
	double relativeLowHeight1 = desc12.mLowHeight - desc11.mLowHeight;
	double relativeLowHeight2 = desc22.mLowHeight - desc21.mLowHeight;
	double relativeHighHeight1 = desc12.mHighHeight - desc11.mHighHeight;
	double relativeHighHeight2 = desc22.mHighHeight - desc21.mHighHeight;
	double distCenterHeight = fabs(relativeCenterHeight1 - relativeCenterHeight2);
	double distLowHeight = fabs(relativeLowHeight1 - relativeLowHeight2);
	double distHighHeight = fabs(relativeHighHeight1 - relativeHighHeight2);

	double relativeRadius1 = desc12.mRadialDistance - desc11.mRadialDistance;
	double relativeRadius2 = desc22.mRadialDistance - desc21.mRadialDistance;
	double distRadius = fabs(relativeRadius1 - relativeRadius2);

	Eigen::Vector3d bbSize11 = desc11.mBoundingBox.sizes().cwiseMax(desc11.mSizeEpsilon);
	Eigen::Vector3d bbSize12 = desc12.mBoundingBox.sizes().cwiseMax(desc12.mSizeEpsilon);
	Eigen::Vector3d bbSize21 = desc21.mBoundingBox.sizes().cwiseMax(desc21.mSizeEpsilon);
	Eigen::Vector3d bbSize22 = desc22.mBoundingBox.sizes().cwiseMax(desc22.mSizeEpsilon);
	Eigen::Vector3d relSize1 = bbSize12.cwiseQuotient(bbSize11);
	Eigen::Vector3d relSize2 = bbSize22.cwiseQuotient(bbSize21);
	Eigen::Vector3d relVariance1 = desc12.mAxisAlignedVariance.cwiseQuotient(desc11.mAxisAlignedVariance);
	Eigen::Vector3d relVariance2 = desc22.mAxisAlignedVariance.cwiseQuotient(desc21.mAxisAlignedVariance);
	Eigen::Vector3d quoSize = relSize1.cwiseQuotient(relSize2);
	Eigen::Vector3d quoVariance = relVariance1.cwiseQuotient(relVariance2);
	Eigen::Vector3d distSize = quoSize.cwiseMax(quoSize.cwiseInverse()) - Eigen::Vector3d::Ones();
	Eigen::Vector3d distVariance = quoVariance.cwiseMax(quoVariance.cwiseInverse()) - Eigen::Vector3d::Ones();

	double relativeAngle1 = cml::acos_safe(fabs(desc11.mMajorOrientation.dot(desc12.mMajorOrientation)));
	double relativeAngle2 = cml::acos_safe(fabs(desc21.mMajorOrientation.dot(desc22.mMajorOrientation)));
	double distOrientation = 1.0 - cos(relativeAngle1 - relativeAngle2);

	// NOTE: how to compute RELATIVE mass distribution distance ?
	vector<double> relativeMassDist1 = desc11.mMassDistribution;	
	vector<double> relativeMassDist2 = desc21.mMassDistribution;
	for (int id = 0; id < (int)relativeMassDist1.size(); id++) relativeMassDist1[id] -= desc12.mMassDistribution[id];
	for (int id = 0; id < (int)relativeMassDist2.size(); id++) relativeMassDist2[id] -= desc22.mMassDistribution[id];
	double distMassDist = FeatureUtil::computeL1D(relativeMassDist1, relativeMassDist2);

	distance.clear();
	distance.push_back(distCenter);
	distance.push_back(distCenterHeight);
	distance.push_back(distLowHeight);
	distance.push_back(distHighHeight);
	distance.push_back(distRadius);
	for (int dim = 0; dim < 3; dim++) distance.push_back(distSize[dim]);
	for (int dim = 0; dim < 3; dim++) distance.push_back(distVariance[dim]);
	distance.push_back(distOrientation);
	distance.push_back(distMassDist);

	return true;
}

bool ContextPartGraphNode::computeSimilarity(
	vector<double> &distance,
	vector<double> &sigma,
	vector<double> &sigmaMultipliers,
	vector<double> &weights,
	double &similarity)
{

	int numTerms = (int)distance.size();
	if ((int)weights.size() != numTerms) {
		cout << "Error: incorrect dimension for weights" << endl;
		return false;
	}

	double sumWeightedValues = 0;
	//double sumWeights = 0;
	for (int dim = 0; dim < numTerms; dim++) {
		sumWeightedValues += weights[dim] * exp( -cml::sqr(distance[dim] / sigma[dim]) * sigmaMultipliers[dim] ); // multiplier outside square
		//sumWeights += weights[dim];
	}
	similarity = sumWeightedValues; // sumWeights;

	return true;
}

bool ContextPartGraphNode::saveHierarchy(ostream &fileStream) {

	// meta data
	fileStream.write((char*)&mPartLevelID, sizeof(mPartLevelID));
	fileStream.write((char*)&mPartSegmentID, sizeof(mPartSegmentID));

	// parent
	fileStream.write((char*)&mParent->mUID, sizeof(int));

	// children
	int numChildren = (int)mChildren.size();
	fileStream.write((char*)&numChildren, sizeof(int));
	for (auto &node : mChildren) fileStream.write((char*)&node->mUID, sizeof(int));

	return true;
}

bool ContextPartGraphNode::loadHierarchy(istream &fileStream, vector<ContextPartGraphNode*> *allNodes, ContextPartGraphNode *rootNode) {

	// meta data
	fileStream.read((char*)&mPartLevelID, sizeof(mPartLevelID));
	fileStream.read((char*)&mPartSegmentID, sizeof(mPartSegmentID));

	// parent
	int parentID;
	fileStream.read((char*)&parentID, sizeof(int));
	mParent = parentID >= 0 ? (*allNodes)[parentID] : rootNode;
	if (parentID < 0) rootNode->mChildren.push_back(this);

	// children
	int numChildren;
	fileStream.read((char*)&numChildren, sizeof(int));
	mChildren.resize(numChildren);
	for (auto &node : mChildren) {
		int childID;
		fileStream.read((char*)&childID, sizeof(int));
		node = (*allNodes)[childID];
	}

	return true;
}

bool ContextPartGraphNode::saveDescriptor(ostream &fileStream) {

	if (!mNodeDescriptors.saveData(fileStream)) return false;

	return true;
}

bool ContextPartGraphNode::loadDescriptor(istream &fileStream) {

	if (!mNodeDescriptors.loadData(fileStream)) return false;

	return true;
}

bool ContextPartGraphNode::saveContext(ostream &fileStream) {

	// symmetry
	int numSymmetry = (int)mSymmetry.size();
	fileStream.write((char*)&numSymmetry, sizeof(int));
	for (auto &node : mSymmetry) fileStream.write((char*)&node->mUID, sizeof(int));

	// co-centric
	int numCocentric = (int)mCocentric.size();
	fileStream.write((char*)&numCocentric, sizeof(int));
	for (auto &node : mCocentric) fileStream.write((char*)&node->mUID, sizeof(int));

	// adjacent
	int numAdjacent = (int)mAdjacent.size();
	fileStream.write((char*)&numAdjacent, sizeof(int));
	for (auto &node : mAdjacent) fileStream.write((char*)&node->mUID, sizeof(int));

	// contact
	int numContact = (int)mContact.size();
	fileStream.write((char*)&numContact, sizeof(int));
	for (auto &node : mContact) fileStream.write((char*)&node->mUID, sizeof(int));

	// support
	int numSupport = (int)mSupport.size();
	fileStream.write((char*)&numSupport, sizeof(int));
	for (auto &node : mSupport) fileStream.write((char*)&node->mUID, sizeof(int));

	return true;
}

bool ContextPartGraphNode::loadContext(istream &fileStream, vector<ContextPartGraphNode*> *allNodes) {

	// symmetry
	int numSymmetry;
	fileStream.read((char*)&numSymmetry, sizeof(int));
	mSymmetry.resize(numSymmetry);
	for (auto &node : mSymmetry) {
		int symmetryID;
		fileStream.read((char*)&symmetryID, sizeof(int));
		node = (*allNodes)[symmetryID];
	}

	// co-centric
	int numCocentric;
	fileStream.read((char*)&numCocentric, sizeof(int));
	mCocentric.resize(numCocentric);
	for (auto &node : mCocentric) {
		int cocentricID;
		fileStream.read((char*)&cocentricID, sizeof(int));
		node = (*allNodes)[cocentricID];
	}

	// adjacent
	int numAdjacent;
	fileStream.read((char*)&numAdjacent, sizeof(int));
	mAdjacent.resize(numAdjacent);
	for (auto &node : mAdjacent) {
		int adjacentID;
		fileStream.read((char*)&adjacentID, sizeof(int));
		node = (*allNodes)[adjacentID];
	}

	// contact
	int numContact;
	fileStream.read((char*)&numContact, sizeof(int));
	mContact.resize(numContact);
	for (auto &node : mContact) {
		int contactID;
		fileStream.read((char*)&contactID, sizeof(int));
		node = (*allNodes)[contactID];
	}

	// support
	int numSupport;
	fileStream.read((char*)&numSupport, sizeof(int));
	mSupport.resize(numSupport);
	for (auto &node : mSupport) {
		int supportID;
		fileStream.read((char*)&supportID, sizeof(int));
		node = (*allNodes)[supportID];
	}

	return true;
}