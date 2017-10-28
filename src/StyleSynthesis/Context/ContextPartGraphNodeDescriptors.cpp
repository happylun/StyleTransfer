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

#include "ContextPartGraphNodeDescriptors.h"

#include "Mesh/MeshUtil.h"

#include "Sample/SampleSimplePoissonDisk.h"
#include "Sample/SampleUtil.h"

#include "Feature/FeatureUtil.h"

#include "Utility/PlyExporter.h"

#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

ContextPartGraphNodeDescriptors::ContextPartGraphNodeDescriptors() {

	mNodeSamples.amount = 0;
	mNodeMesh.amount = 0;
}

ContextPartGraphNodeDescriptors::~ContextPartGraphNodeDescriptors() {
}

bool ContextPartGraphNodeDescriptors::preprocess(TTriangleMesh &mesh, vector<int> &segment) {

	if (mNodeSamples.amount) return true; // aleady pre-processed

	if (!MeshUtil::extractSubMesh(mesh, segment, mNodeMesh)) return false;
	if (!MeshUtil::recomputeNormals(mNodeMesh)) return false;

	if (!preprocess()) return false;

	return true;
}

bool ContextPartGraphNodeDescriptors::preprocess() {

	if (mNodeSamples.amount) return true; // aleady pre-processed
	if (mNodeMesh.amount == 0) return false; // node mesh is not established yet

	// sampling

	int numSamples = 1000; // UNDONE: param desired number of samples
	SampleSimplePoissonDisk sspd(&mNodeMesh);
	if (!sspd.runSampling(numSamples)) return false;
	if (!sspd.exportSample(mNodeSamples)) return false;
	numSamples = mNodeSamples.amount;

	// initialize sample matrices

	if (!SampleUtil::buildMatrices(mNodeSamples, mSampleMatP, mSampleMatN)) return false;

	return true;
}

bool ContextPartGraphNodeDescriptors::compute(TTriangleMesh &mesh, int mode) {

	// special handling for different modes (only affects position descriptors)

	Eigen::Matrix3Xd matPosition = mSampleMatP;

	if (mode == 1) {
		// ceiling lamp
		// flip along Y axis
		matPosition.row(1) = -matPosition.row(1);
	}

	if (mode == 2) {
		// wall lamp
		// rotate to make the "wall" look like the "ground"
		Eigen::Matrix3d rot;
		rot << 1, 0, 0, 0, 0, 1, 0, -1, 0;
		matPosition = rot * matPosition;
	}

	if (true) {

		// position

		mMassCenter = matPosition.rowwise().mean();
		mCenterHeight = mMassCenter[1];
		mLowHeight = matPosition.row(1).minCoeff();
		mHighHeight = matPosition.row(1).maxCoeff();
		mRadialDistance = cml::sqrt_safe(cml::sqr(mMassCenter[0]) + cml::sqr(mMassCenter[2]));
	}

	if (true) {

		// size

		if (!MeshUtil::computeFaceArea(mNodeMesh, mMeshArea)) return false;

		mSizeEpsilon = mNodeSamples.radius;

		vec3 bbMin, bbMax;
		if (!MeshUtil::computeAABB(mNodeMesh, bbMin, bbMax)) return false;
		mBoundingBox = Eigen::AlignedBox3d(
			Eigen::Vector3d(vec3d(bbMin).data()),
			Eigen::Vector3d(vec3d(bbMax).data()));

		Eigen::Matrix3Xd matCentered = mSampleMatP;
		matCentered = matCentered.colwise() - matCentered.rowwise().mean();
		mAxisAlignedVariance = matCentered.array().square().rowwise().mean().max(mSizeEpsilon*mSizeEpsilon);

		Eigen::JacobiSVD< Eigen::Matrix3Xd > svd(matCentered, Eigen::ComputeThinU);
		Eigen::Vector3d vecSglrVal = svd.singularValues();
		mPrincipalVariance = vecSglrVal.cwiseMax(mSizeEpsilon*mSizeEpsilon);

		// orientation
				
		Eigen::Matrix3d matU = svd.matrixU();
		Eigen::Matrix3Xd matLocalCS = matU.transpose() * matCentered;
		
		Eigen::Vector3d vecPrimWeights(1.0, 2.0, 3.0);
		auto &primWeights = StyleSynthesisConfig::mContext_GraphNodePrimitiveAspectWeights.values;
		if(primWeights.size() >= 3) vecPrimWeights = Eigen::Vector3d(primWeights[0], primWeights[1], primWeights[2]);
		Eigen::Vector3d vecExtent = matLocalCS.rowwise().maxCoeff() - matLocalCS.rowwise().minCoeff();
		Eigen::Vector3d vecAspect(vecExtent[0] - vecExtent[1], vecExtent[1] - vecExtent[2], vecExtent[2]);
		vecAspect = vecAspect.cwiseProduct(vecPrimWeights);
		if (vecAspect[0] > vecAspect[1] && vecAspect[0] > vecAspect[2]) {
			// stick-like
			mPrimitive = 0;
			mMajorOrientation = matU.col(0);
		} else if (vecAspect[1] > vecAspect[0] && vecAspect[1] > vecAspect[2]) {
			// plane-like
			mPrimitive = 1;
			mMajorOrientation = matU.col(2);
		} else {
			// sphere-like
			mPrimitive = 2;
			mMajorOrientation = matU.col(0);
		}
		mMajorOrientation.normalize();
	}

	if (true) {

		// distribution

		if (!computeMassDistributions(mesh, mNodeMesh, mMassDistribution)) return false;
	}

	return true;
}

bool ContextPartGraphNodeDescriptors::computeMassDistributions(TTriangleMesh &wholeMesh, TTriangleMesh &nodeMesh, vector<double> &histogram) {

	// mass distribution by counting "inside" points

	int numDenseGrids = NUM_MASS_DIST_GRIDS * 8; // UNDONE: param number of dense grids for mass distribution

	TKDTree nodeTree;
	TKDTreeData nodeTreeData;
	if (!MeshUtil::buildKdTree(nodeMesh, nodeTree, nodeTreeData)) return false;

	vec3 meshBBMin, meshBBMax;
	if (!MeshUtil::computeAABB(wholeMesh, meshBBMin, meshBBMax)) return false;
	vec3 denseGridSize = (meshBBMax - meshBBMin) * (1.0f / numDenseGrids);

	vector<vec3i> gridIndex(0); // (xID, yID, zID) : # of dense grids
	vector<int> gridMap(0); // output grid ID : # of dense grids
	for (int xID = 0; xID < numDenseGrids; xID++) {
		for (int yID = 0; yID < numDenseGrids; yID++) {
			for (int zID = 0; zID < numDenseGrids; zID++) {
				gridIndex.push_back(vec3i(xID, yID, zID));
				int binXID = xID*NUM_MASS_DIST_GRIDS / numDenseGrids;
				int binYID = yID*NUM_MASS_DIST_GRIDS / numDenseGrids;
				int binZID = zID*NUM_MASS_DIST_GRIDS / numDenseGrids;
				int binID = (binXID * NUM_MASS_DIST_GRIDS + binYID) * NUM_MASS_DIST_GRIDS + binZID;
				gridMap.push_back(binID);
			}
		}
	}

	int numTotalPoints = (int)gridIndex.size();
	vector<bool> gridInsideFlag(numTotalPoints, false);
	vector<vec3> gridPos(numTotalPoints);
#pragma omp parallel for
	for (int pointID = 0; pointID < numTotalPoints; pointID++) {

		vec3i gridID = gridIndex[pointID];
		gridPos[pointID] = meshBBMin + vec3(
			denseGridSize[0] * (gridID[0] + 0.5f),
			denseGridSize[1] * (gridID[1] + 0.5f),
			denseGridSize[2] * (gridID[2] + 0.5f));
		bool insideFlag = false;
		if (!MeshUtil::checkInsideMesh(nodeTree, gridPos[pointID], insideFlag)) error("check inside mesh");
		gridInsideFlag[pointID] = insideFlag;
	}

	int numHistogramBins = NUM_MASS_DIST_GRIDS * NUM_MASS_DIST_GRIDS * NUM_MASS_DIST_GRIDS;
	vector<int> allHistCount(numHistogramBins, 0);
	for (int pointID = 0; pointID < numTotalPoints; pointID++) {
		if (gridInsideFlag[pointID]) {
			int binID = gridMap[pointID];
			allHistCount[binID]++;
		}
	}

	histogram.resize(numHistogramBins);
	for (int binID = 0; binID < numHistogramBins; binID++) {
		histogram[binID] = allHistCount[binID] / (double)numTotalPoints;
	}

	/*
	if (true) {
		// visualization
		PlyExporter pe;
		vector<vec3i> colors(numTotalPoints, vec3i(127, 127, 127));
		for (int pointID = 0; pointID < numTotalPoints; pointID++) {
			if (gridInsideFlag[pointID]) colors[pointID] = vec3i(255, 0, 0);
		}
		if (!pe.addPoint(&gridPos, 0, &colors)) return false;
		if (!pe.output("massDist.ply")) return false;
		if (!MeshUtil::saveMesh("massShape.ply", nodeMesh)) return false;
		system("pause");
	}
	*/

	return true;
}

bool ContextPartGraphNodeDescriptors::computeShapeDistributions(TSampleSet &samples, vector<double> &histogram) {

	// simple shape distribution

	int numBins = 64; // UNDONE: param SD bins
	vector<double> distribution(0);
	distribution.reserve(cml::sqr(samples.amount) / 2);
	double maxDist = 0;
	double minDist = DBL_MAX;
	for (int i = 0; i < samples.amount - 1; i++) {
		vec3 pi = samples.positions[i];
		for (int j = i + 1; j < samples.amount; j++) {
			vec3 pj = samples.positions[j];
			double dist = (double)(pi - pj).length();
			distribution.push_back(dist);
			maxDist = max(maxDist, dist);
			minDist = min(minDist, dist);
		}
	}
	histogram.assign(numBins, 1.0 / numBins);
	if (!FeatureUtil::computeHistogram(distribution, histogram, numBins, minDist, maxDist)) return false;

	return true;
}

bool ContextPartGraphNodeDescriptors::saveData(ostream &fileStream) {

	// descriptors

	fileStream.write((char*)mMassCenter.data(), sizeof(Eigen::Vector3d));
	fileStream.write((char*)&mCenterHeight, sizeof(mCenterHeight));
	fileStream.write((char*)&mLowHeight, sizeof(mLowHeight));
	fileStream.write((char*)&mHighHeight, sizeof(mHighHeight));
	fileStream.write((char*)&mRadialDistance, sizeof(mRadialDistance));
	fileStream.write((char*)&mSizeEpsilon, sizeof(mSizeEpsilon));
	fileStream.write((char*)&mMeshArea, sizeof(mMeshArea));
	fileStream.write((char*)mBoundingBox.min().data(), sizeof(Eigen::Vector3d));
	fileStream.write((char*)mBoundingBox.max().data(), sizeof(Eigen::Vector3d));
	fileStream.write((char*)mAxisAlignedVariance.data(), sizeof(Eigen::Vector3d));
	fileStream.write((char*)mPrincipalVariance.data(), sizeof(Eigen::Vector3d));
	fileStream.write((char*)mMajorOrientation.data(), sizeof(Eigen::Vector3d));
	fileStream.write((char*)&mPrimitive, sizeof(mPrimitive));

	int numHistogramBins = NUM_MASS_DIST_GRIDS * NUM_MASS_DIST_GRIDS * NUM_MASS_DIST_GRIDS;
	for (int binID = 0; binID < numHistogramBins; binID++) {
		fileStream.write((char*)&mMassDistribution[binID], sizeof(mMassDistribution[binID]));
	}

	return true;
}

bool ContextPartGraphNodeDescriptors::loadData(istream &fileStream) {

	// descriptors

	fileStream.read((char*)mMassCenter.data(), sizeof(Eigen::Vector3d));
	fileStream.read((char*)&mCenterHeight, sizeof(mCenterHeight));
	fileStream.read((char*)&mLowHeight, sizeof(mLowHeight));
	fileStream.read((char*)&mHighHeight, sizeof(mHighHeight));
	fileStream.read((char*)&mRadialDistance, sizeof(mRadialDistance));
	fileStream.read((char*)&mSizeEpsilon, sizeof(mSizeEpsilon));
	fileStream.read((char*)&mMeshArea, sizeof(mMeshArea));

	Eigen::Vector3d bbMin, bbMax;
	fileStream.read((char*)bbMin.data(), sizeof(Eigen::Vector3d));
	fileStream.read((char*)bbMax.data(), sizeof(Eigen::Vector3d));
	mBoundingBox.setEmpty();
	mBoundingBox.extend(bbMin);
	mBoundingBox.extend(bbMax);

	fileStream.read((char*)mAxisAlignedVariance.data(), sizeof(Eigen::Vector3d));
	fileStream.read((char*)mPrincipalVariance.data(), sizeof(Eigen::Vector3d));
	fileStream.read((char*)mMajorOrientation.data(), sizeof(Eigen::Vector3d));
	fileStream.read((char*)&mPrimitive, sizeof(mPrimitive));

	int numHistogramBins = NUM_MASS_DIST_GRIDS * NUM_MASS_DIST_GRIDS * NUM_MASS_DIST_GRIDS;
	mMassDistribution.resize(numHistogramBins);
	for (int binID = 0; binID < numHistogramBins; binID++) {
		fileStream.read((char*)&mMassDistribution[binID], sizeof(mMassDistribution[binID]));
	}

	return true;
}