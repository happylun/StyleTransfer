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

#include "SimilarityDistance.h"

#include <fstream>
#include <set>

#include "Mesh/MeshUtil.h"

#include "Sample/SampleSimplePoissonDisk.h"
#include "Sample/SampleUtil.h"

#include "Match/MatchSimpleICP.h"

#include "Similarity/SimilarityMetric.h"

#include "Utility/PlyExporter.h"

#include "Data/StyleSynthesisConfig.h"

#define OUTPUT_PROGRESS

using namespace StyleSynthesis;

SimilarityDistance::SimilarityDistance(SimilarityData *data) {

	mpData = data;
	mDistanceData.mElementDistance.resize(0, 0);
}

SimilarityDistance::~SimilarityDistance() {
}

bool SimilarityDistance::process() {

	if (!computeSegmentData()) return false;
	if (!computeElementDistance()) return false;

	return true;
}

bool SimilarityDistance::computeSegmentData() {

#ifdef OUTPUT_PROGRESS
	cout << "Computing segment features..." << endl;
#endif

	int numSourceSegments = (int)mpData->mSourceSegments.size();
	int numTargetSegments = (int)mpData->mTargetSegments.size();

	mSourceSegmentData.resize(numSourceSegments);
	mTargetSegmentData.resize(numTargetSegments);

	// sample points on segments

	if (true) {
#ifdef OUTPUT_PROGRESS
		cout << "sampling points on segments" << endl;
#endif
		for (int srcID = 0; srcID < numSourceSegments; srcID++) {
			if (!extractSegmentSamples(
				mpData->mSourceMesh,
				mpData->mSourceSamples,
				mpData->mSourceSegments[srcID],
				mSourceSegmentData[srcID])) return false;
		}
		for (int tgtID = 0; tgtID < numTargetSegments; tgtID++) {
			if (!extractSegmentSamples(
				mpData->mTargetMesh,
				mpData->mTargetSamples,
				mpData->mTargetSegments[tgtID],
				mTargetSegmentData[tgtID])) return false;
		}
	}

	// sample points on curves

	if (true) {
#ifdef OUTPUT_PROGRESS
		cout << "sampling points on curves" << endl;
#endif
		auto &sourceCurves = mpData->mSourceFeatures.mCurve.mPointClouds;
		auto &targetCurves = mpData->mTargetFeatures.mCurve.mPointClouds;

		int numCurveTypes = 4; // HACK: should check implementation for any update
		for (int srcID = 0; srcID < numSourceSegments; srcID++) {
			mSourceSegmentData[srcID].mCurvePosition.resize(numCurveTypes);
			mSourceSegmentData[srcID].mCurveNormal.resize(numCurveTypes);
		}
		for (int tgtID = 0; tgtID < numTargetSegments; tgtID++) {
			mTargetSegmentData[tgtID].mCurvePosition.resize(numCurveTypes);
			mTargetSegmentData[tgtID].mCurveNormal.resize(numCurveTypes);
		}

		for (int typeID = 0; typeID < numCurveTypes; typeID++) {

			if (sourceCurves.empty() || sourceCurves[typeID].amount == 0) {
				// no source curves extracted
				for (int srcID = 0; srcID < numSourceSegments; srcID++) {
					mSourceSegmentData[srcID].mCurvePosition[typeID] = mSourceSegmentData[srcID].mMatchingPosition.leftCols(1);
					mSourceSegmentData[srcID].mCurveNormal[typeID] = mSourceSegmentData[srcID].mMatchingNormal.leftCols(1);
				}
			} else {
				Eigen::Matrix3Xd sourceCurveMatP, sourceCurveMatN;
				if (!SampleUtil::buildMatrices(sourceCurves[typeID], sourceCurveMatP, sourceCurveMatN)) return false;

				SKDTree sourceTree;
				SKDTreeData sourceTreeData;
				if (!SampleUtil::buildKdTree(sourceCurveMatP, sourceTree, sourceTreeData)) return false;

				for (int srcID = 0; srcID < numSourceSegments; srcID++) {
					Eigen::VectorXi nnIndices;
					if (!SampleUtil::findNearestNeighbors(sourceTree, mSourceSegmentData[srcID].mSegmentPosition, nnIndices)) return false;
					set<int> nnIndexSet(nnIndices.data(), nnIndices.data() + nnIndices.size());
					if (nnIndexSet.size() <= 1) {
						mSourceSegmentData[srcID].mCurvePosition[typeID] = mSourceSegmentData[srcID].mMatchingPosition.leftCols(1);
						mSourceSegmentData[srcID].mCurveNormal[typeID] = mSourceSegmentData[srcID].mMatchingNormal.leftCols(1);
						continue;
					}
					vector<int> nnIndexList(nnIndexSet.begin(), nnIndexSet.end());
					nnIndices = Eigen::Map<Eigen::VectorXi>(nnIndexList.data(), nnIndexList.size());
					if (!SampleUtil::sliceMatrices(sourceCurveMatP, nnIndices, mSourceSegmentData[srcID].mCurvePosition[typeID])) return false;
					if (!SampleUtil::sliceMatrices(sourceCurveMatN, nnIndices, mSourceSegmentData[srcID].mCurveNormal[typeID])) return false;
				}
			}


			if (targetCurves.empty() || targetCurves[typeID].amount == 0) {
				// no target curves extracted
				for (int tgtID = 0; tgtID < numTargetSegments; tgtID++) {
					mTargetSegmentData[tgtID].mCurvePosition[typeID] = mTargetSegmentData[tgtID].mMatchingPosition.leftCols(1);
					mTargetSegmentData[tgtID].mCurveNormal[typeID] = mTargetSegmentData[tgtID].mMatchingNormal.leftCols(1);
				}
			} else {
				Eigen::Matrix3Xd targetCurveMatP, targetCurveMatN;
				if (!SampleUtil::buildMatrices(targetCurves[typeID], targetCurveMatP, targetCurveMatN)) return false;

				SKDTree targetTree;
				SKDTreeData targetTreeData;
				if (!SampleUtil::buildKdTree(targetCurveMatP, targetTree, targetTreeData)) return false;

				for (int tgtID = 0; tgtID < numTargetSegments; tgtID++) {
					Eigen::VectorXi nnIndices;
					if (!SampleUtil::findNearestNeighbors(targetTree, mTargetSegmentData[tgtID].mSegmentPosition, nnIndices)) return false;
					set<int> nnIndexSet(nnIndices.data(), nnIndices.data() + nnIndices.size());
					if (nnIndexSet.size() <= 1) {
						mTargetSegmentData[tgtID].mCurvePosition[typeID] = mTargetSegmentData[tgtID].mMatchingPosition.leftCols(1);
						mTargetSegmentData[tgtID].mCurveNormal[typeID] = mTargetSegmentData[tgtID].mMatchingNormal.leftCols(1);
						continue;
					}
					vector<int> nnIndexList(nnIndexSet.begin(), nnIndexSet.end());
					nnIndices = Eigen::Map<Eigen::VectorXi>(nnIndexList.data(), nnIndexList.size());
					if (!SampleUtil::sliceMatrices(targetCurveMatP, nnIndices, mTargetSegmentData[tgtID].mCurvePosition[typeID])) return false;
					if (!SampleUtil::sliceMatrices(targetCurveMatN, nnIndices, mTargetSegmentData[tgtID].mCurveNormal[typeID])) return false;
				}
			}
		}
	}

	// compute features on segments

	if (true) {
#ifdef OUTPUT_PROGRESS
		cout << "extracting features on segments" << endl;
#endif
		for (int srcID = 0; srcID < numSourceSegments; srcID++) {
			if (!extractSegmentFeatures(
				mpData->mSourceMesh,
				mpData->mSourceFeatures,
				mpData->mSourceSegments[srcID],
				srcID, mSourceSegmentData[srcID])) return false;
		}
		for (int tgtID = 0; tgtID < numTargetSegments; tgtID++) {
			if (!extractSegmentFeatures(
				mpData->mTargetMesh,
				mpData->mTargetFeatures,
				mpData->mTargetSegments[tgtID],
				tgtID, mTargetSegmentData[tgtID])) return false;
		}
	}

	return true;
}

bool SimilarityDistance::computeElementDistance() {

#ifdef OUTPUT_PROGRESS
	cout << "Computing element distance..." << endl;
#endif

	int numSourceSegments = (int)mpData->mSourceSegments.size();
	int numTargetSegments = (int)mpData->mTargetSegments.size();

	int numTotalElementPairs = numSourceSegments * numTargetSegments;
	
	vector<vector<double>> allPairDistances(numTotalElementPairs);
	vector<double> allPairEvaluatedDistances(numTotalElementPairs);
	vector<bool> allPairElementFlags(numTotalElementPairs, false);

	// compute all pair element distance

	double elementDistanceFilteringThreshold = StyleSynthesisConfig::mStyle_ElementDistanceFilteringThreshold;

	int pairCount = 0;
	int filterCount = 0;
//#pragma omp parallel for shared(pairCount)
	for (int pairID = 0; pairID < numTotalElementPairs; pairID++) {		

#ifdef OUTPUT_PROGRESS
#pragma omp atomic
		pairCount++;
		if (pairCount % max(1, numTotalElementPairs / 100) == 0) {
#pragma omp critical
			cout << "\rComputing pair " << pairCount << " / " << numTotalElementPairs << "        ";
#endif
		}
		
		int sourceID = pairID / numTargetSegments;
		int targetID = pairID % numTargetSegments;
		
		TSegmentData &sourceData = mSourceSegmentData[sourceID];
		TSegmentData &targetData = mTargetSegmentData[targetID];

		if (!checkSegmentPairs(sourceData, targetData)) {
#ifdef OUTPUT_PROGRESS
#pragma omp atomic
			filterCount++;
#endif
			continue;
		}

		Eigen::Affine3d surfaceTransformationST, surfaceTransformationTS;
		surfaceTransformationST.setIdentity();
		surfaceTransformationTS.setIdentity();

		// surface ICP distance
		vector<double> surfaceDistance(0);
		if (true) {
			//MatchSimpleICP::visualize("before-ICP-st.ply", sourceData.mMatchingPosition, targetData.mMatchingPosition, surfaceTransformationST);
			//MatchSimpleICP::visualize("before-ICP-ts.ply", targetData.mMatchingPosition, sourceData.mMatchingPosition, surfaceTransformationTS);
			if (!MatchSimpleICP::runAnyShape(20, // UNDONE: param matching iteration
				sourceData.mMatchingPosition, sourceData.mMatchingNormal,
				targetData.mMatchingPosition, targetData.mMatchingNormal,
				surfaceTransformationST)) error("matching surface");
			if (!MatchSimpleICP::runAnyShape(20,
				targetData.mMatchingPosition, targetData.mMatchingNormal,
				sourceData.mMatchingPosition, sourceData.mMatchingNormal,
				surfaceTransformationTS)) error("matching surface");
			//MatchSimpleICP::visualize("after-ICP-st.ply", sourceData.mMatchingPosition, targetData.mMatchingPosition, surfaceTransformationST);
			//MatchSimpleICP::visualize("after-ICP-ts.ply", targetData.mMatchingPosition, sourceData.mMatchingPosition, surfaceTransformationTS);
			Eigen::Matrix3Xd matXSP = surfaceTransformationST * sourceData.mSegmentPosition;
			Eigen::Matrix3Xd matXSN = surfaceTransformationST.rotation() * sourceData.mSegmentNormal;
			Eigen::Matrix3Xd matXTP = surfaceTransformationTS * targetData.mSegmentPosition;
			Eigen::Matrix3Xd matXTN = surfaceTransformationTS.rotation() * targetData.mSegmentNormal;
			double distP1, distN1, distP2, distN2;
			if (!MatchSimpleICP::distance(matXSP, matXSN,
				targetData.mSegmentPosition, targetData.mSegmentNormal,
				distP1, distN1)) error("surface ICP distance");
			if (!MatchSimpleICP::distance(matXTP, matXTN,
				sourceData.mSegmentPosition, sourceData.mSegmentNormal,
				distP2, distN2)) error("surface ICP distance");
			distP1 = distP1 / mpData->mTargetSamples.radius; // make distance invariant to sampling density
			distP2 = distP2 / mpData->mSourceSamples.radius;

			double distP = (distP1 + distP2) / 2; // average from both direction
			double distN = (distN1 + distN2) / 2;
			surfaceDistance.push_back(distP);
			surfaceDistance.push_back(distN);
		}

		// curve ICP distance
		vector<double> curveDistance(0);
		if (true) {
			int numCurveType = (int)sourceData.mCurvePosition.size();

			for (int typeID = 0; typeID < numCurveType; typeID++) {

				Eigen::Matrix3Xd matXSP = surfaceTransformationST * sourceData.mCurvePosition[typeID];
				Eigen::Matrix3Xd matXSN = surfaceTransformationST.rotation() * sourceData.mCurveNormal[typeID];
				Eigen::Matrix3Xd matXTP = surfaceTransformationTS * targetData.mCurvePosition[typeID];
				Eigen::Matrix3Xd matXTN = surfaceTransformationTS.rotation() * targetData.mCurveNormal[typeID];
				double distP1, distN1, distP2, distN2;
				if (!MatchSimpleICP::distance(matXSP, matXSN,
					targetData.mCurvePosition[typeID], targetData.mCurveNormal[typeID],
					distP1, distN1)) error("curve ICP distance");
				if (!MatchSimpleICP::distance(matXTP, matXTN,
					sourceData.mCurvePosition[typeID], sourceData.mCurveNormal[typeID],
					distP2, distN2)) error("curve ICP distance");
				distP1 = distP1 / mpData->mTargetSamples.radius; // make distance invariant to sampling density
				distP2 = distP2 / mpData->mSourceSamples.radius;

				double distP = (distP1 + distP2) / 2; // average from both direction
				double distN = (distN1 + distN2) / 2;
				curveDistance.push_back(distP);
				curveDistance.push_back(distN);
			}
		}

		// scale distance
		vector<double> scaleDistance;
		if (true) {
			Eigen::JacobiSVD<Eigen::Matrix3d> svdST(surfaceTransformationST.linear(), Eigen::ComputeFullU | Eigen::ComputeFullV);
			Eigen::JacobiSVD<Eigen::Matrix3d> svdTS(surfaceTransformationTS.linear(), Eigen::ComputeFullU | Eigen::ComputeFullV);
			Eigen::Vector3d scaleVecST = svdST.singularValues().cwiseAbs();
			Eigen::Vector3d scaleVecTS = svdTS.singularValues().cwiseAbs();
			for (int dim = 0; dim<3; dim++) {
				double distST = scaleVecST[dim] > 1.0 ? scaleVecST[dim] : (1.0 / max(0.01, scaleVecST[dim]));
				double distTS = scaleVecTS[dim] > 1.0 ? scaleVecTS[dim] : (1.0 / max(0.01, scaleVecTS[dim]));
				double dist = (distST + distTS) / 2; // average from both direction
				scaleDistance.push_back(dist);
			}
		}

		// shape distributions distance
		vector<double> sdDistance;
		if (true) {
			if (!FeatureShapeDistributions::compareFeatures(
				sourceData.mFeatures.mSD, targetData.mFeatures.mSD,
				sdDistance)) error("SD distance");
		}

		// curvature distances
		vector<double> curvatureDistance;
		if (true) {
			if (!FeatureSampleCurvature::compareFeatures(
				sourceData.mFeatures.mCurvature, targetData.mFeatures.mCurvature,
				curvatureDistance)) error("curvature distance");
		}


		// SDF distance
		vector<double> sdfDistance;
		if (true) {
			if (!FeatureShapeDiameterFunctions::compareFeatures(
				sourceData.mFeatures.mSDF, targetData.mFeatures.mSDF,
				sdfDistance)) error("SDF distance");
		}

		// LFD distance
		vector<double> lfdDistance;
		if (true) {
			if (!FeatureLightFieldDescriptors::compareFeatures(
				sourceData.mFeatures.mLFD, targetData.mFeatures.mLFD,
				lfdDistance)) error("LFD distance");
		}

		// combine all distances

		vector<double> &allDistances = allPairDistances[pairID];
		allDistances.clear();
		allDistances.insert(allDistances.end(), surfaceDistance.begin(), surfaceDistance.end());
		allDistances.insert(allDistances.end(), curveDistance.begin(), curveDistance.end());
		allDistances.insert(allDistances.end(), scaleDistance.begin(), scaleDistance.end());
		allDistances.insert(allDistances.end(), sdDistance.begin(), sdDistance.end());
		allDistances.insert(allDistances.end(), curvatureDistance.begin(), curvatureDistance.end());
		allDistances.insert(allDistances.end(), sdfDistance.begin(), sdfDistance.end());
		allDistances.insert(allDistances.end(), lfdDistance.begin(), lfdDistance.end());

		// filter element pairs

		Eigen::VectorXd allDistanceVec = Eigen::Map<Eigen::VectorXd>(allDistances.data(), allDistances.size());
		double evaluatedDistance = SimilarityMetric::evaluateDistance(allDistanceVec);
		allPairEvaluatedDistances[pairID] = evaluatedDistance;

		if (evaluatedDistance < elementDistanceFilteringThreshold) {
			allPairElementFlags[pairID] = true;
			
		}
	}
#ifdef OUTPUT_PROGRESS
	cout << endl;
	cout << "Filtered " << filterCount << " pairs" << endl;
#endif

	// keep only the best pair for each segment

	if (false) {
		vector<int> sourceElementMap(numSourceSegments, -1);
		vector<int> targetElementMap(numTargetSegments, -1);
		for (int sourceID = 0; sourceID < numSourceSegments; sourceID++) {
			double minDist = DBL_MAX;
			for (int targetID = 0; targetID < numTargetSegments; targetID++) {
				int pairID = sourceID * numTargetSegments + targetID;
				if (!allPairElementFlags[pairID]) continue;
				double pairDist = allPairEvaluatedDistances[pairID];
				if (pairDist < minDist) {
					minDist = pairDist;
					sourceElementMap[sourceID] = targetID;
				}
			}
		}
		for (int targetID = 0; targetID < numTargetSegments; targetID++) {
			double minDist = DBL_MAX;
			for (int sourceID = 0; sourceID < numSourceSegments; sourceID++) {
				int pairID = sourceID * numTargetSegments + targetID;
				if (!allPairElementFlags[pairID]) continue;
				double pairDist = allPairEvaluatedDistances[pairID];
				if (pairDist < minDist) {
					minDist = pairDist;
					targetElementMap[targetID] = sourceID;
				}
			}
		}
		allPairElementFlags.assign(numTotalElementPairs, false);
		for (int sourceID = 0; sourceID < numSourceSegments; sourceID++) {
			int targetID = sourceElementMap[sourceID];
			if (targetID < 0) continue;
			int pairID = sourceID * numTargetSegments + targetID;
			allPairElementFlags[pairID] = true;
		}
		for (int targetID = 0; targetID < numTargetSegments; targetID++) {
			int sourceID = targetElementMap[targetID];
			if (sourceID < 0) continue;
			int pairID = sourceID * numTargetSegments + targetID;
			allPairElementFlags[pairID] = true;
		}
	}

	// gather element data

	int numElementPairs = 0;
	int numDistanceDimensions = 0;
	for (int pairID = 0; pairID < numTotalElementPairs; pairID++) {
		if (allPairElementFlags[pairID]) {
			numElementPairs++;
			if (!numDistanceDimensions) numDistanceDimensions = (int)allPairDistances[pairID].size();
		}
	}
#ifdef OUTPUT_PROGRESS
	cout << "Extracted " << numElementPairs << " out of " << numTotalElementPairs << " element pairs" << endl;
#endif

	mDistanceData.mElementIndices.resize(numElementPairs);
	mDistanceData.mElementSourcePoints.resize(numElementPairs);
	mDistanceData.mElementTargetPoints.resize(numElementPairs);
	mDistanceData.mElementDistance.resize(numElementPairs, numDistanceDimensions);

	vector<bool> sourcePointFlags(mpData->mSourceSamples.amount, false);
	vector<bool> targetPointFlags(mpData->mTargetSamples.amount, false);

	int elementID = 0;
	for (int pairID = 0; pairID < numTotalElementPairs; pairID++) {
		if (!allPairElementFlags[pairID]) continue;
		int sourceID = pairID / numTargetSegments;
		int targetID = pairID % numTargetSegments;

		mDistanceData.mElementIndices[elementID] = vec2i(sourceID, targetID);
		mDistanceData.mElementSourcePoints[elementID] = mSourceSegmentData[sourceID].mSegmentPointIndices;
		mDistanceData.mElementTargetPoints[elementID] = mTargetSegmentData[targetID].mSegmentPointIndices;
		for (int pointID : mDistanceData.mElementSourcePoints[elementID]) sourcePointFlags[pointID] = true;
		for (int pointID : mDistanceData.mElementTargetPoints[elementID]) targetPointFlags[pointID] = true;

		for (int dim = 0; dim < numDistanceDimensions; dim++) {
			mDistanceData.mElementDistance(elementID, dim) = allPairDistances[pairID][dim];
		}

		elementID++;
	}

	// gather unmatched points

	mDistanceData.mUnmatchSourcePoints.clear();
	mDistanceData.mUnmatchTargetPoints.clear();
	for (int pointID = 0; pointID < mpData->mSourceSamples.amount; pointID++) {
		if (!sourcePointFlags[pointID]) mDistanceData.mUnmatchSourcePoints.push_back(pointID);
	}
	for (int pointID = 0; pointID < mpData->mTargetSamples.amount; pointID++) {
		if (!targetPointFlags[pointID]) mDistanceData.mUnmatchTargetPoints.push_back(pointID);
	}

	return true;
}

bool SimilarityDistance::extractSegmentSamples(
	TTriangleMesh &inMesh,
	TSampleSet &inSample,
	vector<int> &inSegment,
	TSegmentData &outData)
{

	// gather sample points on segment

	int numFaces = (int)inMesh.indices.size();
	vector<bool> faceFlags(numFaces, false);
	for (int faceID : inSegment) faceFlags[faceID] = true;

	outData.mSegmentPointIndices.clear();
	outData.mSegmentPointSet.positions.clear();
	outData.mSegmentPointSet.normals.clear();
	for (int sampleID = 0; sampleID < inSample.amount; sampleID++) {
		if (faceFlags[inSample.indices[sampleID]]) {
			outData.mSegmentPointIndices.push_back(sampleID);
			outData.mSegmentPointSet.positions.push_back(inSample.positions[sampleID]);
			outData.mSegmentPointSet.normals.push_back(inSample.normals[sampleID]);
		}
	}
	outData.mSegmentPointSet.amount = (int)outData.mSegmentPointSet.positions.size();

	if (!SampleUtil::buildMatrices(outData.mSegmentPointSet, outData.mSegmentPosition, outData.mSegmentNormal)) return false;

	int numSamples = outData.mSegmentPointSet.amount;

	if (numSamples >= 20 && numSamples <= 1024) {
		outData.mMatchingPosition = outData.mSegmentPosition;
		outData.mMatchingNormal = outData.mSegmentNormal;
	} else {

		// generate constant amount of sample points for matching

		TTriangleMesh segmentMesh;
		segmentMesh.positions = inMesh.positions;
		segmentMesh.normals = inMesh.normals;
		segmentMesh.amount = inMesh.amount;
		segmentMesh.indices.clear();
		for (int faceID : inSegment) segmentMesh.indices.push_back(inMesh.indices[faceID]);
		if (!MeshUtil::cleanUp(segmentMesh)) return false;

		numSamples = cml::clamp(numSamples, 20, 1024); // UNDONE: param reduced amount of samples
		TSampleSet reducedSamples;
		SampleSimplePoissonDisk sspd(&segmentMesh);
		if (!sspd.runSampling(numSamples)) return false;
		if (!sspd.exportSample(reducedSamples)) return false;
		if (!SampleUtil::buildMatrices(reducedSamples, outData.mMatchingPosition, outData.mMatchingNormal)) return false;
	}

	return true;
}

bool SimilarityDistance::extractSegmentFeatures(
	TTriangleMesh &inMesh,
	FeatureAsset &inFeatures,
	vector<int> &inSegment,
	int &inSegmentID,
	TSegmentData &outData)
{
	// sample feature values on mesh

	if (true) {
		int numPoints = outData.mSegmentPointSet.amount;
		outData.mFeatures.mCurvature.resize(numPoints);
		outData.mFeatures.mSDF.resize(numPoints);
		for (int pointID = 0; pointID < numPoints; pointID++) {
			int meshSampleID = outData.mSegmentPointIndices[pointID];
			outData.mFeatures.mCurvature[pointID] = inFeatures.mCurvature[meshSampleID];
			outData.mFeatures.mSDF[pointID] = inFeatures.mSDF[meshSampleID];
		}
	}

	// extract part features

	if (true) {
		outData.mFeatures.mLFD = inFeatures.mPartLFD[inSegmentID];
		outData.mFeatures.mSD = inFeatures.mPartSD[inSegmentID];
	}

	// compute intrinsic properties

	if (true) {

		Eigen::Matrix3Xd mat = outData.mMatchingPosition;

		// volume
		Eigen::AlignedBox3d boundingBox(mat.rowwise().minCoeff(), mat.rowwise().maxCoeff());
		outData.mVolume = boundingBox.volume();

		// principal variance
		
		Eigen::Matrix3Xd centeredPoints = mat.colwise() - mat.rowwise().mean();
		Eigen::JacobiSVD< Eigen::Matrix3Xd > svd(centeredPoints);
		outData.mVariance = svd.singularValues();
	}

	return true;
}

bool SimilarityDistance::checkSegmentPairs(TSegmentData &source, TSegmentData &target) {

	double varianceFactor = 5.0; // UNDONE: param element pair filtering threshold
	double volumeFactor = 10.0;

	if (source.mSegmentPointSet.amount < 20) return false;
	if (target.mSegmentPointSet.amount < 20) return false;

	if (source.mVolume > target.mVolume * volumeFactor ||
		source.mVolume < target.mVolume / volumeFactor) return false;

	for (int dim = 0; dim < 3; dim++) {
		if (source.mVariance[dim] > target.mVariance[dim] * varianceFactor ||
			source.mVariance[dim] < target.mVariance[dim] / varianceFactor) return false;
	}

	return true;
}

bool SimilarityDistance::output(string folderName) {

	if (!mDistanceData.saveData(folderName)) return false;

	return true;
}

bool SimilarityDistance::visualize(string fileName) {

	if (!visualizeMatchingElements(fileName, *mpData, mDistanceData)) return false;

	return true;
}

bool SimilarityDistance::visualizeMatchingElements(
	string fileName,
	SimilarityData &data,
	SimilarityDistanceData &distanceData)
{

	int numElementPairs = (int)distanceData.mElementSourcePoints.size();

	auto &srcMesh = data.mSourceMesh;
	auto &tgtMesh = data.mTargetMesh;

	// compute spacing

	vec3 srcBBMin, srcBBMax;
	vec3 tgtBBMin, tgtBBMax;
	if (!MeshUtil::computeAABB(srcMesh, srcBBMin, srcBBMax)) return false;
	if (!MeshUtil::computeAABB(tgtMesh, tgtBBMin, tgtBBMax)) return false;

	float hSpace = ((srcBBMax[0] - srcBBMin[0]) + (tgtBBMax[0] - tgtBBMin[0])) * 1.5f;
	float vSpace = max((srcBBMax[1] - srcBBMin[1]), (tgtBBMax[1] - tgtBBMin[1])) * 1.5f;
	float stSpace = (srcBBMax[0] - tgtBBMin[0]) * 1.2f;

	// sort elements by style distance

	Eigen::VectorXd elementDistance;
	Eigen::VectorXd elementSaliency;
	double unmatchSaliency;
	if (!SimilarityMetric::evaluateSimilarity(data, distanceData, elementDistance, elementSaliency, unmatchSaliency)) return false;

	Eigen::VectorXd elementScoreVec = elementDistance.array() * elementSaliency.array();
	//Eigen::VectorXd elementScoreVec = elementSaliency.array();
	//Eigen::VectorXd elementScoreVec = elementDistance.array() / elementSaliency.array();
	double distanceThreshold = StyleSynthesisConfig::mStyle_ElementDistanceFilteringThreshold;
	elementScoreVec = (elementDistance.array() > distanceThreshold).select(0.0, elementScoreVec); // HACK: prune weird elements
	vector<double> elementScore(elementScoreVec.data(), elementScoreVec.data() + elementScoreVec.size());
	vector<int> elementOrder(numElementPairs);
	for (int k = 0; k < numElementPairs; k++) elementOrder[k] = k;
	sort(elementOrder.begin(), elementOrder.end(),
		[&elementScore](int lhs, int rhs){ return elementScore[lhs] > elementScore[rhs]; });

	// visualize elements

	PlyExporter pe;
	numElementPairs = min(numElementPairs, 30); // only visualize first 30 pairs...
	for (int orderID = 0; orderID < numElementPairs; orderID++) {
		int elementID = elementOrder[orderID];
		//cout << elementDistance[elementID] << ", " << elementSaliency[elementID] << ", " << elementScore[elementID] << endl;
		int hID = orderID % 5;
		int vID = orderID / 5;
		vec3 srcOffset(hID*hSpace, -vID*vSpace, 0.0f);
		vec3 tgtOffset = srcOffset + vec3(stSpace, 0.0f, 0.0f);

		vector<vec3i> srcColors(srcMesh.indices.size(), vec3i(127, 127, 127));
		vector<vec3i> tgtColors(tgtMesh.indices.size(), vec3i(127, 127, 127));
		int srcSegID = distanceData.mElementIndices[elementID][0];
		int tgtSegID = distanceData.mElementIndices[elementID][1];
		for (int faceID : data.mSourceSegments[srcSegID]) srcColors[faceID] = vec3i(255, 0, 0);
		for (int faceID : data.mTargetSegments[tgtSegID]) tgtColors[faceID] = vec3i(255, 0, 0);

		if (!pe.addMesh(&srcMesh.positions, &srcMesh.normals, &srcMesh.indices, &srcColors, srcOffset)) return false;
		if (!pe.addMesh(&tgtMesh.positions, &tgtMesh.normals, &tgtMesh.indices, &tgtColors, tgtOffset)) return false;
	}
	if (!pe.output(fileName)) return false;

	return true;
}

bool SimilarityDistance::evaluateFaceContribution(
	double &outShapeDistance,
	vector<double> &outSourceContribution,
	vector<double> &outTargetContribution,
	vector<double> &outSourceSaliency,
	vector<double> &outTargetSaliency)
{

	int numSourceFaces = (int)mpData->mSourceMesh.indices.size();
	int numTargetFaces = (int)mpData->mTargetMesh.indices.size();
	outSourceContribution.assign(numSourceFaces, 0);
	outTargetContribution.assign(numTargetFaces, 0);
	outSourceSaliency.assign(numSourceFaces, 0);
	outTargetSaliency.assign(numTargetFaces, 0);

	// evaluate distance terms

	Eigen::VectorXd elementPairDistance;
	Eigen::VectorXd elementPairSaliency;
	double unmatchSaliency;
	if (!SimilarityMetric::evaluateSimilarity(*mpData, mDistanceData, elementPairDistance, elementPairSaliency, unmatchSaliency)) return false;

	double shapeDistance = elementPairSaliency.dot(elementPairDistance) + unmatchSaliency * SimilarityMetric::mWeightsPrevalence;
	Eigen::VectorXd elementPairContribution = elementPairDistance.array() * elementPairSaliency.array() / shapeDistance;
	outShapeDistance = shapeDistance;

	// account for matched region contribution

	vector<bool> sourceFaceMatchedFlags(numSourceFaces, false);
	vector<bool> targetFaceMatchedFlags(numTargetFaces, false);

	int numElementPairs = (int)elementPairDistance.size();
	for (int pairID = 0; pairID < numElementPairs; pairID++) {
		vec2i elementPair = mDistanceData.mElementIndices[pairID];
		int sourceID = elementPair[0];
		int targetID = elementPair[1];
		double contribution = elementPairContribution[pairID];
		auto &srcPatch = mpData->mSourceSegments[sourceID];
		auto &tgtPatch = mpData->mTargetSegments[targetID];
		double srcContribution = contribution / (double)srcPatch.size();
		double tgtContribution = contribution / (double)tgtPatch.size();
		for (int faceID : srcPatch) {
			outSourceContribution[faceID] += srcContribution;
			sourceFaceMatchedFlags[faceID] = true;
		}
		for (int faceID : tgtPatch) {
			outTargetContribution[faceID] += tgtContribution;
			targetFaceMatchedFlags[faceID] = true;
		}
	}

	// account for unmatched region contribution

	Eigen::VectorXd sourcePointSaliency; // this should already be calculated when evaluating similarity...
	Eigen::VectorXd targetPointSaliency; // but it should be pretty fast to compute this...
	if (!SimilarityMetric::evaluatePointSaliency(*mpData, mDistanceData, sourcePointSaliency, targetPointSaliency)) return false;
	Eigen::VectorXd srcPointUnmatchContribution = sourcePointSaliency * SimilarityMetric::mWeightsPrevalence / shapeDistance;
	Eigen::VectorXd tgtPointUnmatchContribution = targetPointSaliency * SimilarityMetric::mWeightsPrevalence / shapeDistance;

	for (int srcPointID = 0; srcPointID < mpData->mSourceSamples.amount; srcPointID++) {
		int faceID = mpData->mSourceSamples.indices[srcPointID];
		if (sourceFaceMatchedFlags[faceID]) continue;
		outSourceContribution[faceID] += srcPointUnmatchContribution[srcPointID];
	}
	for (int tgtPointID = 0; tgtPointID < mpData->mTargetSamples.amount; tgtPointID++) {
		int faceID = mpData->mTargetSamples.indices[tgtPointID];
		if (targetFaceMatchedFlags[faceID]) continue;
		outTargetContribution[faceID] += tgtPointUnmatchContribution[tgtPointID];
	}

	// compute face saliency

	for (int srcPointID = 0; srcPointID < mpData->mSourceSamples.amount; srcPointID++) {
		int faceID = mpData->mSourceSamples.indices[srcPointID];
		outSourceSaliency[faceID] += sourcePointSaliency[srcPointID];
	}
	for (int tgtPointID = 0; tgtPointID < mpData->mTargetSamples.amount; tgtPointID++) {
		int faceID = mpData->mTargetSamples.indices[tgtPointID];
		outTargetSaliency[faceID] += targetPointSaliency[tgtPointID];
	}

	if (true) {
		// sanity check
		double sumSourceContrib = 0;
		double sumTargetContrib = 0;
		for (int srcID = 0; srcID < numSourceFaces; srcID++) sumSourceContrib += outSourceContribution[srcID];
		for (int tgtID = 0; tgtID < numTargetFaces; tgtID++) sumTargetContrib += outTargetContribution[tgtID];
		cout << "Total source contribution = " << sumSourceContrib << endl;
		cout << "Total target contribution = " << sumTargetContrib << endl;
		cout << "Total contribution (should be 2) = " << (sumSourceContrib + sumTargetContrib) << endl;

		double sumSourceSaliency = 0;
		double sumTargetSaliency = 0;
		for (int srcID = 0; srcID < numSourceFaces; srcID++) sumSourceSaliency += outSourceSaliency[srcID];
		for (int tgtID = 0; tgtID < numTargetFaces; tgtID++) sumTargetSaliency += outTargetSaliency[tgtID];
		cout << "Total source saliency = " << sumSourceSaliency << ", " << sourcePointSaliency.sum() << endl;
		cout << "Total target saliency = " << sumTargetSaliency << ", " << targetPointSaliency.sum() << endl;
	}

	return true;
}