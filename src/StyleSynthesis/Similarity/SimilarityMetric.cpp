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

#include "SimilarityMetric.h"

#include <fstream>

#include "Data/DataUtil.h"

using namespace StyleSynthesis;

Eigen::VectorXd SimilarityMetric::mWeightsDistance;
Eigen::VectorXd SimilarityMetric::mWeightsSaliency;
double          SimilarityMetric::mWeightsPrevalence;

Eigen::VectorXd SimilarityMetric::mScaleDistance;
Eigen::VectorXd SimilarityMetric::mScaleSaliency;

bool SimilarityMetric::initMetric(string weightsFolder, string weigtsAffix) {

	// names

	string weightsDistanceName = weightsFolder + "weights-style-distance" + weigtsAffix + ".txt";
	string weightsSaliencyName = weightsFolder + "weights-style-saliency" + weigtsAffix + ".txt";
	string weightsPrevalenceName = weightsFolder + "weights-style-prevalence" + weigtsAffix + ".txt";

	string scaleDistanceName = weightsFolder + "scale-style-distance" + weigtsAffix + ".txt";
	string scaleSaliencyName = weightsFolder + "scale-style-saliency" + weigtsAffix + ".txt";

	// data

	if (!DataUtil::loadVectorASCII(weightsDistanceName, mWeightsDistance)) return false;
	if (!DataUtil::loadVectorASCII(weightsSaliencyName, mWeightsSaliency)) return false;
	if (!DataUtil::loadVectorASCII(scaleDistanceName, mScaleDistance)) return false;
	if (!DataUtil::loadVectorASCII(scaleSaliencyName, mScaleSaliency)) return false;

	vector<double> valueList;
	if (!DataUtil::loadValueListASCII(weightsPrevalenceName, valueList)) return false;
	mWeightsPrevalence = valueList[0];

	// normalize weights (this should already be done)

	double denom = mWeightsDistance.sum() + mWeightsPrevalence;
	if (denom) {
		mWeightsDistance /= denom;
		mWeightsPrevalence /= denom;
	}

	return true;
}

bool SimilarityMetric::evaluateSimilarity(
	SimilarityData &inData,
	SimilarityDistanceData &inDistanceData,
	double &outDistance)
{

	Eigen::VectorXd elementDistance;
	Eigen::VectorXd elementSaliency;
	double unmatchSaliency;

	if (!evaluateSimilarity(inData, inDistanceData, elementDistance, elementSaliency, unmatchSaliency)) return false;
	outDistance = elementSaliency.dot(elementDistance) + unmatchSaliency * mWeightsPrevalence;
	/*
	cout << elementSaliency.sum() << ", " << elementDistance.sum() << ", " << unmatchSaliency << endl;
	cout << outDistance << endl;
	system("pause");
	*/

	return true;
}

bool SimilarityMetric::evaluateSimilarity(
	SimilarityData &inData,
	SimilarityDistanceData &inDistanceData,
	Eigen::VectorXd &outElementDistance,
	Eigen::VectorXd &outElementSaliency,
	double &outUnmatchSaliency)
{

	// compute point saliency

	Eigen::VectorXd sourcePointSaliency;
	Eigen::VectorXd targetPointSaliency;
	if (!evaluatePointSaliency(inData, inDistanceData, sourcePointSaliency, targetPointSaliency)) return false;

	// compute element saliency and unmatched region saliency

	int numElements = (int)inDistanceData.mElementIndices.size();
	outElementSaliency.resize(numElements);
#pragma omp parallel for
	for (int elementID = 0; elementID < numElements; elementID++) {
		double sourcePatchSaliency = 0;
		double targetPatchSaliency = 0;
		for (int pointID : inDistanceData.mElementSourcePoints[elementID]) {
			sourcePatchSaliency += sourcePointSaliency[pointID];
		}
		for (int pointID : inDistanceData.mElementTargetPoints[elementID]) {
			targetPatchSaliency += targetPointSaliency[pointID];
		}
		outElementSaliency[elementID] = (sourcePatchSaliency + targetPatchSaliency) * 0.5;
		//outElementSaliency[elementID] = min(sourcePatchSaliency, targetPatchSaliency); // HACK: use min saliency
	}

	double sourceUnmatchSaliency = 0;
	double targetUnmatchSaliency = 0;
	for (int pointID : inDistanceData.mUnmatchSourcePoints) {
		sourceUnmatchSaliency += sourcePointSaliency[pointID];
	}
	for (int pointID : inDistanceData.mUnmatchTargetPoints) {
		targetUnmatchSaliency += targetPointSaliency[pointID];
	}
	outUnmatchSaliency = (sourceUnmatchSaliency + targetUnmatchSaliency) * 0.5;

	// compute element distance

	outElementDistance.resize(numElements);
#pragma omp parallel for
	for (int elementID = 0; elementID < numElements; elementID++) {
		Eigen::VectorXd distanceVec = inDistanceData.mElementDistance.row(elementID);
		outElementDistance[elementID] = evaluateDistance(distanceVec);
	}

	return true;
}

bool SimilarityMetric::evaluatePointSaliency(
	SimilarityData &inData,
	SimilarityDistanceData &inDistanceData,
	Eigen::VectorXd &outSourcePointSaliency,
	Eigen::VectorXd &outTargetPointSaliency)
{
	// compute point saliency

	int numSourcePoints = (int)inData.mSourceSaliency.rows();
	int numTargetPoints = (int)inData.mTargetSaliency.rows();
	outSourcePointSaliency.resize(numSourcePoints);
	outTargetPointSaliency.resize(numTargetPoints);

#pragma omp parallel for
	for (int srcPointID = 0; srcPointID < numSourcePoints; srcPointID++) {
		Eigen::VectorXd saliencyVec = inData.mSourceSaliency.row(srcPointID).transpose();
		outSourcePointSaliency[srcPointID] = evaluateSaliency(saliencyVec);
	}
	outSourcePointSaliency /= outSourcePointSaliency.sum();

#pragma omp parallel for
	for (int tgtPointID = 0; tgtPointID < numTargetPoints; tgtPointID++) {
		Eigen::VectorXd saliencyVec = inData.mTargetSaliency.row(tgtPointID).transpose();
		outTargetPointSaliency[tgtPointID] = evaluateSaliency(saliencyVec);
	}
	outTargetPointSaliency /= outTargetPointSaliency.sum();

	// normalize point saliency by element count

	vector<int> sourcePointElementCount(numSourcePoints, 0);
	vector<int> targetPointElementCount(numTargetPoints, 0);
	for (auto &patch : inDistanceData.mElementSourcePoints) {
		for (int pointID : patch) sourcePointElementCount[pointID]++;
	}
	for (auto &patch : inDistanceData.mElementTargetPoints) {
		for (int pointID : patch) targetPointElementCount[pointID]++;
	}
	for (int srcPointID = 0; srcPointID < numSourcePoints; srcPointID++) {
		outSourcePointSaliency[srcPointID] /= (double)max(1, sourcePointElementCount[srcPointID]);
	}
	for (int tgtPointID = 0; tgtPointID < numTargetPoints; tgtPointID++) {
		outTargetPointSaliency[tgtPointID] /= (double)max(1, targetPointElementCount[tgtPointID]);
	}

	return true;
}