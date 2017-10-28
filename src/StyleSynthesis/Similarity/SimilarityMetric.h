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

#include <Eigen/Dense>

#include "Similarity/SimilarityData.h"
#include "Similarity/SimilarityDistanceData.h"

using namespace std;

namespace StyleSynthesis {

	class SimilarityMetric {

	private:

		SimilarityMetric() {}
		~SimilarityMetric() {}

	public:

		static bool initMetric(string weightsFolder, string weigtsAffix = "");

		static bool evaluateSimilarity(
			SimilarityData &inData,
			SimilarityDistanceData &inDistanceData,
			double &outDistance);

		static bool evaluateSimilarity(
			SimilarityData &inData,
			SimilarityDistanceData &inDistanceData,
			Eigen::VectorXd &outElementDistance,
			Eigen::VectorXd &outElementSaliency,
			double &outUnmatchSaliency);

		static bool evaluatePointSaliency(
			SimilarityData &inData,
			SimilarityDistanceData &inDistanceData,
			Eigen::VectorXd &outSourcePointSaliency,
			Eigen::VectorXd &outTargetPointSaliency);

		inline static double evaluateDistance(Eigen::VectorXd &distanceVec) {
			Eigen::VectorXd normalizedVec = distanceVec.cwiseQuotient(SimilarityMetric::mScaleDistance);
			normalizedVec = normalizedVec.cwiseMin(1.0);
			//cout << endl << "Style distance = " << (normalizedVec.array()*SimilarityMetric::mWeightsDistance.array()).transpose() << endl;
			//cout << "Total = " << normalizedVec.dot(SimilarityMetric::mWeightsDistance) << endl;
			return normalizedVec.dot(SimilarityMetric::mWeightsDistance);
		}

		inline static double evaluateSaliency(Eigen::VectorXd &saliencyVec) {
			Eigen::VectorXd normalizedVec = saliencyVec.cwiseQuotient(SimilarityMetric::mScaleSaliency);
			normalizedVec = normalizedVec.cwiseMax(-1.0);
			normalizedVec = normalizedVec.cwiseMin(1.0);
			Eigen::VectorXd homoVec(normalizedVec.size() + 1);
			homoVec << normalizedVec, 1.0;
			double s = homoVec.dot(SimilarityMetric::mWeightsSaliency);
			return 1.0 / (1.0 + exp(-s));
		}

	public:

		static Eigen::VectorXd mWeightsDistance;
		static Eigen::VectorXd mWeightsSaliency;
		static double          mWeightsPrevalence;

		static Eigen::VectorXd mScaleDistance;
		static Eigen::VectorXd mScaleSaliency;
	};
}