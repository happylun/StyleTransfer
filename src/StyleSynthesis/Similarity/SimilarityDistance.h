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

#include "Feature/FeatureAsset.h"

#include "Similarity/SimilarityData.h"
#include "Similarity/SimilarityDistanceData.h"

using namespace std;

namespace StyleSynthesis {

	class SimilarityDistance {

	public:

		SimilarityDistance(SimilarityData *data);
		~SimilarityDistance();

	public:

		bool process();
		bool output(string folderName);
		bool visualize(string fileName);

		static bool visualizeMatchingElements(
			string fileName,
			SimilarityData &data,
			SimilarityDistanceData &distanceData);

		bool evaluateFaceContribution(
			double &outShapeDistance,
			vector<double> &outSourceContribution,
			vector<double> &outTargetContribution,
			vector<double> &outSourceSaliency,
			vector<double> &outTargetSaliency); // this should be done after having learned weights

	private:

		bool computeSegmentData();
		bool computeElementDistance();

	private:

		inline static bool error(string s) { cout << "Error: " << s << endl; return false; }

	private:

		SimilarityData *mpData;
		SimilarityDistanceData mDistanceData;

	///////////// for internal use /////////////

	private:

		struct TSegmentData {
			vector<int> mSegmentPointIndices;
			TPointSet mSegmentPointSet;
			Eigen::Matrix3Xd mSegmentPosition;
			Eigen::Matrix3Xd mSegmentNormal;
			vector<Eigen::Matrix3Xd> mCurvePosition;
			vector<Eigen::Matrix3Xd> mCurveNormal;
			FeatureAsset mFeatures;
			// reduced sample points -- used for matching
			Eigen::Matrix3Xd mMatchingPosition;
			Eigen::Matrix3Xd mMatchingNormal;
			// intrinsic properties -- used for pruning pairs
			double mVolume;
			Eigen::Vector3d mVariance;
		};

		bool extractSegmentSamples(
			TTriangleMesh &inMesh,
			TSampleSet &inSample,
			vector<int> &inSegment,
			TSegmentData &outData);

		bool extractSegmentFeatures(
			TTriangleMesh &inMesh,
			FeatureAsset &inFeatures,
			vector<int> &inSegment,
			int &inSegmentID,
			TSegmentData &outData);

		bool checkSegmentPairs(
			TSegmentData &source,
			TSegmentData &target);

	private:

		vector<TSegmentData> mSourceSegmentData;
		vector<TSegmentData> mTargetSegmentData;

	};
}