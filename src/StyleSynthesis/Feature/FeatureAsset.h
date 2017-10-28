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

#include <vector>
#include <string>

#include <Eigen/Eigen>

#include "Feature/FeatureSampleCurvature.h"
#include "Feature/FeatureShapeDiameterFunctions.h"
#include "Feature/FeatureLightFieldDescriptors.h"
#include "Feature/FeatureShapeDistributions.h"
#include "Feature/FeatureCurve.h"

#include "Feature/FeatureUtil.h"

using namespace std;

namespace StyleSynthesis {

	class FeatureAsset {

	public:

		FeatureAsset();
		~FeatureAsset();

		friend class SimilarityDistance;
		friend class FeatureSaliency;

	public:

		bool loadAllFeatures(string folderName);

		static bool saveFeature(string fileName, vector<double> &feature);
		static bool loadFeature(string fileName, vector<double> &feature);

		static bool savePartFeature(string fileName, vector<vector<double>> &feature);
		static bool loadPartFeature(string fileName, vector<vector<double>> &feature);

	protected:

		vector<FeatureSampleCurvature::TCurvature> mCurvature;
		vector<double> mSDF;
		vector<double> mLFD;
		vector<double> mSD;
		FeatureCurve mCurve;

		vector<vector<double>> mPartLFD;
		vector<vector<double>> mPartSD;
	};
}