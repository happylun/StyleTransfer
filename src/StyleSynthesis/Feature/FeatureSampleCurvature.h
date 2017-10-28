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

#include "Data/StyleSynthesisTypes.h"

using namespace std;

namespace StyleSynthesis {

	class FeatureSampleCurvature {

		// Curvature

	public:

		struct TCurvature {
			float k1; // max curvature
			float k2; // min curvature
			vec3 d1; // max curvature direction
			vec3 d2; // min curvature direction
			// {n, d1, d2} span an orthonormal right-hand-side coordinate system
		};

	public:

		FeatureSampleCurvature(TSampleSet *samples, vector<TCurvature> *curvatures);
		~FeatureSampleCurvature();
	
	public:

		bool calculate();
		bool visualize(string fileName);

	public:
		
		// 13 metrics
		static bool getMetrics(TCurvature inCurvature, vector<double> &outMetrics);

		// 4 hist X 13 metrics
		static bool compareFeatures(vector<TCurvature> &curvature1, vector<TCurvature> &curvature2, vector<double> &distance);

		static bool saveFeature(string fileName, vector<TCurvature> &curvatures);
		static bool loadFeature(string fileName, vector<TCurvature> &curvatures);

	private:

		bool buildGeodesicGraph();
		bool getPatchNeighbors();
		bool runPatchFitting();

	private:

		TSampleSet *mpSamples;
		vector<TCurvature> *mpCurvatures;

		vector<vector<int>> mGraph; // neighbor ID : # of neighbors : # of samples
		vector<vector<int>> mPatch; // point ID : # of points on patch : # of samples
	};
}