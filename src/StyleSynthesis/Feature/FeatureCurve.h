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

	class FeatureCurve {

		// feature curve (generate point clouds to be used in element distance)

	public:

		FeatureCurve();
		~FeatureCurve();
	
	public:

		bool loadCurves(string curveFolder);
		bool extractPointClouds(double sampleRadius);

		bool visualize(string fileName);
		bool loadFeature(string fileName);
		bool saveFeature(string fileName);

	private:

		vector<vector<vec3>> mCurveRV;        // ridges & valleys
		vector<vector<vec3>> mCurveB;         // boundaries
		vector<vector<vector<vec3>>> mCurveC; // contours from multiple view points

	public:

		vector<TPointSet> mPointClouds;
	};
}