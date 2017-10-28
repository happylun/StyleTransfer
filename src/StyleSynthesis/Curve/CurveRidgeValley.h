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
#include <set>

#include "Data/StyleSynthesisTypes.h"
#include "Feature/FeatureMeshCurvature.h"

using namespace std;

namespace StyleSynthesis {

	class CurveRidgeValley {

	public:

		CurveRidgeValley(TTriangleMesh &mesh);
		~CurveRidgeValley();

	public:

		bool extractCurve();
		bool output(vector<vector<vec3>> &curves);
		bool visualize(string fileName);

	private:

		bool computeCurvatures();
		bool extractZeroCrossings();
		bool chainCurves();
		bool snapCurves();

		static vec2i makeKey(int i1, int i2) { return i1 < i2 ? vec2i(i1, i2) : vec2i(i2, i1); }
		static vec2i makeKey(vec2i key) { return makeKey(key[0], key[1]); }

	private:

		TTriangleMesh mMesh;

		vector<FeatureMeshCurvature::TCurvature> mCurvatures; // curvature info : # of vertices
		vector<vec4> mCurvatureDerivatives; // 4 unique elements of the 2x2x2 tensor : # of vertices

		vector<vec3> mCurvePoints;
		vector<vec2i> mCurvePointSource; // zero-crossing edge : # of curve points
		vector<double> mCurvePointAlpha; // alpha value for interpolation : # of curve points
		vector<float> mCurvePointThickness; // thickness = abs(KMax)
		vector<set<int>> mCurvePointGraph; // neighbor point ID set : # of curve points

		vector<vector<vec3>> mChainedCurves;
		vector<vector<int>> mChainedCurvePoints;
	};
}