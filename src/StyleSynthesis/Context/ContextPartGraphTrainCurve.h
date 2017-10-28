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

#include <Eigen/Dense>

#include "ContextPartGraph.h"

#include "Data/StyleSynthesisTypes.h"

using namespace std;

namespace StyleSynthesis {

	class ContextPartGraphTrainCurve {

	public:

		ContextPartGraphTrainCurve();
		~ContextPartGraphTrainCurve();

	public:

		bool loadSceneCurves(vector<string> &sceneCurveFolders, vector<vec3> &viewPoints);
		bool process();
		bool visualizeCurvePairs(string fileName, vector<TTriangleMesh> &meshes);
		bool outputCurvePairs(string fileName);


	private:

		bool initCurves();
		bool extractCurvePairs();

		bool checkCurvePair(vector<vec3> &curve1, vector<vec3> &curve2);

		inline static bool error(string s) { cout << "Error: " << s << endl; return false; }

	private:

		// curve type: 0 - ridge/valley; 1 - boundary; 2+ - contours
		vector<vector<vec2i>> mSceneCurveIndices; // (start index, end index) : # of curve types : # of meshes
		vector<vector<vector<vec3>>> mSceneCurves; // curve point : # of curve points : # of curves (mixed types) : # of meshes
		vector<vector<bool>> mSceneCurveFlags; // valid flag : # of curves (mixed types) : # of meshes

		vector<vec3> mCurveViewPoints;

		vector<vector<vec2i>> mCurvePairs; // (curve ID, curve ID) : # of curve pairs : # of mesh pairs
		vector<vec2i> mMeshPairs; // (mesh ID, mesh ID) : # of mesh pairs
	};
}