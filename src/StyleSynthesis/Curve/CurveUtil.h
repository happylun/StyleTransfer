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

#include "Data/StyleSynthesisTypes.h"

namespace StyleSynthesis {

	class CurveUtil {

	private:

		// make it non-instantiable
		CurveUtil() {}
		~CurveUtil() {}

	public:

		static bool makeTube(vector<vec3> &curve, TTriangleMesh &tube, float radius);
		static bool cleanCurve(vector<vec3> &inCurve, vector<vec3> &outCurve);
		static bool buildMatrix(vector<vec3> &inCurve, Eigen::Matrix3Xd &outMatrix);
		static bool computeCurveDirections(vector<vec3> &inCurve, vector<vec3> &outDirections);
		static bool computeCurveLength(vector<vec3> &inCurve, double &outLength);
		static bool checkStraightCurve(vector<vec3> &inCurve, bool &outFlag);

		static bool chainLines(vector<vec2i> &edges, vector<vec3> &vertices, vector<vector<int>> &indices, double maxAngle = 30.0);
		static bool chainCurves(vector<vec3> &vertices, vector<vector<int>> &indices, double maxAngle = 30.0);
		static bool extractLines(vector<vector<int>> &indices, vector<vec3> &vertices, vector<vector<vec3>> &lines);
		static bool removeDuplicatePoints(vector<vec2i> &edges, vector<vec3> &vertices);
		static bool removeDuplicateLines(vector<vector<vec3>> &curves);
		static bool sampleLine(vector<vec3> &curve, vector<vec3> &sample, float radius);
		static bool sampleLines(vector<vector<vec3>> &curves, vector<vector<vec3>> &samples, float radius);
		static bool filterShortLines(vector<vector<vec3>> &curves, float length);

		static bool saveCurves(string fileName, vector<vector<vec3>> &curves);
		static bool saveCurves(ostream &fileStream, vector<vector<vec3>> &curves);
		static bool loadCurves(string fileName, vector<vector<vec3>> &curves);
		static bool loadCurves(istream &fileStream, vector<vector<vec3>> &curves);
		static bool saveCurveGroups(string fileName, vector<vector<vector<vec3>>> &curves);
		static bool loadCurveGroups(string fileName, vector<vector<vector<vec3>>> &curves);
		static bool visualizeCurves(string fileName, vector<vector<vec3>> &curves, int style = 0);

		static bool generateViewPoints(vector<vec3> &viewPoints);
		static bool saveViewPoints(string fileName, vector<vec3> &viewPoints);
		static bool loadViewPoints(string fileName, vector<vec3> &viewPoints);
		
	};
}