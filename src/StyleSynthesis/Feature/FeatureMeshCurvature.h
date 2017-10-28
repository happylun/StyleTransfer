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

	class FeatureMeshCurvature {

		// Curvature

	public:

		struct TCurvature {
			float k1; // max curvature
			float k2; // min curvature
			vec3 d1; // max curvature direction
			vec3 d2; // min curvature direction
			// {n, d1, d2} span an orthonormal right-hand-side coordinate system
		};

	private:

		FeatureMeshCurvature() {}
		~FeatureMeshCurvature() {}

	public:

		static bool computeCurvature(TTriangleMesh &mesh, vector<TCurvature> &curvature);
		static bool computeCurvatureDerivative(TTriangleMesh &mesh, vector<TCurvature> &curvature, vector<vec4> &derivative);

		static bool smoothNormal(TTriangleMesh &mesh, vector<TCurvature> &curvature);
		static bool smoothCurvature(TTriangleMesh &mesh, vector<TCurvature> &curvature);
		static bool smoothCurvatureDerivative(TTriangleMesh &mesh, vector<TCurvature> &curvature, vector<vec4> &derivative);

	public:

		// utility methods

		static bool computeNormalCurvature(TCurvature &tensor, vec3 &direction, float &curvature);

	private:

		static bool computeFeatureSize(TTriangleMesh &mesh, vector<TCurvature> &curvature, float &featureSize);
		static bool computePointArea(TTriangleMesh &mesh, vector<float> &area);

		static bool projectCurvatureTensor(const vec3d(&oldCS)[3], vec3d oldTensor, const vec3d(&newCS)[3], vec3d &newTensor);
		static bool projectCurvatureDerivative(const vec3d(&oldCS)[3], vec4d oldDerv, const vec3d(&newCS)[3], vec4d &newDerv);
		static bool computeCornerArea(const vec3d(&faceEdge)[3], vec3d &cornerArea);
	};
}