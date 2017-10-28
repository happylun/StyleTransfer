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

	class CurveDeformVertex {

	public:

		CurveDeformVertex();
		~CurveDeformVertex();

	public:

		bool loadMesh(TTriangleMesh &mesh);
		bool loadCurves(
			vector<vector<vector<vec3>>> &exemplarCurvesGroup,
			vector<vector<vector<vec3>>> &candidateCurvesGroup);
		bool loadNames(string visualizeFolder);

		bool process();
		bool visualize();
		bool output(TTriangleMesh &mesh);

	private:

		bool alignCurves();
		bool extractHandles();
		bool deformMesh();

		bool alignGroup(vector<vector<vec3>> &group);

		bool alignCurve(
			vector<vec3> &exemCurve,
			vector<vec3> &candCurve,
			vector<vec3d> &alignedExemPoints,
			vector<vec3d> &alignedCandPoints);

		bool interpolateGroup(
			vector<vector<vec3d>> &exemGroup,
			vector<vector<vec3d>> &candGroup);

		bool interpolateCurve(
			double alpha,
			vector<vec3d> &sourceCurve,
			vector<vec3d> &targetCurve,
			vector<vec3d> &interpolatedCurve);

		bool smoothMesh(TTriangleMesh &mesh, vector<bool> &constraintFlags);

		static inline bool error(string info) { cout << "\nError: " << info << endl; return false; }

	private:

		TTriangleMesh *mpMesh;
		vector<vector<vector<vec3>>> mExemplarCurvesGroup;
		vector<vector<vector<vec3>>> mCandidateCurvesGroup;

		string mVisualizeFolder;

		bool mNoDeformation;

		vector<vector<vec3d>> mExemplarCurvePoints;
		vector<vector<vec3d>> mCandidateCurvePoints;
		vector<bool> mInterpolatedCurveFlags; // is interpolated curve : # of all curves

		vector<bool> mDeformCurveHandle; // is handle : # of vertices
		vector<bool> mDeformHardHandle; // is "hard" handle : # of vertices
		vector<bool> mDeformFixedHandle; // is fixed : # of vertices
		vector<vec3d> mDeformDisplacement; // displacement : # of vertices

		TTriangleMesh mDeformedMesh;
	};
}