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

#include <Eigen/Sparse>

#include "Data/StyleSynthesisTypes.h"

using namespace std;

namespace StyleSynthesis {

	class SegmentMeshApxCvx {

	public:

		SegmentMeshApxCvx();
		~SegmentMeshApxCvx();

	public:

		bool initSegmentation(TTriangleMesh &mesh);
		bool runSegmentation(double visibility);
		bool exportMesh(TTriangleMesh &outMesh);
		bool exportSegmentation(vector<vector<int>> &outSegments);
		bool visualizeSegmentation(string fileName);

	private:

		bool extractPatch();
		bool initPatchGraph();
		bool computePatchVisibility();
		bool processConvexNeighbors();

		bool splitPlanarPatch(vector<int> &inPatch, vector<vector<int>> &outPatches);
		bool computeSegmentVisibility(vector<int> &seg1, vector<int> &seg2, double &visibility);

		inline static vec2i makeKey(int a, int b) { return a < b ? vec2i(a, b) : vec2i(b, a); }
		inline static bool error(string s) { cout << "Error: " << s << endl; return false; }

	private:

		TTriangleMesh mMesh;

		int mNumPatches;
		int mNumSegments;

		vector<vector<int>> mFaceGraph; // neighbor face ID : # of neighbors : # of faces
		vector<vector<int>> mPatches; // face ID : # of faces : # of patches
		vector<set<int>> mPatchGraph; // neighboring patch ID set : # of patches
		vector<double> mPatchSize; // total face area : # of patches
		Eigen::MatrixXd mPatchVisibility; // visibility in [0,1] : (# of patches) X (# of patches)
		Eigen::MatrixXi mPatchAdjacency; // (0: not adjacent; 1: adjacent) : (# of patches) X (# of patches)
		Eigen::MatrixXi mPatchConvexAdjacency; // 0/1 : (# of patches) X (# of patches)

		vector<vector<int>> mSegments; // patch ID : # of patches : # of segments
	};
}