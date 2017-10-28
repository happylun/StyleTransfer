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

#include "Library/TheaKDTreeHelper.h"

#include "Data/StyleSynthesisTypes.h"

using namespace std;

namespace StyleSynthesis {

	class DeformVolumetricGraph {

	public:

		DeformVolumetricGraph();
		~DeformVolumetricGraph();

	public:

		bool loadMesh(TTriangleMesh &mesh);
		bool loadNames(string visualizeFolder);
		bool process();
		bool output(TTetrahedralMesh &mesh);

	private:

		bool latticeSampling();
		bool generateShell();
		bool tetrahedralization();
		bool extractAlphaShape();
		bool simplifyGraph();
		bool smoothGraph();
		bool washoutInternal();

		bool checkInside(vec3 &point);
		int getVertexMap(vector<int> &vmap, int id);
		inline vec2i makeKey(int a, int b) { return a < b ? vec2i(a, b) : vec2i(b, a); }
		inline vec2i makeKey(vec2i i) { return i[0] < i[1] ? i : vec2i(i[1], i[0]); }

		inline bool error(string s) { cout << "Error: " << s << endl; return false; }

	private:

		TTriangleMesh *mpMesh;

		string mVisualizeFolder;

		vec3 mBBMin, mBBMax; // bounding box size
		double mAverageEdgeLength;
		TKDTree mMeshTree;
		TKDTreeData mMeshTreeData;

		vector<vec3> mPointCloud;
		TTetrahedralMesh mTetraMesh;
	};
}