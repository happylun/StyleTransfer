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

	namespace AbstractionSubvolumeNamespace {
		struct TPrimitive;
	}
	typedef AbstractionSubvolumeNamespace::TPrimitive TPrim; // shorter name

	class SegmentMeshPriFit {

	public:

		SegmentMeshPriFit();
		~SegmentMeshPriFit();

	public:

		bool loadData(TTriangleMesh &mesh, vector<vector<int>> &patches);
		bool runSegmentation();
		bool exportSegmentation(vector<vector<int>> &segments);
		bool visualizeSegmentation(string fileName);
		bool visualizePrimitive(string fileName);

	private:

		bool initSegmentation();
		bool iterativeMerging();
		bool cleanupGroups();

	private:

		bool samplePrimitiveVertices(vec4i &sample, set<int> &groups);
		bool checkPrimitive(TPrim *primitive);
		bool checkFitting(
			TPrim *inPrimitive, set<int> &inGroups,
			TPointSet &outFittingPoints, double &outFittingDistance,
			bool &outValidFlag);
		bool computeFittingDistance(TPrim *primitive, TPointSet &points, double &distance);
		bool refitPrimitive(
			TPrim *&primitive, set<int> &groups,
			TPointSet &fittingPoints, double &fittingDistance);
		bool checkMerging(TPrim *primitive, set<int> &groups);
		bool mergeGroups(set<int> &groups);

		inline static bool error(string s) { cout << "Error: " << s << endl; return false; }

	private:

		TTriangleMesh mMesh;
		vector<vector<int>> mPatches; // face ID : # of faces : # of patches

		int mNumVertices;
		int mNumFaces;
		int mNumPatches;

		vector<set<int>> mVertexPatches; // patch ID : # of touching patches : # of vertices
		vector<TPointSet> mPatchPoints; // vertex set within patch : # of patches

		vector<set<int>> mPatchGroups; // patch ID : patch set : # of groups
		vector<int> mPatchMap; // group ID : # of patches

		TTriangleMesh mPrimitiveMesh;
	};
}