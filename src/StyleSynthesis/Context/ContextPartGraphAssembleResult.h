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

#include "ContextPartGraph.h"

#include "Data/StyleSynthesisTypes.h"

using namespace std;

namespace StyleSynthesis {

	class ContextPartGraphAssemble;

	class ContextPartGraphAssembleResult {

	public:

		ContextPartGraphAssembleResult();
		~ContextPartGraphAssembleResult();

		typedef ContextPartGraph TGraph;
		typedef ContextPartGraphNode TNode;
		typedef ContextPartGraphNodeGenerator TNodeGen;
		typedef ContextPartGraphAssemble TAssemble;

	public:

		bool identity(TGraph *graph);
		bool generate(TAssemble *assemble);
		bool chain(ContextPartGraphAssembleResult *previousResult);
		bool visualize(string fileName);
		bool visualize(string fileName, ContextPartGraphAssembleResult *prevSolution);
		bool saveData(string folderName);
		bool loadData(string folderName, TNodeGen &nodeGen);

	private:

		inline static bool error(string s) { cout << "Error: " << s << endl; return false; }

	public:

		TGraph mGraph;
		TTriangleMesh mMesh;
		vector<vector<vector<int>>> mSegment;

		vector<int> mNodeMapping; // node ID in new graph (-1 if replaced or removed or modified) : # of nodes in target graph
		vector<Eigen::Affine3d> mNodeTransformation; // transformation applied on node (identity if replaced or removed or modified) : # of nodes in target graph

		vector<vec2i> mReplaceMapping; // (source node ID, target node ID) (both ID are -1 if not "replaced"; target ID is -1 if "added") : # of nodes in new graph
		vector<Eigen::Affine3d> mReplaceTransformation; // transformation applied on source node to generate this new node (identity if not replaced) : # of nodes in new graph
	};
}