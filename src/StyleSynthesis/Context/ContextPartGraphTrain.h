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

	class ContextPartGraphTrain {

	public:

		ContextPartGraphTrain();
		~ContextPartGraphTrain();

		typedef ContextPartGraphNode TNode;
		typedef ContextPartGraph TGraph;

	public:

		bool loadSceneMesh(vector<TTriangleMesh> &sceneMeshes);
		bool loadSceneSegments(vector<vector<vector<vector<int>>>> &sceneSegments);
		bool loadSceneGraph(vector<TGraph> &sceneGraphs);
		bool process();
		bool visualizeElements(string fileName);
		bool visualizeSymmetries(string fileName);
		bool visualizePartPairs(string fileName);
		bool outputElements(string fileName);
		bool outputPartPairs(string fileName);
		bool outputSigma(string nodeSigmaName, string edgeSigmaName);

	private:

		bool generateElementGroups();
		bool generateSymmetryGroups();
		bool extractTrainPartPairs();

		bool checkNodes(TNode *sourceNode, TNode *targetNode, bool &isElement);

		inline static bool error(string s) { cout << "Error: " << s << endl; return false; }

	private:

		vector<TTriangleMesh> *mpSceneMeshes;
		vector<vector<vector<vector<int>>>> *mpSceneSegments;
		vector<TGraph> *mpSceneGraphs;

		vector<vector<vec2i>> mElementGroup; // (mesh ID, node ID) : # of individual elements : # of element groups
		vector<vector<vec2i>> mSymmetryGroup; // (mesh ID, node ID) : # of individual elements : # of element groups

		vector<vector<vec2i>> mTrainPartPairs; // (node ID, node ID) : # of part pairs : # of mesh pairs
		vector<vec2i> mTrainMeshPairs; // (mesh ID, mesh ID) : # of mesh pairs

		vector<double> mSceneNodeSigma;
		vector<double> mSceneEdgeSigma;
	};
}