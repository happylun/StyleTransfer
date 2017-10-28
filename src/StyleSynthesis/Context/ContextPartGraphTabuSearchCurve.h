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

#include <Eigen/Dense>

#include "Data/StyleSynthesisTypes.h"

using namespace std;

namespace StyleSynthesis {

	class ContextPartGraph;
	class ContextPartGraphNode;
	class ContextPartGraphAssembleResult;

	class ContextPartGraphTabuSearchCurve {

	public:

		ContextPartGraphTabuSearchCurve();
		~ContextPartGraphTabuSearchCurve();

		typedef ContextPartGraph TGraph;
		typedef ContextPartGraphNode TNode;
		typedef ContextPartGraphAssembleResult TSolution;

	public:

		bool loadGraph(TGraph *source, TGraph *target);
		bool loadAssemble(TSolution *solution);
		bool loadCurve(
			vector<vector<vector<vec3>>> &srcCurves,
			vector<vector<vector<vec3>>> &tgtCurves,
			vector<vector<int>> &srcNodes,
			vector<vector<int>> &tgtNodes);
		bool loadScore(
			vector<double> &sourceContrib,
			vector<double> &targetContrib,
			vector<double> &sourceSaliency,
			vector<double> &targetSaliency);
		bool loadNames(string solutionFolder);

		bool process();

		inline bool getSolutionScore(vector<double> &score) { score = mSolutionScore; return true; }

	private:

		bool initialize();
		bool runTabuSearch();

		bool subdivideMesh(TTriangleMesh &inMesh, TTriangleMesh &outMesh);
		bool computeSolutionScore(set<int> &srcNodes, set<int> &tgtNodes, double &outScore);

		inline static bool error(string s) { cout << "Error: " << s << endl; return false; }

	protected:

		TGraph *mpSourceGraph;
		TGraph *mpTargetGraph;
		
		TSolution *mpAssembleSolution;

		vector<vector<vector<vec3>>> mMatchingSourceCurves; // curve point : # points on source curve : # of curve pairs : # of matching group
		vector<vector<vector<vec3>>> mMatchingTargetCurves; // curve point : # points on target curve : # of curve pairs : # of matching group
		vector<vector<int>> mMatchingSourceNodes; // source node ID : # of nodes in this group  : # of matching group
		vector<vector<int>> mMatchingTargetNodes; // target node ID : # of nodes in this group  : # of matching group

		vector<double> mSourceNodeContribution;
		vector<double> mTargetNodeContribution;
		vector<double> mSourceNodeSaliency;
		vector<double> mTargetNodeSaliency;

		vector<double> mSolutionScore;
		string mSolutionFolder;
	};
}