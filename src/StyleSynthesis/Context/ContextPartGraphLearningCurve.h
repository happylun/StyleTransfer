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

#include "Data/StyleSynthesisTypes.h"

#include "Library/cppoptlibHelper.h"

using namespace std;

namespace StyleSynthesis {

	class ContextPartGraph;
	class ContextPartGraphNodeGenerator;
	class ContextPartGraphMatchCurve;

	class ContextPartGraphLearningCurve : public cppoptlib::Problem < double > {

	public:
		ContextPartGraphLearningCurve();
		~ContextPartGraphLearningCurve();

	public:
		bool loadData(
			vector<string> &trainScenes,
			vector<vector<string>> &trainMeshes,
			map<string, int> &meshLabelMap,
			vec2i labelPair,
			int &outNumMeshPairs,
			int &outNumCurvePairs);
		bool process();
		bool computeThreshold(double &threshold);
		bool outputSigmas(string fileName);
		bool outputWeights(string fileName);

	private:

		bool initData();
		bool initLearning();
		bool runSolver();

		bool getNumCurves(string inMeshName, vec3i &outNumCurves);
		bool initGraph( // there are a bit too many arguments; but I can't do it more elegantly...
			string sourceName, string targetName,
			TTriangleMesh &sourceMesh, TTriangleMesh &targetMesh,
			vector<vector<vector<int>>> &sourceSegments,
			vector<vector<vector<int>>> &targetSegments,
			ContextPartGraph *sourceGraph, ContextPartGraph *targetGraph,
			ContextPartGraphNodeGenerator *nodeGen,
			Eigen::MatrixXd &simMat);
		bool updateWeights();
		int numErrors();

	private:

		// overriding methods
		double value(const Eigen::VectorXd &x);
		void gradient(Eigen::MatrixXd &dx);

		inline bool error(string s) { cout << "Error: " << s << endl; return false; }

	protected:

		Eigen::VectorXd mSigmas;
		Eigen::VectorXd mWeights;

	private:

		vector<pair<string, string>> mMeshPairs; // (mesh name 1, mesh name 2) : # of mesh pairs
		vector<vector<vec2i>> mCurvePairs; // (curveID1, curveID2) : # of curve pairs : # of mesh pairs
		vector<ContextPartGraphMatchCurve*> mCurveMatchings; // CPG Match Curve instance : # of mesh pairs
	};
}