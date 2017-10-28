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
	class ContextPartGraphNodeGenerator;
	class ContextPartGraphAssembleResult;

	class ContextPartGraphTabuSearchPart {

	public:

		ContextPartGraphTabuSearchPart();
		~ContextPartGraphTabuSearchPart();

		typedef ContextPartGraph TGraph;
		typedef ContextPartGraphNode TNode;
		typedef ContextPartGraphNodeGenerator TNodeGen;
		typedef ContextPartGraphAssembleResult TSolution;

	public:

		bool loadGraph(TNodeGen *generator, TGraph *source, TGraph *target);
		bool loadMatchings(vector<vector<vec2i>> &matchings);
		bool loadContributions(
			vector<double> &sourceContribution,
			vector<double> &targetContribution,
			vector<double> &sourceSaliency,
			vector<double> &targetSaliency);
		bool loadKeyGroups(vector<vector<int>> &sourceGroups, vector<vector<int>> &targetGroups);
		bool loadNames(string weightsFolder, string resultFolder);
		bool process();

		int getNumSolutions() { return (int)mSolutionPool.size(); }

	private:

		bool initialize();
		bool runTabuSearch();
		bool analyzeSolutions();
		bool cleanUp();

		bool finalizeSolution(int solutionID);
		bool runAddingRemoving(
			vector<int> &addingList, vector<int> &removingList,
			TSolution* currentSolution, TSolution* &newSolution);
		bool runPostProcessing(TSolution* currentSolution, TSolution* &newSolution);

		bool computeFunctionalSimilarity(TSolution *solution, double &similarity);
		bool computeStyleScore(vector<int> &sourceContributors, vector<int> &targetContributors, double &score);

		inline static bool error(string s) { cout << "Error: " << s << endl; return false; }

	protected:

		TNodeGen *mpNodeGenerator;
		TGraph *mpSourceGraph;
		TGraph *mpTargetGraph;
		vector<vector<vec2i>> mMatchings;
		vector<double> mSourceContribution;
		vector<double> mTargetContribution;
		vector<double> mSourceSaliency;
		vector<double> mTargetSaliency;
		vector<vector<int>> mSourceKeyGroups;
		vector<vector<int>> mTargetKeyGroups;

		string mWeightsFolder;
		string mResultFolder;

		vector<TSolution*> mSolutionPool;
		vector<vector<int>> mSolutionSourceContributors;
		vector<vector<int>> mSolutionTargetContributors;
		vector<double> mSolutionStyleScore;

		Eigen::MatrixXd mOriginalIdentitySimilarity;
	};
}