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

#include <string>
#include <vector>

#include <Eigen/Dense>

#include "Data/StyleSynthesisTypes.h"

using namespace std;

namespace StyleSynthesis {

	class ContextPartGraph;
	class ContextPartGraphNodeGenerator;

	class ContextPartGraphCrossValidation {

		// not only cross validation
		// basically all sorts of validations we have tried...
		// I am just too lazy to change this class name...

	public:

		ContextPartGraphCrossValidation();
		~ContextPartGraphCrossValidation();

	public:

		bool train(int numSplits);
		bool test(int numSplits);
		bool evaluate();
		bool compare();

	private:

		bool validatePair(
			string sourceFolder,
			string targetFolder,
			string weightsFolder,
			vector<vec2i> &nodePairs,
			int &numCorrectPairs);

		bool validatePair(
			ContextPartGraph *sourceGraph,
			ContextPartGraph *targetGraph,
			string weightsFolder,
			vector<vec2i> &nodePairs,
			int &numCorrectPairs);

		bool computeSimilarity(
			ContextPartGraph *sourceGraph,
			ContextPartGraph *targetGraph,
			string weightsFolder,
			Eigen::MatrixXd &similarityMatrix);

		bool getDifferentMatching(
			ContextPartGraph *sourceGraph,
			ContextPartGraph *targetGraph,
			Eigen::MatrixXd &simMatOurs,
			Eigen::MatrixXd &simMatLaga,
			vector<int> &validSourceNodes,
			vector<vec3i> &validTriplets,
			vector<vec3i> &diffTriplets,
			vec3i &stats);

		bool loadMatchingData(
			string sourceName, string targetName,
			ContextPartGraphNodeGenerator &nodeGen,
			TTriangleMesh &sourceMesh,
			TTriangleMesh &targetMesh,
			vector<vector<vector<int>>> &sourceSegment,
			vector<vector<vector<int>>> &targetSegment,
			ContextPartGraph &sourceGraph,
			ContextPartGraph &targetGraph);

		bool loadMatchingPair(
			string fileName,
			vector<string> &sourceList,
			vector<string> &targetList,
			vector<vector<vec2i>> &pairList);
	};
}