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

#include "IO/BaseIO.h"

using namespace std;

namespace StyleSynthesis {

	class ContextPartGraph;

	class PipelineMatchPartIO : public BaseIO {

	public:

		static bool process();

	private:

		static bool generatePairs(string pairFileName);

		static bool runModelPairs(string sourceName, string targetName, string pairName);
		static bool organizeResults(string sourceName, string targetName, string pairName);

	public:

		static bool computeNodeContribution(
			string sourceName, string targetName,
			ContextPartGraph *sourceGraph,
			ContextPartGraph *targetGraph,
			vector<vector<vector<int>>> *sourceSegments,
			vector<vector<vector<int>>> *targetSegments,
			double &shapeStyleDistance,
			vector<double> &sourceNodeContribution,
			vector<double> &targetNodeContribution,
			vector<double> &sourceNodeSaliency,
			vector<double> &targetNodeSaliency);

		static bool getKeyGroups(
			ContextPartGraph *graph,
			vector<double> &nodeContributions,
			vector<vector<int>> &outGroups,
			vector<double> &outGroupContributions);

		static bool error(string s);
	};

}