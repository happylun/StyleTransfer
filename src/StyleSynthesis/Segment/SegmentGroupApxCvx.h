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

#include "Data/StyleSynthesisTypes.h"

using namespace std;

namespace StyleSynthesis {

	class SegmentGroupApxCvx {

		// segmentation wrapper to be used in pipeline

	private:

		SegmentGroupApxCvx() {}
		~SegmentGroupApxCvx() {}

	public:

		static bool runSegmentation(
			vector<TTriangleMesh> &inGroups,
			TTriangleMesh &outMesh,
			vector<vector<vector<int>>> &outSegments); // face ID : # of faces : # of segments : # of levels

		static bool extractUniqueSegments(
			TTriangleMesh &inMesh,
			vector<vector<vector<int>>> &inSegments,
			vector<vector<int>> &outSegments);

		static bool findElementMapping(
			vector<vector<vector<int>>> &inSegments,
			vector<vector<int>> &inElements,
			vector<vector<int>> &outMapping); // element ID : # of segments : # of levels

		static bool visualize(
			string fileName,
			TTriangleMesh &mesh,
			vector<vector<vector<int>>> &segments);

		static bool saveSegments(string fileName, vector<vector<vector<int>>> &segments);
		static bool loadSegments(string fileName, vector<vector<vector<int>>> &segments);
	};
}