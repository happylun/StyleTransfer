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

#include <Library/CMLHelper.h>

#include "Data/StyleSynthesisTypes.h"

using namespace std;

namespace StyleSynthesis {

	class SegmentUtil {

	private:

		// make it non-instantiable
		SegmentUtil() {}
		~SegmentUtil() {}

	public:

		static bool saveSegmentationData(string fileName, vector<vector<int>> &segmentData);
		static bool loadSegmentationData(string fileName, vector<vector<int>> &segmentData);

		static bool savePatchData(string fileName, vector<vector<int>> &patchData, vector<vector<int>> &patchGraph);
		static bool loadPatchData(string fileName, vector<vector<int>> &patchData, vector<vector<int>> &patchGraph);

		static bool saveLabelData(string fileName, vector<int> &labelData);
		static bool loadLabelData(string fileName, vector<int> &labelData);

		static bool visualizeSegmentHierarchy(string fileName, TSampleSet &sample, vector<vector<int>> &patch, vector<vector<int>> &segment);
		static bool visualizeSegmentOnMesh(string fileName, TTriangleMesh &mesh, vector<vector<int>> &segment);

		static bool buildKNNGraph(
			TSampleSet &inSamples,
			vector<vector<int>> &outGraph, // sample ID : # of neighbors : # of samples
			vector<bool> &outFlag,         // sample is inliner : # of samples
			double optRadius = 1.5,        // optional neighbor radius
			bool optNormalCheck = true     // optional normal checking
			);

		static vec3i colorMapping(int index);
	};
}