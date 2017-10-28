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

#include "Library/CMLHelper.h"
#include "Library/TheaKDTreeHelper.h"

#include "Data/StyleSynthesisTypes.h"

using namespace std;

namespace StyleSynthesis {

	class FeatureShapeDiameterFunctions {

		// Shape Diameter Function

	public:

		FeatureShapeDiameterFunctions(TSampleSet *samples, TTriangleMesh *mesh, vector<double> *features);
		FeatureShapeDiameterFunctions(TSampleSet *samples, TKDTree *meshTree, vector<double> *features);
		~FeatureShapeDiameterFunctions();
	
	public:

		bool calculate();
		bool visualize(string fileName);

		// 4 hists
		static bool compareFeatures(vector<double> &feature1, vector<double> &feature2, vector<double> &distance);

	private:

		bool calculateSDF();
		bool smoothSDF();
		bool normalizeSDF();

	private:

		TSampleSet *mpSamples;
		TKDTree *mpMeshTree;
		vector<double> *mpFeatures;

		TKDTreeData mTreeData;
	};
}