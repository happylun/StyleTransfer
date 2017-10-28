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

#include "Data/StyleSynthesisTypes.h"

using namespace std;

namespace StyleSynthesis {

	class FeatureShapeDistributions {

		// Shape Distributions

	private:

		const static int NUM_SAMPLES;
		const static int HIST_BINS[];

	public:

		FeatureShapeDistributions(TPointSet *points, vector<double> *features);
		~FeatureShapeDistributions();
	
	public:

		bool calculate();
		bool visualize(string fileName);

		// 4 hists
		static bool compareFeatures(vector<double> &feature1, vector<double> &feature2, vector<double> &distance);

	private:

		TPointSet *mpPoints;
		vector<double> *mpFeatures;
	};
}