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

	class PipelineSimilarityIO : public BaseIO {

	public:

		static bool process();

	private:

		static bool learningPhase();
		static bool evaluatingPhase();
		static bool rankingPhase();

		static bool runModelPairs(string sourceName, string targetName);
		static bool runModelTriplets(string meshNameA, string meshNameB, string meshNameC, int tripletID);

		static bool evaluateSingleModel(string meshName);
		static bool evaluateModelPairs(string sourceName, string targetName, double &distance);

		static bool error(string s);
	};

}