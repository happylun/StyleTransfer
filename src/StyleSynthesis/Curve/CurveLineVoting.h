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

#include "Data/StyleSynthesisTypes.h"

namespace StyleSynthesis {

	class CurveLineVoting {

		// align curves to lines

	public:

		CurveLineVoting(vector<vector<vec3>> &curves, vector<vector<vec3>> &lines);
		~CurveLineVoting();

	public:

		bool process();
		bool output(vector<vector<vec3>> &alignedLines);
		bool visualize(string fileName);

	private:

		vector<vector<vec3>> *mpCurves;
		vector<vector<vec3>> *mpLines;

		int mNumCurves;
		int mNumLines;

		vector<vector<vec3>> mAlignedLines;
	};
}