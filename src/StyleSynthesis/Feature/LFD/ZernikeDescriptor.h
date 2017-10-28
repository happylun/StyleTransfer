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

using namespace std;

namespace LFD {
	class ZernikeDescriptor {

	private:

		static const int BASIS_RADIUS = 50;
		static const int BASIS_LUT_SIZE = BASIS_RADIUS * 2 + 1;
		static const int BASIS_ANGULAR = 12;
		static const int BASIS_RADIAL = 3;

	public:

		static const int DESCRIPTOR_SIZE = BASIS_ANGULAR*BASIS_RADIAL - 1;

	private:

		ZernikeDescriptor() {}
		~ZernikeDescriptor() {}

	public:

		static bool init();

		static bool radius(
			vector<unsigned char> &inImage,
			pair<double, double> &inCenter,
			double &outRadius);

		static bool calculate(
			vector<unsigned char> &inImage,
			pair<double,double> &inCenter,
			double &inRadius,
			vector<double> &outDescriptor);

	private:

		static double mBasisR[BASIS_ANGULAR][BASIS_RADIAL][BASIS_LUT_SIZE][BASIS_LUT_SIZE];
		static double mBasisI[BASIS_ANGULAR][BASIS_RADIAL][BASIS_LUT_SIZE][BASIS_LUT_SIZE];
	};
}