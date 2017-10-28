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
#include <algorithm>

using namespace std;

namespace LFD {
	class FourierDescriptor {

	public:

		static const int DESCRIPTOR_SIZE = 10;

	private:

		FourierDescriptor() {}
		~FourierDescriptor() {}

	public:

		static bool calculate(
			vector<unsigned char> &inImage,
			pair<double, double> &inCenter,
			vector<double> &outDescriptor);

	private:

		static bool extractContour(
			vector<unsigned char> &inImage,
			vector<int> &outContour);

		static bool getDescriptor(
			vector<int> &inContour,
			pair<double, double> &inCenter,
			vector<double> &outDescriptor);

		// morphological processing
		static bool isConnected(vector<unsigned char> &image);
		static bool connect(vector<unsigned char> &image);
		static bool erosion(vector<unsigned char> &image, vector<int> &mask);
		static bool dilation(vector<unsigned char> &image, vector<int> &mask);
	};
}