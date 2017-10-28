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

#include "EccentricityDescriptor.h"

#include <iostream>

#include "constants.h"

using namespace LFD;

bool EccentricityDescriptor::calculate(
	vector<unsigned char> &inImage,
	pair<double, double> &inCenter,
	double &outDescriptor)
{
	double i11 = 0;
	double i02 = 0;
	double i20 = 0;

	for (int y = 0; y < WINDOW_HEIGHT; y++) {
		for (int x = 0; x < WINDOW_WIDTH; x++) {
			if (inImage[y*WINDOW_WIDTH + x] < 255) {
				double dx = x - inCenter.first;
				double dy = y - inCenter.second;
				i02 += dy*dy;
				i11 += dx*dy;
				i20 += dx*dx;
			}
		}
	}

	if (i20+i02 != 0) {
		outDescriptor = ((i20 - i02)*(i20 - i02) + 4 * i11*i11) / ((i20 + i02)*(i20 + i02));
		if (outDescriptor > 1) outDescriptor = 1;
	} else {
		outDescriptor = 0;
	}

	return true;
}