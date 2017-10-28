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

#include "CircularityDescriptor.h"

#include <iostream>

#include "constants.h"

using namespace LFD;

bool CircularityDescriptor::calculate(
	vector<unsigned char> &inImage,
	double &outDescriptor)
{
	double a = calculateArea(inImage);
	double p = calculatePerimeter(inImage);

	if (p > 0) {
		outDescriptor = 4 * MY_PI * a / (p * p);
		if (outDescriptor > 1) outDescriptor = 1;
	} else {
		outDescriptor = 0;
	}

	return true;
}

double CircularityDescriptor::calculateArea(vector<unsigned char> &inImage) {

	int count = 0;
	for (int y = 0; y < WINDOW_HEIGHT; y++) {
		for (int x = 0; x < WINDOW_WIDTH; x++) {
			if (inImage[y*WINDOW_WIDTH + x] < 255) {
				count++;
			}
		}
	}

	return (double)count;
}

double CircularityDescriptor::calculatePerimeter(vector<unsigned char> &inImage) {

	int offset[8] = { -1, 1,
		-WINDOW_WIDTH, WINDOW_WIDTH,
		-WINDOW_WIDTH - 1, -WINDOW_WIDTH + 1,
		WINDOW_WIDTH - 1, WINDOW_WIDTH + 1 };

	int count = 0;
	for (int y = 1; y < WINDOW_HEIGHT-1; y++) {
		for (int x = 1; x < WINDOW_WIDTH-1; x++) {
			int pos = y*WINDOW_WIDTH + x;
			if (inImage[pos] < 255) {
				for (int k = 0; k < 8; k++) {
					if (inImage[pos + offset[k]] == 255) {
						count++;
						break;
					}
				}
			}
		}
	}

	return (double)count;
}