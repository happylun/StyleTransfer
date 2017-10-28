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

#include "ZernikeDescriptor.h"

#include <iostream>

#include "constants.h"

using namespace LFD;

double ZernikeDescriptor::mBasisR[BASIS_ANGULAR][BASIS_RADIAL][BASIS_LUT_SIZE][BASIS_LUT_SIZE];
double ZernikeDescriptor::mBasisI[BASIS_ANGULAR][BASIS_RADIAL][BASIS_LUT_SIZE][BASIS_LUT_SIZE];

bool ZernikeDescriptor::init() {

	int maxRadius = BASIS_RADIUS;
	for (int y = 0; y < BASIS_LUT_SIZE; y++) {
		for (int x = 0; x < BASIS_LUT_SIZE; x++) {
			double radius = sqrt((x - maxRadius)*(x - maxRadius) + (y - maxRadius)*(y - maxRadius));
			if (radius < maxRadius) {
				double angle = atan2(y - maxRadius, x - maxRadius);
				for (int p = 0; p < BASIS_ANGULAR; p++) {
					for (int r = 0; r < BASIS_RADIAL; r++) {
						double temp = cos(radius*MY_PI*r / maxRadius);
						mBasisR[p][r][x][y] = temp*cos(angle*p);
						mBasisI[p][r][x][y] = temp*sin(angle*p);
					}
				}
			} else {
				for (int p = 0; p < BASIS_ANGULAR; p++) {
					for (int r = 0; r < BASIS_RADIAL; r++) {
						mBasisR[p][r][x][y] = 0;
						mBasisI[p][r][x][y] = 0;
					}
				}
			}
		}
	}

	return true;
}

bool ZernikeDescriptor::radius(
	vector<unsigned char> &inImage,
	pair<double, double> &inCenter,
	double &outRadius)
{

	outRadius = 0;
	for (int y = 0; y < WINDOW_HEIGHT; y++) {
		for (int x = 0; x < WINDOW_WIDTH; x++) {
			if (inImage[y*WINDOW_WIDTH + x] < 255) {
				double dx = x - inCenter.first;
				double dy = y - inCenter.second;
				double radius = sqrt(dx*dx + dy*dy);
				if (radius > outRadius) outRadius = radius;
			}
		}
	}

	return true;
}

bool ZernikeDescriptor::calculate(
	vector<unsigned char> &inImage,
	pair<double, double> &inCenter,
	double &inRadius,
	vector<double> &outDescriptor)
{

	double coeffR[BASIS_ANGULAR][BASIS_RADIAL];
	double coeffI[BASIS_ANGULAR][BASIS_RADIAL];
	memset(coeffR, 0, BASIS_ANGULAR * BASIS_RADIAL * sizeof(coeffR[0][0]));
	memset(coeffI, 0, BASIS_ANGULAR * BASIS_RADIAL * sizeof(coeffI[0][0]));

	int count = 0;
	for (int y = 0; y < WINDOW_HEIGHT; y++) {
		for (int x = 0; x < WINDOW_WIDTH; x++) {
			if (inImage[y*WINDOW_WIDTH + x] < 255) {
				double dx = x - inCenter.first;
				double dy = y - inCenter.second;
				double tx = dx/inRadius * BASIS_RADIUS + BASIS_RADIUS;
				double ty = dy/inRadius * BASIS_RADIUS + BASIS_RADIUS;
				int ix = (int)tx; // integer part
				int iy = (int)ty;
				double fx = tx - ix; // fractional part
				double fy = ty - iy;

				for (int p = 0; p < BASIS_ANGULAR; p++) {
					for (int r = 0; r < BASIS_RADIAL; r++) {
						double x1R = mBasisR[p][r][ix][iy] + (mBasisR[p][r][ix + 1][iy] - mBasisR[p][r][ix][iy]) * fx;
						double x2R = mBasisR[p][r][ix][iy + 1] + (mBasisR[p][r][ix + 1][iy + 1] - mBasisR[p][r][ix][iy + 1]) * fx;
						coeffR[p][r] += (x1R + (x2R - x1R) * fy);

						double x1I = mBasisI[p][r][ix][iy] + (mBasisI[p][r][ix + 1][iy] - mBasisI[p][r][ix][iy]) * fx;
						double x2I = mBasisI[p][r][ix][iy + 1] + (mBasisI[p][r][ix + 1][iy + 1] - mBasisI[p][r][ix][iy + 1]) * fx;
						coeffI[p][r] -= (x1I + (x2I - x1I) * fy);
					}
				}

				count++;
			}
		}
	}

	if (count > 0) {
		double coeff[BASIS_ANGULAR][BASIS_RADIAL];
		for (int p = 0; p < BASIS_ANGULAR; p++) {
			for (int r = 0; r < BASIS_RADIAL; r++) {
				double vr = coeffR[p][r] / count;
				double vi = coeffI[p][r] / count;
				coeff[p][r] = sqrt(vr*vr + vi*vi);
			}
		}
		outDescriptor.clear();
		for (int p = 0; p < BASIS_ANGULAR; p++) {
			for (int r = 0; r < BASIS_RADIAL; r++) {
				if (p != 0 || r != 0) outDescriptor.push_back(coeff[p][r] / coeff[0][0]);
			}
		}
	} else {
		outDescriptor.clear();
		outDescriptor.resize(DESCRIPTOR_SIZE, 0);
	}

	return true;
}