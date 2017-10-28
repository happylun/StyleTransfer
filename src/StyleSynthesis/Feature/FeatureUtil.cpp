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

#include "FeatureUtil.h"

#include <fstream>

#include "Library/CMLHelper.h"

using namespace StyleSynthesis;

bool FeatureUtil::computeHistogram(
	vector<double> &distribution, vector<double> &histogram,
	int numBins, double minValue, double maxValue, bool smooth)
{
	if (minValue > maxValue) {
		cout << "Error: incorrect min/max" << endl;
		return false;
	}
	if (minValue == maxValue) {
		double weight = 1.0 / numBins;
		histogram.assign(numBins, weight);
		return true;
	}

	const double SIGMA = 2.0; // Gaussian deviation (number of bins)
	const double DENOM = 2*SIGMA*SIGMA;

	histogram.clear();
	histogram.resize(numBins, 0);
	double totalWeight = 0;

	// compute smooth histogram
	for (double value : distribution) {
		double pos = (value - minValue) / (maxValue - minValue) * numBins;
		if (smooth) {			
			int lbin = max(0, (int)(pos - SIGMA*3)); // only count contribution within 3 standard deviation
			int rbin = min(numBins-1, (int)(pos + SIGMA*3));
			for (int bin = lbin; bin <= rbin; bin++) {
				double binCenter = bin + 0.5;
				double weight = exp(-cml::sqr(pos - binCenter) / DENOM);
				histogram[bin] += weight;
				totalWeight += weight;
			}
		} else {
			int bin = (int)pos;
			histogram[bin] += 1.0;
		}
	}
	if (!smooth) totalWeight = (int)distribution.size();

	// normalize histogram
	if (totalWeight) {
		for (double &weight : histogram) {
			weight /= totalWeight;
		}
	} else {
		double value = 1.0 / numBins;
		for (double &weight : histogram) {
			weight  = value;
		}
	}

	return true;
}

double FeatureUtil::computeEMD(vector<double> &histogram1, vector<double> &histogram2) {

	// Earth Mover's Distance

	if (histogram1.size() != histogram2.size()) {
		cout << "Error: unmatched sizes of histogram" << endl;
		return false;
	}

	int n = (int)histogram1.size();

	// compute CDF
	vector<double> cdf1(n, 0);
	vector<double> cdf2(n, 0);
	cdf1[0] = histogram1[0];
	cdf2[0] = histogram2[0];
	for (int i = 1; i < n; i++) {
		cdf1[i] = cdf1[i - 1] + histogram1[i];
		cdf2[i] = cdf2[i - 1] + histogram2[i];
	}

	// compute L1 distance of CDFs
	return computeL1D(cdf1, cdf2);
}

double FeatureUtil::computeL1D(vector<double> &histogram1, vector<double> &histogram2) {

	// L1 distance

	if (histogram1.size() != histogram2.size()) {
		cout << "Error: unmatched sizes of histogram" << endl;
		return false;
	}

	int n = (int)histogram1.size();

	double dist = 0;
	for (int i = 0; i < n; i++) {
		dist += fabs(histogram1[i] - histogram2[i]);
	}

	return dist;
}

vec3i FeatureUtil::colorMapping(double v) {

	// ref: http://stackoverflow.com/questions/7706339/grayscale-to-red-green-blue-matlab-jet-color-scale

	// 0 => 1
	// blue => cyan => yellow => red

	int b = v < 0.125 ? int((v + 0.125) * 4 * 255) : v < 0.375 ? 255 : v<0.625 ? int((0.625 - v) * 4 * 255) : 0;
	int g = v < 0.125 ? 0 : v < 0.375 ? int((v - 0.125) * 4 * 255) : v < 0.625 ? 255 : v<0.875 ? int((0.875 - v) * 4 * 255) : 0;
	int r = v < 0.375 ? 0 : v < 0.625 ? int((v - 0.375) * 4 * 255) : v < 0.875 ? 255 : int((1.125 - v) * 4 * 255);

	return vec3i(r, g, b);
}