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

#include "FeatureShapeDistributions.h"

#include <fstream>
#include <iostream>
#include <omp.h>

#include "Utility/PlyExporter.h"

#include "Feature/FeatureUtil.h"

using namespace StyleSynthesis;

const int FeatureShapeDistributions::NUM_SAMPLES = 1024;
const int FeatureShapeDistributions::HIST_BINS[] = {16, 32, 64, 128};

FeatureShapeDistributions::FeatureShapeDistributions(TPointSet *points, vector<double> *features) {

	mpPoints = points;
	mpFeatures = features;
}

FeatureShapeDistributions::~FeatureShapeDistributions() {
}

bool FeatureShapeDistributions::calculate() {

	// sub-sample points
	
	int n = mpPoints->amount;
	vector<int> indices(n);
	for (int k = 0; k < n; k++) indices[k] = k;
	if (NUM_SAMPLES) {
		random_shuffle(indices.begin(), indices.end());
		n = min(n, NUM_SAMPLES);
	}

	// init for omp

#ifdef _OPENMP
	int numThreads = omp_get_max_threads();
#else
	int numThreads = 1;
#endif

	vector<vector<double>> ompDistributions(numThreads);
	for (auto &it : ompDistributions) {
		it.clear();
		it.reserve(n / numThreads * n * 2);
	}
	vector<double> ompMaxD2(numThreads, 0);

	// get all pairs

	vector<vec2i> pairs;
	pairs.reserve(n*n/2);
	for (int i1 = 0; i1 < n-1; i1++) {
		for (int i2 = i1+1; i2 < n; i2++) {
			pairs.push_back(vec2i(i1, i2));
		}
	}
	int totalCount = (int)pairs.size();

	// get shape distributions

#pragma omp parallel for num_threads(numThreads)
	for (int pairID = 0; pairID < totalCount; pairID++) {		

#ifdef _OPENMP
		int threadID = omp_get_thread_num();
#else
		int threadID = 0;
#endif

		vec3 p1 = mpPoints->positions[indices[pairs[pairID][0]]];
		vec3 p2 = mpPoints->positions[indices[pairs[pairID][1]]];
		double d = (p2 - p1).length();
		ompDistributions[threadID].push_back(d);
		ompMaxD2[threadID] = max(ompMaxD2[threadID], d);
	}

	double maxD2 = 0;
	for (int k = 0; k < numThreads; k++) {
		maxD2 = max(maxD2, ompMaxD2[k]);
	}

	// export features (concat histograms)

	vector<double> distributions;
	for (int k = 0; k < numThreads; k++) {
		distributions.insert(distributions.end(), ompDistributions[k].begin(), ompDistributions[k].end());
	}
	int numHists = sizeof(HIST_BINS) / sizeof(HIST_BINS[0]);
	mpFeatures->clear();
	for (int histID = 0; histID < numHists; histID++) {
		vector<double> histogram;
		if (!FeatureUtil::computeHistogram(distributions, histogram, HIST_BINS[histID], 0.0, maxD2)) return false;
		mpFeatures->insert(mpFeatures->end(), histogram.begin(), histogram.end());
	}

	return true;
}

bool FeatureShapeDistributions::visualize(string fileName) {

	PlyExporter pe;

	int numHists = sizeof(HIST_BINS) / sizeof(HIST_BINS[0]);
	int startPos = 0;
	for (int histID = 0; histID < numHists; histID++) {

		int numBins = HIST_BINS[histID];
		vector<double> histogram(mpFeatures->begin() + startPos, mpFeatures->begin() + startPos + numBins);

		float w = 1.0f / numBins;
		float h = numBins / 10.0f;

		vector<vec3> zeroLine;
		vector<vec3> histLine;
		vec3 currentZ(0, 0, -0.5f*histID);
		vec3 currentH(0, 0, -0.5f*histID);
		for (int j = 0; j < (int)histogram.size(); j++) {
			vec3 nextZ = currentZ + vec3(w, 0.0f, 0.0f);
			vec3 nextH = nextZ + vec3(0.0f, (float)histogram[j] * h, 0.0f);
			zeroLine.push_back(currentZ);
			zeroLine.push_back(nextZ);
			histLine.push_back(currentH);
			histLine.push_back(nextH);
			currentZ = nextZ;
			currentH = nextH;
		}

		if (!pe.addLine(&zeroLine)) return false;
		if (!pe.addLine(&histLine)) return false;

		startPos += numBins;
	}

	if (!pe.output(fileName)) return false;

	return true;
}

bool FeatureShapeDistributions::compareFeatures(vector<double> &feature1, vector<double> &feature2, vector<double> &distance) {

	distance.clear();

	int numHists = sizeof(HIST_BINS) / sizeof(HIST_BINS[0]);
	int startPos = 0;
	for (int histID = 0; histID < numHists; histID++) {

		int numBins = HIST_BINS[histID];

		vector<double> hist1(feature1.begin() + startPos, feature1.begin() + startPos + numBins);
		vector<double> hist2(feature2.begin() + startPos, feature2.begin() + startPos + numBins);
		double dist = FeatureUtil::computeEMD(hist1, hist2);
		distance.push_back(dist);
	}

	return true;
}