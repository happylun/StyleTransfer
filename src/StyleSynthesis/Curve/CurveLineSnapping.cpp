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

#include "CurveLineSnapping.h"

#include <fstream>
#include <set>

#include "Mesh/MeshUtil.h"
#include "Sample/SampleUtil.h"
#include "Segment/SegmentUtil.h"
#include "Curve/CurveUtil.h"

#include "Utility/PlyExporter.h"

#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

CurveLineSnapping::CurveLineSnapping(vector<vector<vec3>> &curves, vector<vector<vec3>> &lines) {

	mpCurves = &curves;
	mpLines = &lines;

	mNumCurves = (int)curves.size();
	mNumLines = (int)lines.size();
}

CurveLineSnapping::~CurveLineSnapping() {
}

bool CurveLineSnapping::process() {

	// compute thresholds

	vec3 bbMin(FLT_MAX, FLT_MAX, FLT_MAX);
	vec3 bbMax(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	for (auto &curve : *mpCurves) {
		for (vec3 &v : curve) {
			bbMin.minimize(v);
			bbMax.maximize(v);
		}
	}
	for (auto &line : *mpLines) {
		for (vec3 &v : line) {
			bbMin.minimize(v);
			bbMax.maximize(v);
		}
	}
	float bbLen = (bbMax - bbMin).length();

	double sampleRadius = bbLen * 0.001; // UNDONE: param line sampling radius
	double snapMaxDistance = bbLen * 0.01; // UNDONE: param max average sample to curve distance for snapping
	double searchMaxDistance = bbLen * 0.1;

	// build Kd tree for curve points

	vector<vector<vec3>> curvePoints;
	if (!CurveUtil::sampleLines(*mpCurves, curvePoints, (float)sampleRadius)) return false;
	vector<vec3> allCurvePoints;
	for (auto &curve : curvePoints) {
		for (vec3 &v : curve) {
			allCurvePoints.push_back(v);
		}
	}
	SKDTree curveTree;
	SKDTreeData curveTreeData;
	if (!SampleUtil::buildKdTree(allCurvePoints, curveTree, curveTreeData)) return false;

	// test for each line

	vector<bool> lineFlags(mNumLines, false);
	vector<vector<vec3>> linePoints;
	if (!CurveUtil::sampleLines(*mpLines, linePoints, (float)sampleRadius)) return false;

	for (int lineID = 0; lineID < mNumLines; lineID++) {

		vector<vec3> &lineSamples = linePoints[lineID];
		int numLineSamples = (int)lineSamples.size();

		// compute sample score

		vector<double> lineSampleDist(numLineSamples);
#pragma omp parallel for
		for (int lineSampleID = 0; lineSampleID < numLineSamples; lineSampleID++) {

			double dist = searchMaxDistance;

			// get nearest sample distance
			vec3 point = lineSamples[lineSampleID];
			SKDT::NamedPoint queryPoint(point[0], point[1], point[2]);
			Thea::BoundedSortedArray<SKDTree::Neighbor> queryResult(1);
			curveTree.kClosestElements<Thea::MetricL2>(queryPoint, queryResult, searchMaxDistance);
			if (!queryResult.isEmpty()) {
				int curvePointID = (int)curveTree.getElements()[queryResult[0].getIndex()].id;
				dist = (double)(allCurvePoints[curvePointID] - lineSamples[lineSampleID]).length();
			}

			lineSampleDist[lineSampleID] = dist;
		}

		// compute line score

		double lineDist = 0;
		for (int k = 0; k < numLineSamples; k++) {
			lineDist += lineSampleDist[k];
		}
		lineDist /= numLineSamples;

		if (lineDist < snapMaxDistance) {
			lineFlags[lineID] = true;
		}
	}

	
	// extract snapped lines

	mSnapLines.clear();
	for (int lineID = 0; lineID < mNumLines; lineID++) {
		if (lineFlags[lineID]) {
			mSnapLines.push_back((*mpLines)[lineID]);
		}
	}

	return true;
}

bool CurveLineSnapping::output(vector<vector<vec3>> &snapLines) {

	snapLines = mSnapLines;

	return true;
}

bool CurveLineSnapping::visualize(string fileName) {

	if (!CurveUtil::visualizeCurves(fileName, mSnapLines, 0)) return false;

	return true;
}