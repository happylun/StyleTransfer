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

#include "CurveLineVoting.h"

#include <fstream>
#include <set>

#include "Mesh/MeshUtil.h"
#include "Sample/SampleUtil.h"
#include "Segment/SegmentUtil.h"
#include "Curve/CurveUtil.h"

#include "Utility/PlyExporter.h"

#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

CurveLineVoting::CurveLineVoting(vector<vector<vec3>> &curves, vector<vector<vec3>> &lines) {

	mpCurves = &curves;
	mpLines = &lines;

	mNumCurves = (int)curves.size();
	mNumLines = (int)lines.size();
}

CurveLineVoting::~CurveLineVoting() {
}

bool CurveLineVoting::process() {

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

	double sampleRadius = bbLen * 0.001; // UNDONE: param sample radius (also max alignment radius)
	double alignThreshold = 0.2; // UNDONE: param minimum alignment percentage threshold

	// build Kd tree for line points

	vector<vector<vec3>> linePoints;
	if (!CurveUtil::sampleLines(*mpLines, linePoints, (float)sampleRadius)) return false;
	vector<vec3> allLinePoints(0);
	vector<int> allLinePointSource(0); // line ID : # of line points
	for (int lineID = 0; lineID < mNumLines; lineID++) {
		auto &line = linePoints[lineID];
		for (vec3 &v : line) {
			allLinePoints.push_back(v);
			allLinePointSource.push_back(lineID);
		}
	}
	SKDTree lineTree;
	SKDTreeData lineTreeData;
	if (!SampleUtil::buildKdTree(allLinePoints, lineTree, lineTreeData)) return false;

	// vote by curve points

	vector<bool> allLinePointFlag(allLinePoints.size(), false);

	vector<vector<vec3>> curvePoints;
	if (!CurveUtil::sampleLines(*mpCurves, curvePoints, (float)sampleRadius)) return false;

	for (auto &curve : curvePoints) {
		int numCurveSamples = (int)curve.size();
#pragma omp parallel for
		for (int sampleID = 0; sampleID < numCurveSamples; sampleID++) {
			vec3 point = curve[sampleID];
			SKDT::NamedPoint queryPoint(point[0], point[1], point[2]);
			Thea::BoundedSortedArray<SKDTree::Neighbor> queryResult(100);
			lineTree.kClosestElements<Thea::MetricL2>(queryPoint, queryResult, sampleRadius);
			for (int qID = 0; qID < queryResult.size(); qID++) {
				int lineSampleID = (int)lineTree.getElements()[queryResult[qID].getIndex()].id;
				allLinePointFlag[lineSampleID] = true;
			}
		}
	}

	// select lines by thresholding votes

	vector<vec2i> lineVotes(mNumLines, vec2i(0, 0)); // (aligned votes, unaligned votes) : # of lines

	for (int pointID = 0; pointID < (int)allLinePoints.size(); pointID++) {
		int lineID = allLinePointSource[pointID];
		if (allLinePointFlag[pointID]) {
			lineVotes[lineID][0]++;
		} else {
			lineVotes[lineID][1]++;
		}
	}

	mAlignedLines.clear();
	for (int lineID = 0; lineID < mNumLines; lineID++) {
		vec2i votes = lineVotes[lineID];
		if (votes[0] >(votes[0] + votes[1]) * alignThreshold) {
			mAlignedLines.push_back((*mpLines)[lineID]);
		}
	}


	return true;
}

bool CurveLineVoting::output(vector<vector<vec3>> &alignedLines) {

	alignedLines = mAlignedLines;

	return true;
}

bool CurveLineVoting::visualize(string fileName) {

	if (!CurveUtil::visualizeCurves(fileName, mAlignedLines, 0)) return false;

	return true;
}