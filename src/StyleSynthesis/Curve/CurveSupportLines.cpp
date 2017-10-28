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

#include "CurveSupportLines.h"

#include <iostream>
#include <fstream>
#include <map>
#include <set>

#include "Curve/CurveUtil.h"

#include "Data/StyleSynthesisConfig.h"

using namespace std;
using namespace StyleSynthesis;

CurveSupportLines::CurveSupportLines(TTriangleMesh &mesh) {

	mpMesh = &mesh;
}

CurveSupportLines::~CurveSupportLines() {

}

bool CurveSupportLines::extractLines() {

	if (!extractSupportEdges()) return false;
	if (!chainSupportEdges()) return false;

	return true;
}

bool CurveSupportLines::extractSupportEdges() {

	// find support edges (boundaries or sharing faces have different normals)

	map<vec2i, pair<vec3,bool>> allEdgeSet; // edge key => (any neighbor face normal, is boundary)
	set<vec2i> supportEdgeSet;

	double cosineEps = 1e-3;

	// process all edges

	for (vec3i faceIdx : mpMesh->indices) {
		vec3 faceP[3];
		for (int k = 0; k < 3; k++) faceP[k] = mpMesh->positions[faceIdx[k]];
		vec3 faceN = cml::normalize(cml::cross(faceP[1] - faceP[0], faceP[2] - faceP[0]));
		for (int k = 0; k < 3; k++) {
			vec2i edgeKey = makeKey(faceIdx[k], faceIdx[(k + 1) % 3]);
			auto it = allEdgeSet.find(edgeKey);
			if (it == allEdgeSet.end()) {
				allEdgeSet[edgeKey] = make_pair(faceN, true); // initially it is boundary
			} else {
				it->second.second = false; // shared faces -- not boundary
				vec3 edgeN = it->second.first;
				double cosine = fabs(cml::dot(faceN, edgeN));
				if (fabs(cosine-1.0) > cosineEps) { // different normal
					supportEdgeSet.insert(edgeKey);
				}
			}
		}
	}

	// add boundaries

	for (auto &it : allEdgeSet) {
		if (it.second.second) {
			supportEdgeSet.insert(it.first);
		}
	}

	mSupportEdges.assign(supportEdgeSet.begin(), supportEdgeSet.end());

	return true;
}

bool CurveSupportLines::chainSupportEdges() {

	double maxChainAngle = StyleSynthesisConfig::mCurve_MaximumChainingAngle;
	vector<vector<int>> chainIndices;
	vector<vector<vec3>> initialCurves;
	if (!CurveUtil::chainLines(mSupportEdges, mpMesh->positions, chainIndices, maxChainAngle)) return false;
	if (!CurveUtil::extractLines(chainIndices, mpMesh->positions, initialCurves)) return false;
	if (!CurveUtil::removeDuplicateLines(initialCurves)) return false;

	// split highly non-linear curve
	mSupportLines.clear();
	for (int lineID = 0; lineID < (int)initialCurves.size(); lineID++) {
		vector<vector<vec3>> curveSegments;
		vector<vector<vec3>> reversedCurveSegments;
		vector<vec3> &curve = initialCurves[lineID];
		vector<vec3> reversedCurve(curve.rbegin(), curve.rend());
		if (!splitCurve(curve, curveSegments)) return false;
		if (!splitCurve(reversedCurve, reversedCurveSegments)) return false;
		if (curveSegments.size() >= reversedCurveSegments.size()) {
			mSupportLines.insert(mSupportLines.end(), curveSegments.begin(), curveSegments.end());
		} else {
			mSupportLines.insert(mSupportLines.end(), reversedCurveSegments.begin(), reversedCurveSegments.end());
		}
	}

	// HACK: prune curves if there are lots of them...
	if ((int)mSupportLines.size() > 10000) {
		vector<double> lineLen(mSupportLines.size());
		double maxLen = 0;
		for (int lineID = 0; lineID < (int)mSupportLines.size(); lineID++) {
			if (!CurveUtil::computeCurveLength(mSupportLines[lineID], lineLen[lineID])) return false;
			if (lineLen[lineID] > maxLen) maxLen = lineLen[lineID];
		}
		vector<vector<vec3>> prunedLines;
		for (int lineID = 0; lineID < (int)mSupportLines.size(); lineID++) {
			if (lineLen[lineID] > maxLen * 0.2) { // UNDONE: param pruning length threshold
				prunedLines.push_back(mSupportLines[lineID]);
			}
		}
		mSupportLines.swap(prunedLines);
	}

	return true;
}

bool CurveSupportLines::splitCurve(vector<vec3> &inCurve, vector<vector<vec3>> &outCurves) {

	double curveLen;
	if (!CurveUtil::computeCurveLength(inCurve, curveLen)) return false;

	// check whether it's necessary to split this curve
	if (true) {
		if (inCurve.size() <= 2) {
			// straight line...
			outCurves.assign(1, inCurve);
			return true;
		}
		vec3d firstPoint = inCurve.front();
		vec3d lastPoint = inCurve.back();
		bool skip = true;
		for (int pointID = 1; pointID < (int)inCurve.size() - 1; pointID++) {
			vec3d midPoint = inCurve[pointID];
			if (cml::unsigned_angle(midPoint - firstPoint, lastPoint - midPoint) > cml::rad(30.0)) {
				// not quite "linear" -- split it!
				skip = false;
				break;
			}
		}
		if (skip) {
			outCurves.assign(1, inCurve);
			return true;
		}
	}

	vector<vec3> curve = inCurve;

	// determine processing direction

	if (true) {
		vec3d startPoint = curve.front();
		vec3d finishPoint = curve.back();
		bool flip = false;
		if (startPoint == finishPoint) {
			vec3d dir1 = cml::normalize(curve[1] - curve.front());
			vec3d dir2 = cml::normalize(curve[curve.size() - 2] - curve.back());
			vec3d axis(1.0, 1.0, 1.0);
			if (cml::dot(dir1, axis) < cml::dot(dir2, axis)) flip = true;
		} else {
			vec3d center(0.0, 0.0, 0.0);
			for (vec3 &point : curve) {
				center += point;
			}
			center /= (double)curve.size();
			double startLen = (center - startPoint).length();
			double finishLen = (center - finishPoint).length();
			if (fabs(startLen - finishLen) < 1e-6) {
				// probably symmetry
				if (cml::dot(finishPoint - startPoint, vec3d(1.0, 1.0, 1.0)) < 0) flip = true;
			} else if (startLen < finishLen) {
				flip = true;
			}
		}
		if (flip) reverse(curve.begin(), curve.end());
	}

	double minSegLen = curveLen * StyleSynthesisConfig::mCurve_MinimumSegmentLength;

	vector<vector<vec3>> splittedCurves(0);
	vector<vector<vec3>> splittedCurvesRev(0);

	int numPoints = (int)curve.size();
	int startID = 0;
	double accumLen = 0;
	bool revFlag = false;
	while (true) {

		double minError = DBL_MAX;
		int minSplitID = -1;
		double minSplitLen = 0;

		double segLen = 0;
		double segAccumWeight = 0;
		vec3d segAccumMidPoint(0.0, 0.0, 0.0);
		for (int splitID = startID + 1; splitID < numPoints; splitID++) {

			double pointSegLen = (double)(curve[splitID] - curve[splitID - 1]).length();
			segLen += pointSegLen;

			double weight = pointSegLen; // length weight
			segAccumWeight += weight;
			segAccumMidPoint += vec3d(curve[splitID] + curve[splitID - 1]) / 2 * weight;
			if (segLen < minSegLen && splitID + 1 < numPoints) continue;

			vec3 segMidPoint = segAccumMidPoint / segAccumWeight;
						
			vec3d lineVec = curve[splitID] - curve[startID];
			double lineVecLen = lineVec.length();
			if (lineVecLen == 0) break;
			lineVec /= lineVecLen;

			// compute error
			double currentError = 0; // weighted average square distance
			double innerAccumWeight = 0;;
			for (int innerID = startID; innerID < splitID; innerID++) {
				vec3d innerPoint = vec3d(curve[innerID] + curve[innerID + 1]) / 2;
				double lineDist = fabs(cml::dot(innerPoint - segMidPoint, lineVec));
				lineDist = cml::sqrt_safe((innerPoint - segMidPoint).length_squared() - lineDist*lineDist);

				vec3d innerSeg = curve[innerID + 1] - curve[innerID];
				double innerSegLen = innerSeg.length();
				double innerWeight = innerSegLen; // length weight
				//double innerWeight = fabs(cml::dot(innerSeg, lineVec)); // projected length weight

				currentError += cml::sqr(lineDist) * innerWeight;
				innerAccumWeight += innerWeight;
			}
			currentError /= innerAccumWeight;

			if (currentError < minError) {
				minError = currentError;
				minSplitID = splitID;
				minSplitLen = segLen;
			}
		}

		double remLen = curveLen - minSplitLen - accumLen;
		if (remLen <= minSegLen || minSplitID <= startID) {
			// remaining segment is shorter than minSegLen; merge with current segment
			minSplitID = numPoints - 1;
			minSplitLen = curveLen - accumLen;
		}

		vector<vec3> currentSegment(curve.begin() + startID, curve.begin() + minSplitID + 1);
		if (revFlag) {
			vector<vec3> revSegmeent(currentSegment.rbegin(), currentSegment.rend());
			splittedCurvesRev.push_back(revSegmeent);
		} else {
			splittedCurves.push_back(currentSegment);
		}

		startID = minSplitID;
		accumLen += minSplitLen;
		if (startID + 1 >= numPoints) break;

		// flip remaining section (so that it proceeds from two sides towards middle)
		reverse(curve.begin() + startID, curve.end());
		revFlag = !revFlag;
	}
	splittedCurves.insert(splittedCurves.end(), splittedCurvesRev.rbegin(), splittedCurvesRev.rend());

	//outCurves.swap(splittedCurves);
	//return true;

	// chain splitted curves

	outCurves.clear();
	if (!splittedCurves.empty()) {
		double maxAngle = cml::rad(StyleSynthesisConfig::mCurve_MaximumChainingAngle);
		int numSplittedCurves = (int)splittedCurves.size();
		vector<vec3> currentCurve = splittedCurves[0];
		for (int curveID = 1; curveID < numSplittedCurves; curveID++) {
			auto &nextCurve = splittedCurves[curveID];
			vec3d currentDir = currentCurve.back() - currentCurve.front();
			vec3d nextDir = nextCurve.back() - nextCurve.front();
			vec3d mergeDir = nextCurve.back() - currentCurve.front();
			double connectAngle = cml::unsigned_angle(currentDir, nextDir);
			double extendAngle = cml::unsigned_angle(currentDir, mergeDir);
			if (connectAngle < maxAngle && extendAngle < maxAngle) {
				currentCurve.insert(currentCurve.end(), nextCurve.begin() + 1, nextCurve.end());
			} else {
				outCurves.push_back(currentCurve);
				currentCurve = nextCurve;
			}
		}
		if(!currentCurve.empty()) outCurves.push_back(currentCurve);
	}

	return true;
}

bool CurveSupportLines::output(vector<vector<vec3>> &curves) {

	curves = mSupportLines;

	return true;
}

bool CurveSupportLines::visualize(string fileName) {

	if (!CurveUtil::visualizeCurves(fileName, mSupportLines, 0)) return false;

	return true;
}
