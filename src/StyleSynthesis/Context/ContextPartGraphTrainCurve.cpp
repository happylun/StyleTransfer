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

#include "ContextPartGraphTrainCurve.h"

#include <iostream>
#include <fstream>

#include "Match/MatchRigidICP.h"

#include "Mesh/MeshUtil.h"
#include "Curve/CurveUtil.h"
#include "Segment/SegmentUtil.h"

#include "Utility/PlyExporter.h"

#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

ContextPartGraphTrainCurve::ContextPartGraphTrainCurve() {
}

ContextPartGraphTrainCurve::~ContextPartGraphTrainCurve() {
}

bool ContextPartGraphTrainCurve::loadSceneCurves(vector<string> &sceneCurveFolders, vector<vec3> &viewPoints) {

	mCurveViewPoints = viewPoints;

	int numMeshes = (int)sceneCurveFolders.size();

	mSceneCurveIndices.resize(numMeshes);
	mSceneCurves.resize(numMeshes);

	for (int meshID = 0; meshID < numMeshes; meshID++) {

		vector<vector<vec3>> &curves = mSceneCurves[meshID];
		vector<vec2i> &curveIndices = mSceneCurveIndices[meshID];
		curves.clear();
		curveIndices.clear();

		string curveFolder = sceneCurveFolders[meshID];

		string curveRVName = curveFolder + "data-snap-rv.txt";
		string curveBName = curveFolder + "data-snap-b.txt";
		string curveCName = curveFolder + "data-snap-c.txt";

		vector<vector<vec3>> curveRV;
		if (!CurveUtil::loadCurves(curveRVName, curveRV)) return false;
		curveIndices.push_back(vec2i((int)curves.size(), (int)(curves.size() + curveRV.size())));
		curves.insert(curves.end(), curveRV.begin(), curveRV.end());

		vector<vector<vec3>> curveB;
		if (!CurveUtil::loadCurves(curveBName, curveB)) return false;
		curveIndices.push_back(vec2i((int)curves.size(), (int)(curves.size() + curveB.size())));
		curves.insert(curves.end(), curveB.begin(), curveB.end());

		int numViewPoints;
		ifstream curveCFile(curveCName, ios::binary);
		if (!curveCFile.is_open()) return error("cannot open file " + curveCName);
		curveCFile.read((char*)&numViewPoints, sizeof(numViewPoints));
		if (numViewPoints != (int)mCurveViewPoints.size()) return error("incompatible contour file " + curveCName);
		for (int viewID = 0; viewID < numViewPoints; viewID++) {
			vector<vector<vec3>> curveC;
			if (!CurveUtil::loadCurves(curveCFile, curveC)) return false;
			curveIndices.push_back(vec2i((int)curves.size(), (int)(curves.size() + curveC.size())));
			curves.insert(curves.end(), curveC.begin(), curveC.end());
		}
		curveCFile.close();
	}

	return true;
}

bool ContextPartGraphTrainCurve::process() {

	if (!initCurves()) return false;
	if (!extractCurvePairs()) return false;

	return true;
}

bool ContextPartGraphTrainCurve::initCurves() {

	// filter out uninteresting curves

	double minCurveLength = StyleSynthesisConfig::mDeform_MinimumMatchedCurveLength;

	vector<vec2i> curveList(0); // (meshID, curveID) : # of curves

	int numMeshes = (int)mSceneCurves.size();
	mSceneCurveFlags.resize(numMeshes);

	for (int meshID = 0; meshID < numMeshes; meshID++) {
		int numCurves = (int)mSceneCurves[meshID].size();
		mSceneCurveFlags[meshID].resize(numCurves, true);

		for (int curveID = 0; curveID < numCurves; curveID++) {
			curveList.push_back(vec2i(meshID, curveID));
		}
	}
	int numCurves = (int)curveList.size();

	int totalCount = 0;
	int validCount = 0;
#pragma omp parallel for
	for (int id = 0; id < numCurves; id++) {
#pragma omp atomic
		totalCount++;
#pragma omp critical
		if (totalCount % 100 == 0) cout << "\rPre-processing curve " << totalCount << " / " << numCurves << "      ";

		int meshID = curveList[id][0];
		int curveID = curveList[id][1];
		vector<vec3> &curve = mSceneCurves[meshID][curveID];

		// check curve length
		double curveLen;
		if (!CurveUtil::computeCurveLength(curve, curveLen)) error("curve length");
		if (curveLen < minCurveLength) {
			mSceneCurveFlags[meshID][curveID] = false;
			continue;
		}

		// check straight curve
		bool straightFlag;
		if (!CurveUtil::checkStraightCurve(curve, straightFlag)) error("straight curve");
		if (straightFlag) {
			mSceneCurveFlags[meshID][curveID] = false;
			continue;
		}

#pragma omp atomic
		validCount++;
	}
	cout << endl;
	cout << "Filtered " << (totalCount - validCount) << " out of " << totalCount << " curves" << endl;

	return true;
}

bool ContextPartGraphTrainCurve::extractCurvePairs() {

	int numMeshes = (int)mSceneCurves.size();

	// extract mesh pairs

	mMeshPairs.clear();
	map<vec2i, int> meshPairMap;
	for (int i = 0; i < numMeshes - 1; i++) {
		for (int j = i + 1; j < numMeshes; j++) {
			vec2i key(i, j);
			meshPairMap[key] = (int)mMeshPairs.size();
			mMeshPairs.push_back(key);
		}
	}
	int numMeshPairs = (int)mMeshPairs.size();

	// extract potential curve pairs

	vector<vec4i> curvePairs; // (srcMeshID, srcCurveID, tgtMeshID, tgtCurveID) : # of potential curve pairs

	for (int pairID = 0; pairID < numMeshPairs; pairID++) {
		int srcMeshID = mMeshPairs[pairID][0];
		int tgtMeshID = mMeshPairs[pairID][1];
		
		// match non-contour to non-contour
		for (int srcCurveID = mSceneCurveIndices[srcMeshID][0][0]; srcCurveID < mSceneCurveIndices[srcMeshID][1][1]; srcCurveID++) {
			if (!mSceneCurveFlags[srcMeshID][srcCurveID]) continue;
			for (int tgtCurveID = mSceneCurveIndices[tgtMeshID][0][0]; tgtCurveID < mSceneCurveIndices[tgtMeshID][1][1]; tgtCurveID++) {
				if (!mSceneCurveFlags[tgtMeshID][tgtCurveID]) continue;
				curvePairs.push_back(vec4i(srcMeshID, srcCurveID, tgtMeshID, tgtCurveID));
			}
		}

		// match contour to contour
		for (int srcCurveID = mSceneCurveIndices[srcMeshID][2][0]; srcCurveID < (int)mSceneCurves[srcMeshID].size(); srcCurveID++) {
			if (!mSceneCurveFlags[srcMeshID][srcCurveID]) continue;
			for (int tgtCurveID = mSceneCurveIndices[tgtMeshID][2][0]; tgtCurveID < (int)mSceneCurves[tgtMeshID].size(); tgtCurveID++) {
				if (!mSceneCurveFlags[tgtMeshID][tgtCurveID]) continue;
				curvePairs.push_back(vec4i(srcMeshID, srcCurveID, tgtMeshID, tgtCurveID));
			}
		}
	}

	// check each curve pairs

	int numCurvePairs = (int)curvePairs.size();
	vector<bool> curveFlags(numCurvePairs, false);
	int count = 0;
#pragma omp parallel for
	for (int pairID = 0; pairID < numCurvePairs; pairID++) {
#pragma omp atomic
		count++;
#pragma omp critical
		if (count % 100000 == 0) cout << "\rChecking curve " << count << " / " << numCurvePairs << "      ";

		vec4i idx = curvePairs[pairID];
		vector<vec3> &srcCurve = mSceneCurves[idx[0]][idx[1]];
		vector<vec3> &tgtCurve = mSceneCurves[idx[2]][idx[3]];
		if (checkCurvePair(srcCurve, tgtCurve)) {
			curveFlags[pairID] = true;
		}
	}
	cout << endl;

	// organize matched curve pairs

	mCurvePairs.resize(numMeshPairs);
	int numMatchedCurvePairs = 0;
	for (int pairID = 0; pairID < numCurvePairs; pairID++) {
		if (!curveFlags[pairID]) continue;
		vec4i idx = curvePairs[pairID];
		int meshPairID = meshPairMap[vec2i(idx[0], idx[2])];
		mCurvePairs[meshPairID].push_back(vec2i(idx[1], idx[3]));
		numMatchedCurvePairs++;
	}
	cout << "Extracted " << numMatchedCurvePairs << " matched curve pairs" << endl;

	return true;
}

bool ContextPartGraphTrainCurve::checkCurvePair(vector<vec3> &curve1, vector<vec3> &curve2) {

	// prune by angle
	if (true) {
		double maxCurveAngle = 10.0; // UNDONE: param maximum angle between matched curve direction

		vec3d dir1 = curve1.front() - curve1.back();
		vec3d dir2 = curve2.front() - curve2.back();
		double dirLen1 = dir1.length();
		double dirLen2 = dir2.length();
		if (dirLen1 && dirLen2) {
			dir1 /= dirLen1;
			dir2 /= dirLen2;
			double angle = cml::deg(acos(cml::clamp(fabs(cml::dot(dir1, dir2)), 0.0, 1.0)));
			if (angle > maxCurveAngle) return false; // pruned by angle
		} else if (dirLen1 || dirLen2) {
			return false; // one curve is a loop while the othe one is not
		} else {
			// both curves are loops
		}
	}

	// prune by length
	if (true) {
		double maxLenRatio = 1.5; // UNDONE: param maximum length ratio between matched curves

		double curveLen1, curveLen2;
		if (!CurveUtil::computeCurveLength(curve1, curveLen1)) return false;
		if (!CurveUtil::computeCurveLength(curve2, curveLen2)) return false;
		double ratio = curveLen1 / curveLen2;
		if (ratio > maxLenRatio || ratio < 1 / maxLenRatio) return false;
	}

	// check ICP
	if (true) {
		double maxError = cml::sqr(StyleSynthesisConfig::mCurve_SamplingRadius * 1.0); // UNDONE: param max match error
		double sampleRadius = StyleSynthesisConfig::mCurve_SamplingRadius;

		// build matrices
		vector<vec3> curvePos1, curvePos2;
		vector<vec3> curveDir1, curveDir2;
		if (!CurveUtil::sampleLine(curve1, curvePos1, (float)sampleRadius)) return false;
		if (!CurveUtil::sampleLine(curve2, curvePos2, (float)sampleRadius)) return false;
		if (!CurveUtil::computeCurveDirections(curvePos1, curveDir1)) return false;
		if (!CurveUtil::computeCurveDirections(curvePos2, curveDir2)) return false;
		Eigen::Matrix3Xd matP1, matN1, matP2, matN2;
		if (!CurveUtil::buildMatrix(curvePos1, matP1)) return false;
		if (!CurveUtil::buildMatrix(curvePos2, matP2)) return false;
		if (!CurveUtil::buildMatrix(curveDir1, matN1)) return false;
		if (!CurveUtil::buildMatrix(curveDir2, matN2)) return false;

		// run ICP
		Eigen::Affine3d xform;
		if (!MatchRigidICP::run(3, matP1, matN1, matP2, matN2, xform)) return false;

		// check error
		double error1;
		Eigen::Matrix3Xd matXP1 = xform * matP1;
		if (!MatchRigidICP::error(matXP1, matP2, error1)) return false;
		if (error1 > maxError) return false;
		double error2;
		Eigen::Matrix3Xd matXP2 = xform.inverse() * matP2;
		if (!MatchRigidICP::error(matXP2, matP1, error2)) return false;
		if (error2 > maxError) return false;
	}

	return true;
}

bool ContextPartGraphTrainCurve::visualizeCurvePairs(string fileName, vector<TTriangleMesh> &meshes) {

	int numVisCurves = 10;
	int numMeshPairs = (int)mMeshPairs.size();

	double tubeRadius = StyleSynthesisConfig::mCurve_VisualizationTubeRadius;

	int numMeshes = (int)meshes.size();
	vector<vec3> meshBBSizes(numMeshes);
	for (int meshID = 0; meshID < numMeshes; meshID++) {
		vec3 bbMin, bbMax;
		if (!MeshUtil::computeAABB(meshes[meshID], bbMin, bbMax)) return false;
		meshBBSizes[meshID] = bbMax - bbMin;
	}

	PlyExporter pe;
	vec3 offset1(0.0f, 0.0f, 0.0f);
	vec3 offset2 = offset1;
	float lastOffsetH = 0.0f;
	for (int meshPairID = 0; meshPairID < numMeshPairs; meshPairID++) {

		int meshID1 = mMeshPairs[meshPairID][0];
		int meshID2 = mMeshPairs[meshPairID][1];
		TTriangleMesh &mesh1 = meshes[meshID1];
		TTriangleMesh &mesh2 = meshes[meshID2];

		float offsetH = max(meshBBSizes[meshID1][0], meshBBSizes[meshID2][0]) * 1.5f;
		float offsetV = max(meshBBSizes[meshID1][1], meshBBSizes[meshID2][1]) * 1.5f;
		offset1 += vec3((lastOffsetH+offsetH)/2, 0.0f, 0.0f);
		offset2 = offset1 + vec3(0.0f, -offsetV, 0.0f);
		lastOffsetH = offsetH;

		if (!pe.addMesh(&mesh1.positions, &mesh1.normals, &mesh1.indices, vec3i(127, 127, 127), offset1)) return false;
		if (!pe.addMesh(&mesh2.positions, &mesh2.normals, &mesh2.indices, vec3i(127, 127, 127), offset2)) return false;

		int numCurvePairs = min(numVisCurves, (int)mCurvePairs[meshPairID].size());
		for (int curvePairID = 0; curvePairID < numCurvePairs; curvePairID++) {
			int curveID1 = mCurvePairs[meshPairID][curvePairID][0];
			int curveID2 = mCurvePairs[meshPairID][curvePairID][1];

			vector<vec3> &curve1 = mSceneCurves[meshID1][curveID1];
			vector<vec3> &curve2 = mSceneCurves[meshID2][curveID2];

			TTriangleMesh tube1, tube2;
			if (!CurveUtil::makeTube(curve1, tube1, (float)tubeRadius)) return false;
			if (!CurveUtil::makeTube(curve2, tube2, (float)tubeRadius)) return false;

			vec3i color = SegmentUtil::colorMapping(curvePairID);
			if (!pe.addMesh(&tube1.positions, &tube1.normals, &tube1.indices, color, offset1)) return false;
			if (!pe.addMesh(&tube2.positions, &tube2.normals, &tube2.indices, color, offset2)) return false;
		}
	}

	if (!pe.output(fileName)) return false;

	return true;
}

bool ContextPartGraphTrainCurve::outputCurvePairs(string fileName) {

	ofstream file(fileName);
	if (!file.is_open()) return error("cannot write to file " + fileName);

	int numMeshPairs = (int)mMeshPairs.size();
	file << numMeshPairs << endl;
	for (int meshPairID = 0; meshPairID < numMeshPairs; meshPairID++) {
		vec2i meshPair = mMeshPairs[meshPairID];
		file << meshPair[0] << " " << meshPair[1];
		file << " " << mCurvePairs[meshPairID].size() << endl;
		for (vec2i curvePair : mCurvePairs[meshPairID]) {
			file << curvePair[0] << " " << curvePair[1] << " ";
		}
		file << endl;
	}

	file.close();

	return true;
}