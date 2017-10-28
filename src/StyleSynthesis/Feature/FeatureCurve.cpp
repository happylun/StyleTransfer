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

#include "FeatureCurve.h"

#include <fstream>
#include <iostream>

#include "Sample/SampleUtil.h"
#include "Curve/CurveUtil.h"

#include "Utility/PlyExporter.h"

using namespace StyleSynthesis;

FeatureCurve::FeatureCurve() {
}

FeatureCurve::~FeatureCurve() {
}

bool FeatureCurve::loadCurves(string curveFolder) {

	string curveRVName = curveFolder + "data-snap-rv.txt";
	string curveBName = curveFolder + "data-snap-b.txt";
	string curveCName = curveFolder + "data-snap-c.txt";

	if (!CurveUtil::loadCurves(curveRVName, mCurveRV)) return false;

	if (!CurveUtil::loadCurves(curveBName, mCurveB)) return false;

	ifstream curveCFile(curveCName, ios::binary);
	int numViews;
	curveCFile.read((char*)&numViews, sizeof(numViews));
	mCurveC.resize(numViews);
	for (int viewID = 0; viewID < numViews; viewID++) {
		auto &viewCurves = mCurveC[viewID];
		if (!CurveUtil::loadCurves(curveCFile, viewCurves)) return false;
	}
	curveCFile.close();

	return true;
}

bool FeatureCurve::extractPointClouds(double sampleRadius) {

	// 0: ridges & valleys
	// 1: boundaries
	// 2: contours from all view points
	// 3: all curves


	// combine all contours

	vector<vector<vec3>> allContours(0);
	for (auto &contour : mCurveC) {
		allContours.insert(allContours.end(), contour.begin(), contour.end());
	}

	// combine all curves

	vector<vector<vec3>> allCurves(0);
	allCurves.insert(allCurves.end(), mCurveRV.begin(), mCurveRV.end());
	allCurves.insert(allCurves.end(), mCurveB.begin(), mCurveB.end());
	allCurves.insert(allCurves.end(), allContours.begin(), allContours.end());

	vector< vector<vector<vec3>>* > allCurveTypes(0);
	allCurveTypes.push_back(&mCurveRV);
	allCurveTypes.push_back(&mCurveB);
	allCurveTypes.push_back(&allContours);
	allCurveTypes.push_back(&allCurves);

	int numCurveTypes = (int)allCurveTypes.size();
	mPointClouds.resize(numCurveTypes);
	for (int typeID = 0; typeID < numCurveTypes; typeID++) {

		vector<vector<vec3>> *curves = allCurveTypes[typeID];
		if (!CurveUtil::removeDuplicateLines(*curves)) return false;
		int numCurves = (int)curves->size();

		// get sample points on curves

		vector<vector<vec3>> sampleP, sampleN; // position & direction
		if (!CurveUtil::sampleLines(*curves, sampleP, (float)sampleRadius)) return false;
		sampleN.resize(numCurves);
		for (int curveID = 0; curveID < numCurves; curveID++) {
			if (!CurveUtil::computeCurveDirections(sampleP[curveID], sampleN[curveID])) return false;
		}

		// export point cloud

		auto &pc = mPointClouds[typeID];
		pc.positions.clear();
		pc.normals.clear();
		for (int curveID = 0; curveID < numCurves; curveID++) {
			for (vec3 &v : sampleP[curveID]) {
				pc.positions.push_back(v);
			}
			for (vec3 &n : sampleN[curveID]) {
				pc.normals.push_back(n);
			}
		}
		pc.amount = (int)pc.positions.size();
	}

	return true;
}

bool FeatureCurve::visualize(string fileName) {

	int numPointClouds = (int)mPointClouds.size();

	vec3 allBBMin(FLT_MAX, FLT_MAX, FLT_MAX);
	vec3 allBBMax(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	for (auto &pc : mPointClouds) {
		if (pc.amount == 0) continue;
		vec3 bbMin, bbMax;
		if (!SampleUtil::computeAABB(pc, bbMin, bbMax)) return false;
		allBBMin.minimize(bbMin);
		allBBMax.maximize(bbMax);
	}

	float spacing = (allBBMax[0] - allBBMin[0]) * 1.5f;

	PlyExporter pe;
	for (int pcID = 0; pcID < numPointClouds; pcID++) {
		auto &pc = mPointClouds[pcID];
		vec3 offset(spacing*pcID, 0.0f, 0.0f);
		if (!pe.addPoint(&pc.positions, &pc.normals, offset)) return false;
	}
	if (!pe.output(fileName)) return false;

	return true;
}

bool FeatureCurve::loadFeature(string fileName) {

	ifstream file(fileName, ios::binary);

	int numTypes;
	file.read((char*)&numTypes, sizeof(numTypes));
	mPointClouds.resize(numTypes);
	for (auto &pc : mPointClouds) {
		int numPoints;
		file.read((char*)&numPoints, sizeof(numPoints));
		pc.positions.resize(numPoints);
		pc.normals.resize(numPoints);
		for (int pointID = 0; pointID < numPoints; pointID++) {
			vec3 v, n;
			file.read((char*)v.data(), sizeof(v));
			file.read((char*)n.data(), sizeof(n));
			pc.positions[pointID] = v;
			pc.normals[pointID] = n;
		}
		pc.amount = (int)pc.positions.size();
	}

	file.close();

	return true;
}

bool FeatureCurve::saveFeature(string fileName) {

	ofstream file(fileName, ios::binary);

	int numTypes = (int)mPointClouds.size();
	file.write((char*)&numTypes, sizeof(numTypes));
	for (auto &pc : mPointClouds) {
		int numPoints = pc.amount;
		file.write((char*)&numPoints, sizeof(numPoints));
		for (int pointID = 0; pointID < numPoints; pointID++) {
			vec3 &v = pc.positions[pointID];
			vec3 &n = pc.normals[pointID];
			file.write((char*)v.data(), sizeof(v));
			file.write((char*)n.data(), sizeof(n));
		}
	}

	file.close();

	return true;
}