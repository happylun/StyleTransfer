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

#include "CurveRidgeValley.h"

#include <iostream>
#include <fstream>
#include <map>

#include "Mesh/MeshUtil.h"
#include "Curve/CurveUtil.h"

#include "Data/StyleSynthesisConfig.h"

using namespace std;
using namespace StyleSynthesis;

CurveRidgeValley::CurveRidgeValley(TTriangleMesh &mesh) {

	mMesh = mesh;
}

CurveRidgeValley::~CurveRidgeValley() {

}

bool CurveRidgeValley::extractCurve() {

	if (!computeCurvatures()) return false;
	if (!extractZeroCrossings()) return false;
	if (!chainCurves()) return false;
	if (!snapCurves()) return false;

	return true;
}

bool CurveRidgeValley::computeCurvatures() {
	/*
	// subdivide mesh

	int numIterations = 5;
	double threshold = 0.01;

	vec3 bbMin, bbMax;
	if (!MeshUtil::computeAABB(mMesh, bbMin, bbMax)) return false;
	double minSubdivLen = (bbMax - bbMin).length() * threshold;

	for (int iterID = 0; iterID < numIterations; iterID++) {
		vector<int> faceIndices;
		if (!MeshUtil::subdivideMeshMidPoint(mMesh, mMesh, faceIndices, minSubdivLen)) return false;
		cout << "Iteration " << (iterID + 1) << ": # faces = " << mMesh.indices.size() << endl;
	}
	if (!MeshUtil::recomputeNormals(mMesh)) return false;

	if (!MeshUtil::saveMesh("submesh.ply", mMesh)) return false;
	*/

	// compute curvature

	if (StyleSynthesisConfig::mCurve_SmoothNormal) {
		if (!FeatureMeshCurvature::computeCurvature(mMesh, mCurvatures)) return false;
		if (!FeatureMeshCurvature::smoothNormal(mMesh, mCurvatures)) return false;
	}

	if (!FeatureMeshCurvature::computeCurvature(mMesh, mCurvatures)) return false;
	if (StyleSynthesisConfig::mCurve_SmoothCurvature) {
		if (!FeatureMeshCurvature::smoothCurvature(mMesh, mCurvatures)) return false;
	}

	if (!FeatureMeshCurvature::computeCurvatureDerivative(mMesh, mCurvatures, mCurvatureDerivatives)) return false;
	if (StyleSynthesisConfig::mCurve_SmoothCurvatureDerivative) {
		if (!FeatureMeshCurvature::smoothCurvatureDerivative(mMesh, mCurvatures, mCurvatureDerivatives)) return false;
	}

	return true;
}

bool CurveRidgeValley::extractZeroCrossings() {

	// compute vertex info

	struct TVertexInfo {
		vec3 P; // position
		vec3 TMax; // curvature direction towards increasing curvature value
		vec3 TMin;
		float KMax; // curvature value
		float KMin;
		float EMax; // curvature derivative along curvature direction
		float EMin;
	};

	int numVertices = mMesh.amount;
	int numFaces = (int)mMesh.indices.size();

	vector<TVertexInfo> vertexInfo(numVertices);
#pragma omp parallel for
	for (int vertID = 0; vertID < numVertices; vertID++) {
		TVertexInfo &info = vertexInfo[vertID];
		info.P = mMesh.positions[vertID];
		info.KMax = mCurvatures[vertID].k1;
		info.KMin = mCurvatures[vertID].k2;
		info.EMax = mCurvatureDerivatives[vertID][0];
		info.EMin = mCurvatureDerivatives[vertID][3];
		info.TMax = mCurvatures[vertID].d1 * info.EMax;
		info.TMin = mCurvatures[vertID].d2 * info.EMin;
	}

	// extract all edges

	int numEdges = 0;

	map<vec2i, int> edgeMap; // edge point ID => edge ID
	for (vec3i idx : mMesh.indices) {
		for (int k = 0; k < 3; k++) {
			vec2i edgeKey = makeKey(idx[k], idx[(k + 1) % 3]);
			if (edgeMap.find(edgeKey) == edgeMap.end()) {
				edgeMap[edgeKey] = numEdges;
				numEdges++;
			}
		}
	}
	vector<vec2i> edgeList(numEdges);
	for (auto &it : edgeMap) {
		edgeList[it.second] = it.first;
	}

	// compute zero crossing points along edges

	vector<int> edgePointType(numEdges); // -1: no point; 0: ridge; 1: valley
	vector<vec3> edgePoints(numEdges);
	vector<double> edgePointAlpha(numEdges);
	vector<float> edgePointThickness(numEdges);

	//vector<int> meshPointType(numVertices, 0); // 0: no type; 1: ridge; 2: valley
	//vector<float> meshPointThickness(numVertices, 0.0f);

#pragma omp parallel for
	for (int edgeID = 0; edgeID < numEdges; edgeID++) {

		vec2i edgeIdx = edgeList[edgeID];
		TVertexInfo &info1 = vertexInfo[edgeIdx[0]];
		TVertexInfo &info2 = vertexInfo[edgeIdx[1]];

		float alpha = -1; // it won't find ridge point and valley at the same time

		// check ridge point
		if (info1.KMax > fabs(info1.KMin) && // ridge
			info2.KMax > fabs(info2.KMin) && // ridge
			cml::dot(info1.TMax, info2.TMax) <= 0 && // zero crosing
			(cml::dot(info2.P - info1.P, info1.TMax) >= 0 || // maximum
			cml::dot(info1.P - info2.P, info2.TMax) >= 0) // maximum
			)
		{
			alpha = fabs(info2.EMax) / (fabs(info1.EMax) + fabs(info2.EMax));
			edgePointType[edgeID] = 0;
		}

		// check valley point
		if (info1.KMax < -fabs(info1.KMin) && // valley
			info2.KMax < -fabs(info2.KMin) && // valley
			cml::dot(info1.TMax, info2.TMax) <= 0 && // zero crosing
			(cml::dot(info2.P - info1.P, info1.TMax) <= 0 || // minimum
			cml::dot(info1.P - info2.P, info2.TMax) <= 0) // minimum
			)
		{
			alpha = fabs(info2.EMax) / (fabs(info1.EMax) + fabs(info2.EMax));
			edgePointType[edgeID] = 1;
		}

		// interpolate
		if (alpha >= 0.0f && alpha <= 1.0f) {
			vec3 zeroP = info1.P * alpha + info2.P * (1 - alpha);
			float zeroKMax = info1.KMax * alpha + info2.KMax * (1 - alpha);
			edgePoints[edgeID] = zeroP;
			edgePointAlpha[edgeID] = alpha;
			edgePointThickness[edgeID] = fabs(zeroKMax);

			//int snapID = alpha > 0.5f ? edgeIdx[0] : edgeIdx[1];
			//int snapType = edgePointType[edgeID] == 0 ? 1 : 2;
			//float snapThickness = fabs(vertexInfo[snapID].KMax);
			//meshPointType[snapID] = meshPointType[snapID] | snapType;
			//meshPointThickness[snapID] = max(meshPointThickness[snapID], snapThickness);
		} else {
			edgePointType[edgeID] = -1;
		}
	}

	///////////////// use zero crossing points /////////////////
	
	// insert edge points as initial curve points

	mCurvePoints.clear();
	mCurvePointSource.clear();
	mCurvePointAlpha.clear();
	mCurvePointThickness.clear();
	vector<int> edgePointMap(numEdges, -1); // curve point ID : # of edges
	for (int edgeID = 0; edgeID < numEdges; edgeID++) {
		if (edgePointType[edgeID] >= 0) {
			edgePointMap[edgeID] = (int)mCurvePoints.size();
			mCurvePoints.push_back(edgePoints[edgeID]);
			mCurvePointSource.push_back(edgeList[edgeID]);
			mCurvePointAlpha.push_back(edgePointAlpha[edgeID]);
			mCurvePointThickness.push_back(edgePointThickness[edgeID]);
		}
	}
	mCurvePointGraph.assign(mCurvePoints.size(), set<int>());

	// connect points within triangles

	for (int faceID = 0; faceID < numFaces; faceID++) {
		vec3i faceIdx = mMesh.indices[faceID];

		// find edges

		vec3i edgeIdx;
		for (int k = 0; k < 3; k++) {
			vec2i edgeKey = makeKey(faceIdx[k], faceIdx[(k + 1) % 3]);
			edgeIdx[k] = edgeMap[edgeKey];
		}

		// find points to connect

		int ridgePointCount = 0;
		int valleyPointCount = 0;
		for (int k = 0; k < 3; k++) {
			int pointType = edgePointType[edgeIdx[k]];
			if (pointType == 0) ridgePointCount++;
			else if (pointType == 1) valleyPointCount++;
		}

		if (ridgePointCount == 2 || valleyPointCount == 2) {

			// connect two points

			int connectPointType = ridgePointCount > valleyPointCount ? 0 : 1;
			vector<int> connectPointID(0);
			for (int k = 0; k < 3; k++) {
				if (edgePointType[edgeIdx[k]] == connectPointType) {
					int edgeID = edgeIdx[k];
					connectPointID.push_back(edgePointMap[edgeID]);
				}
			}
			mCurvePointGraph[connectPointID[0]].insert(connectPointID[1]);
			mCurvePointGraph[connectPointID[1]].insert(connectPointID[0]);

		} else if (ridgePointCount == 3 || valleyPointCount == 3) {

			// connect three points to center point

			vector<int> connectPointID(3);
			vector<vec3> connectPoints(3);
			for (int k = 0; k < 3; k++) {
				int edgeID = edgeIdx[k];
				connectPoints[k] = edgePoints[edgeID];
				connectPointID[k] = edgePointMap[edgeID];
			}

			vec3 centerPoint = (connectPoints[0] + connectPoints[1] + connectPoints[2]) / 3;
			float centerPointThickness = 0;
			for (int k = 0; k < 3; k++) {
				centerPointThickness += mCurvePointThickness[connectPointID[k]];
			}
			centerPointThickness /= 3;

			int centerPointID = (int)mCurvePoints.size();
			mCurvePoints.push_back(centerPoint);
			mCurvePointSource.push_back(vec2i(-1, -1));
			mCurvePointAlpha.push_back(0.0);
			mCurvePointThickness.push_back(centerPointThickness);
			mCurvePointGraph.push_back(set<int>(connectPointID.begin(), connectPointID.end()));

			for (int k = 0; k < 3; k++) {
				mCurvePointGraph[connectPointID[k]].insert(centerPointID);
			}
		}
	}	

	/*
	///////////////// use vertex points /////////////////

	// insert mesh points as initial curve points

	mCurvePoints.resize(numVertices);
	mCurvePointThickness.resize(numVertices);
	for (int vertID = 0; vertID < numVertices; vertID++) {		
		mCurvePoints[vertID] = mMesh.positions[vertID];
		mCurvePointThickness[vertID] = meshPointThickness[vertID];
	}
	mCurvePointGraph.assign(mCurvePoints.size(), set<int>());

	// connect points within triangles

	for (vec2i edge : edgeList) {
		if (meshPointType[edge[0]] & meshPointType[edge[1]]) {
			mCurvePointGraph[edge[0]].insert(edge[1]);
			mCurvePointGraph[edge[1]].insert(edge[0]);
		}
	}
	*/

	return true;
}

bool CurveRidgeValley::chainCurves() {

	vec3 bbMin, bbMax;
	if (!MeshUtil::computeAABB(mMesh, bbMin, bbMax)) return false;
	float bbLen = (bbMax - bbMin).length();

	float strengthThreshold = (float)StyleSynthesisConfig::mCurve_RidgeValleyStrength;
	float lengthThreshold = bbLen * (float)StyleSynthesisConfig::mCurve_RidgeValleyLength;

	int numCurvePoints = (int)mCurvePoints.size();

	// determine start point checking order (end point first then inner point (for looping curve))

	vector<int> endPointSet;
	vector<int> innerPointSet;
	for (int pointID = 0; pointID < numCurvePoints; pointID++) {
		int numNeighbors = (int)mCurvePointGraph[pointID].size();
		if (numNeighbors == 2) innerPointSet.push_back(pointID);
		else if (numNeighbors > 0) endPointSet.push_back(pointID);
	}
	vector<int> checkPointList(0);
	checkPointList.insert(checkPointList.end(), endPointSet.begin(), endPointSet.end());
	checkPointList.insert(checkPointList.end(), innerPointSet.begin(), innerPointSet.end());

	// start chaining

	mChainedCurvePoints.clear();
	set<vec2i> checkedEdgeSet;

	for (int pointID : checkPointList) {
		for (int neighborID : mCurvePointGraph[pointID]) {

			vector<int> chain(1, pointID);
			float chainLength = 0;
			float chainStrength = 0;

			vec2i currentSegment(pointID, neighborID);

			while (checkedEdgeSet.find(makeKey(currentSegment)) == checkedEdgeSet.end()) {

				// add current segment

				checkedEdgeSet.insert(makeKey(currentSegment));
				int lastPointID = currentSegment[0];
				int currentPointID = currentSegment[1];
				vec3 lastPoint = mCurvePoints[lastPointID];
				vec3 currentPoint = mCurvePoints[currentPointID];
				float lastThickness = mCurvePointThickness[lastPointID];
				float currentThickness = mCurvePointThickness[currentPointID];
				float segLength = (currentPoint - lastPoint).length();
				float segStrength = (currentThickness + lastThickness) / 2 * segLength;
				chain.push_back(currentPointID);
				chainLength += segLength;
				chainStrength += segStrength;

				// get next segment

				auto &nextNbSet = mCurvePointGraph[currentPointID];
				if ((int)nextNbSet.size() != 2) break;
				int nextPointID = -1;
				for (int id : nextNbSet) {
					if (id != lastPointID) nextPointID = id;
				}
				currentSegment = vec2i(currentPointID, nextPointID);
			}

			// check chain

			chainStrength /= chainLength;
			if (chainLength > lengthThreshold && chainStrength > strengthThreshold) {
				mChainedCurvePoints.push_back(chain);
			}
		}
	}

	return true;
}

bool CurveRidgeValley::snapCurves() {

	// snap curve points to vertices
	/* // HACK: skip snapping
	int numCurvePoints = (int)mCurvePoints.size();
	for (int pointID = 0; pointID < numCurvePoints; pointID++) {
		vec2i edge = mCurvePointSource[pointID];
		if (edge[0] < 0) continue;
		double alpha = mCurvePointAlpha[pointID];
		mCurvePoints[pointID] = alpha > 0.5 ? mMesh.positions[edge[0]] : mMesh.positions[edge[1]];
	}
	*/

	// update chained curves

	int numCurves = (int)mChainedCurvePoints.size();
	mChainedCurves.resize(numCurves);
	for (int curveID = 0; curveID < numCurves; curveID++) {
		auto &curvePoints = mChainedCurvePoints[curveID];
		auto &curve = mChainedCurves[curveID];
		curve.resize(curvePoints.size());
		for (int id = 0; id < (int)curvePoints.size(); id++) {
			curve[id] = mCurvePoints[curvePoints[id]];
		}
	}

	return true;
}

bool CurveRidgeValley::output(vector<vector<vec3>> &curves) {

	curves = mChainedCurves;

	return true;
}

bool CurveRidgeValley::visualize(string fileName) {

	if (!CurveUtil::visualizeCurves(fileName, mChainedCurves, 2)) return false;

	return true;
}
