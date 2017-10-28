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

#include "CurveContour.h"

#include <iostream>
#include <fstream>
#include <map>

#include "Curve/CurveUtil.h"
#include "Mesh/MeshUtil.h"

#include "Utility/PlyExporter.h"

#include "Data/StyleSynthesisConfig.h"

using namespace std;
using namespace StyleSynthesis;

CurveContour::CurveContour(TTriangleMesh &mesh) {

	mpMesh = &mesh;
}

CurveContour::~CurveContour() {

}

bool CurveContour::extractCurve(vec3 viewPoint) {

	//if (!extractContourEdgesByZeroCrossing(viewPoint)) return false;
	if (!extractContourEdgesByFaceOrientation(viewPoint)) return false;
	if (!chainContourEdges()) return false;

	return true;
}

bool CurveContour::extractContourEdgesByZeroCrossing(vec3 viewPoint) {

	int numFaces = (int)mpMesh->indices.size();
	int numVertives = mpMesh->amount;

	// compute n dot v per vertex

	vector<float> ndotv(numVertives);
#pragma omp parallel for
	for (int vertID = 0; vertID < numVertives; vertID++) {
		vec3 p = mpMesh->positions[vertID];
		vec3 n = mpMesh->normals[vertID];
		vec3 v = cml::normalize(viewPoint - p);
		ndotv[vertID] = cml::dot(n, v);
	}

	// extract all edges

	int numEdges = 0;

	map<vec2i, int> edgeMap; // edge key => edge ID
	for (vec3i idx : mpMesh->indices) {
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

	vector<bool> edgePointFlag(numEdges, false);
	vector<vec3> edgePoints(numEdges);

#pragma omp parallel for
	for (int edgeID = 0; edgeID < numEdges; edgeID++) {

		vec2i edgeIdx = edgeList[edgeID];
		float value1 = ndotv[edgeIdx[0]];
		float value2 = ndotv[edgeIdx[1]];
		if (value1*value2 < 0) {
			vec3 point1 = mpMesh->positions[edgeIdx[0]];
			vec3 point2 = mpMesh->positions[edgeIdx[1]];
			float alpha = fabs(value2) / (fabs(value1) + fabs(value2));
			edgePoints[edgeID] = point1 * alpha + point2 * (1 - alpha);
			edgePointFlag[edgeID] = true;
		}
	}


	// insert valid edge points
	
	vector<int> pointMap(numEdges, -1); // curve point ID : # of edges
	mCurvePoints.clear();
	for (int edgeID = 0; edgeID < numEdges; edgeID++) {
		if (edgePointFlag[edgeID]) {
			pointMap[edgeID] = (int)mCurvePoints.size();
			mCurvePoints.push_back(edgePoints[edgeID]);
		}
	}
	mCurvePointGraph.assign(mCurvePoints.size(), set<int>());

	// build curve point graph

	for (int faceID = 0; faceID < numFaces; faceID++) {

		vec3i faceIdx = mpMesh->indices[faceID];
		vec3 faceP[3];
		for (int k = 0; k < 3; k++) faceP[k] = mpMesh->positions[faceIdx[k]];
		vec3 faceN = cml::normalize(cml::cross(faceP[1] - faceP[0], faceP[2] - faceP[0]));
		//if (cml::dot(faceN, viewPoint) < 0) continue; // skip back faces

		// find edges

		vector<int> zeroEdges(0);
		for (int k = 0; k < 3; k++) {
			vec2i edgeKey = makeKey(faceIdx[k], faceIdx[(k + 1) % 3]);
			int edgeID = edgeMap[edgeKey];
			if (edgePointFlag[edgeID]) zeroEdges.push_back(edgeID);
		}

		// connect zero-crossing points

		if ((int)zeroEdges.size() == 2) {
			int id1 = pointMap[zeroEdges[0]];
			int id2 = pointMap[zeroEdges[1]];
			mCurvePointGraph[id1].insert(id2);
			mCurvePointGraph[id2].insert(id1);
		}
	}

	return true;
}

bool CurveContour::extractContourEdgesByFaceOrientation(vec3 viewPoint) {

	int numFaces = (int)mpMesh->indices.size();
	int numVertives = mpMesh->amount;

	// mark edge flags

	map<vec2i, int> edgeMap; // edge indices => flag (1: positive flag; 2: negative flag)

	for (vec3i &faceIdx : mpMesh->indices) {
		vec3 faceP[3];
		for (int k = 0; k < 3; k++) faceP[k] = mpMesh->positions[faceIdx[k]];
		vec3 faceNormal = cml::cross(faceP[1] - faceP[0], faceP[2] - faceP[0]);
		vec3 faceCenter = (faceP[0] + faceP[1] + faceP[2]) / 3;
		float sign = cml::dot(viewPoint - faceCenter, faceNormal);
		int flag = sign>0 ? 1 : sign<0 ? 2 : 0;

		for (int k = 0; k < 3; k++) {
			vec2i edgeKey = makeKey(faceIdx[k], faceIdx[(k + 1) % 3]);
			auto it = edgeMap.find(edgeKey);
			if (it == edgeMap.end()) {
				edgeMap[edgeKey] = flag;
			} else {
				it->second |= flag;
			}
		}
	}

	// extract all contour edges

	mCurvePoints.clear();
	vector<int> vertexMap(numVertives, -1); // curve point ID : # of vertices
	vector<vec2i> contourEdges; // curve point ID pair : # of edges
	for (auto it : edgeMap) {
		if (it.second == 3) {
			vec2i edge = it.first;
			if (vertexMap[edge[0]] < 0) {
				vertexMap[edge[0]] = (int)mCurvePoints.size();
				mCurvePoints.push_back(mpMesh->positions[edge[0]]);
			}
			if (vertexMap[edge[1]] < 0) {
				vertexMap[edge[1]] = (int)mCurvePoints.size();
				mCurvePoints.push_back(mpMesh->positions[edge[1]]);
			}
			contourEdges.push_back(vec2i(vertexMap[edge[0]], vertexMap[edge[1]]));
		}
	}
	mCurvePointGraph.assign(mCurvePoints.size(), set<int>());
	for (vec2i &edge : contourEdges) {
		mCurvePointGraph[edge[0]].insert(edge[1]);
		mCurvePointGraph[edge[1]].insert(edge[0]);
	}
	
	// merge duplicated points

	if (true) {
		vector<int> vertexIndices;
		vector<vec3> mergedCurvePoints;
		vector<set<int>> mergedCurvePointGraph;
		if (!MeshUtil::removeDuplicateVertices(mCurvePoints, mergedCurvePoints, vertexIndices, 1e-5)) return false;
		mergedCurvePointGraph.assign(mergedCurvePoints.size(), set<int>());
		for (int id = 0; id < (int)mCurvePoints.size(); id++) {
			int newID = vertexIndices[id];
			auto &newGraph = mergedCurvePointGraph[newID];
			for (int nbID : mCurvePointGraph[id]) {
				int newNBID = vertexIndices[nbID];
				if (newNBID != newID) newGraph.insert(newNBID);
			}
		}

		mCurvePoints.swap(mergedCurvePoints);
		mCurvePointGraph.swap(mergedCurvePointGraph);
	}

	return true;
}

bool CurveContour::chainContourEdges() {

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

	mContourChains.clear();
	set<vec2i> checkedEdgeSet;

	for (int pointID : checkPointList) {
		for (int neighborID : mCurvePointGraph[pointID]) {

			vector<vec3> chain(1, mCurvePoints[pointID]);

			vec2i currentSegment(pointID, neighborID);

			while (checkedEdgeSet.find(makeKey(currentSegment)) == checkedEdgeSet.end()) {
				checkedEdgeSet.insert(makeKey(currentSegment));

				// add current segment				
				int lastPointID = currentSegment[0];
				int currentPointID = currentSegment[1];
				vec3 lastPoint = mCurvePoints[lastPointID];
				vec3 currentPoint = mCurvePoints[currentPointID];
				chain.push_back(currentPoint);

				// get next segment
				
				auto &nextNbSet = mCurvePointGraph[currentPointID];
				if ((int)nextNbSet.size() != 2) break;
				int nextPointID = -1;
				for (int id : nextNbSet) {
					if (id != lastPointID) nextPointID = id;
				}
				vec3 nextPoint = mCurvePoints[nextPointID];
				vec3d lastDir = currentPoint - lastPoint;
				vec3d nextDir = nextPoint - currentPoint;
				if (cml::unsigned_angle(lastDir, nextDir) > cml::rad(180.0)) break; // UNDONE: param
				currentSegment = vec2i(currentPointID, nextPointID);
			}

			// add chain

			if ((int)chain.size() > 1) {
				mContourChains.push_back(chain);
			}
		}
	}

	return true;
}

bool CurveContour::output(vector<vector<vec3>> &curves) {

	curves = mContourChains;

	return true;
}

bool CurveContour::visualize(string fileName) {

	if (!CurveUtil::visualizeCurves(fileName, mContourChains, 0)) return false;

	return true;
}
