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

#include "CurveBoundary.h"

#include <iostream>
#include <fstream>
#include <map>
#include <set>

#include "Curve/CurveUtil.h"

#include "Utility/PlyExporter.h"

#include "Data/StyleSynthesisConfig.h"

using namespace std;
using namespace StyleSynthesis;

CurveBoundary::CurveBoundary(TTriangleMesh &mesh) {

	mpMesh = &mesh;
}

CurveBoundary::~CurveBoundary() {

}

bool CurveBoundary::extractCurve() {

	if (!extractBoundaryEdges()) return false;
	if (!chainBoundaryEdges()) return false;

	return true;
}

bool CurveBoundary::extractBoundaryEdges() {

	int numFaces = (int)mpMesh->indices.size();

	map<vec2i, vec3> edgeSet; // edge key => any neighbor face normal
	set<vec2i> boundaryEdgeSet;

	for (vec3i faceIdx : mpMesh->indices) {
		vec3 faceP[3];
		for (int k = 0; k < 3; k++) faceP[k] = mpMesh->positions[faceIdx[k]];
		vec3 faceN = cml::normalize(cml::cross(faceP[1] - faceP[0], faceP[2] - faceP[0]));
		for (int k = 0; k < 3; k++) {
			vec2i edgeKey = makeKey(faceIdx[k], faceIdx[(k + 1) % 3]);
			auto it = edgeSet.find(edgeKey);
			if (it == edgeSet.end()) {
				edgeSet[edgeKey] = faceN;
			} else {
				vec3 edgeN = it->second;
				float deg = cml::deg(cml::acos_safe(fabs(cml::dot(faceN, edgeN))));
				if (deg > 20.0f) { // UNDONE: param boundary edge normal angle threshold
					boundaryEdgeSet.insert(edgeKey);
				}
			}
		}
	}

	mBoundaryEdges.assign(boundaryEdgeSet.begin(), boundaryEdgeSet.end());

	return true;
}

bool CurveBoundary::chainBoundaryEdges() {

	// build graph

	int numPoints = mpMesh->amount;
	vector<set<int>> pointGraph(numPoints, set<int>());
	for (vec2i &edgeKey : mBoundaryEdges) {
		pointGraph[edgeKey[0]].insert(edgeKey[1]);
		pointGraph[edgeKey[1]].insert(edgeKey[0]);
	}

	// determine start point checking order (end point first then inner point (for looping curve))

	vector<int> endPointSet;
	vector<int> innerPointSet;
	for (int pointID = 0; pointID < numPoints; pointID++) {
		int numNeighbors = (int)pointGraph[pointID].size();
		if (numNeighbors == 2) innerPointSet.push_back(pointID);
		else if (numNeighbors > 0) endPointSet.push_back(pointID);
	}
	vector<int> checkPointList(0);
	checkPointList.insert(checkPointList.end(), endPointSet.begin(), endPointSet.end());
	checkPointList.insert(checkPointList.end(), innerPointSet.begin(), innerPointSet.end());

	// start chaining

	mBoundaryChains.clear();
	set<vec2i> checkedEdgeSet;

	for (int pointID : checkPointList) {
		for (int neighborID : pointGraph[pointID]) {

			vector<vec3> chain(1, mpMesh->positions[pointID]);

			vec2i currentSegment(pointID, neighborID);

			while (checkedEdgeSet.find(makeKey(currentSegment)) == checkedEdgeSet.end()) {
				checkedEdgeSet.insert(makeKey(currentSegment));

				// add current segment				
				int lastPointID = currentSegment[0];
				int currentPointID = currentSegment[1];
				vec3 lastPoint = mpMesh->positions[lastPointID];
				vec3 currentPoint = mpMesh->positions[currentPointID];
				chain.push_back(currentPoint);

				// get next segment
				
				auto &nextNbSet = pointGraph[currentPointID];
				if ((int)nextNbSet.size() != 2) break;
				int nextPointID = -1;
				for (int id : nextNbSet) {
					if (id != lastPointID) nextPointID = id;
				}
				vec3 nextPoint = mpMesh->positions[nextPointID];
				vec3d lastDir = currentPoint - lastPoint;
				vec3d nextDir = nextPoint - currentPoint;
				if (cml::unsigned_angle(lastDir, nextDir) > cml::rad(30.0)) break;
				currentSegment = vec2i(currentPointID, nextPointID);
			}

			// add chain

			if ((int)chain.size() > 1) {
				mBoundaryChains.push_back(chain);
			}
		}
	}

	if (!CurveUtil::removeDuplicateLines(mBoundaryChains)) return false;

	return true;
}

bool CurveBoundary::output(vector<vector<vec3>> &curves) {

	curves = mBoundaryChains;

	return true;
}

bool CurveBoundary::visualize(string fileName) {

	if (!CurveUtil::visualizeCurves(fileName, mBoundaryChains, 0)) return false;

	return true;
}
