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

#include "CurveUtil.h"

#include <fstream>
#include <set>

#include "Utility/PlyExporter.h"
#include "Mesh/MeshUtil.h"
#include "Segment/SegmentUtil.h"

#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

bool CurveUtil::makeTube(vector<vec3> &rawCurve, TTriangleMesh &tube, float radius) {

	int tubeSubDivision = 10;

	tube.positions.clear();
	tube.normals.clear();
	tube.indices.clear();

	vector<vec3> curve;
	if (!cleanCurve(rawCurve, curve)) return false;

	int numPoints = (int)curve.size();
	for (int pointID = 0; pointID < numPoints; pointID++) {
		vec3 pointPos = curve[pointID];
		vec3 pointDir;		
		if (pointID == 0) {
			pointDir = curve[pointID + 1] - curve[pointID];
		} else if (pointID == numPoints - 1) {
			pointDir = curve[pointID] - curve[pointID - 1];
		} else {			
			pointDir = curve[pointID + 1] - curve[pointID - 1];
			if (pointDir.length_squared() == 0) {
				pointDir = curve[pointID] - curve[pointID - 1];
			}
		}
		pointDir.normalize();

		vec3 side = fabs(pointDir[0]) < 0.5f ? vec3(1.0f, 0.0f, 0.0f) : vec3(0.0f, 1.0f, 0.0f);
		vec3 norm1 = cml::normalize(cml::cross(side, pointDir));
		vec3 norm2 = cml::normalize(cml::cross(pointDir, norm1));
		for (int divID = 0; divID < tubeSubDivision; divID++) {
			float angle = divID * cml::constantsf::two_pi() / tubeSubDivision;
			vec3 vn = (cos(angle) * norm1 + sin(angle) * norm2).normalize();
			vec3 vp = vn * radius + pointPos;
			tube.positions.push_back(vp);
			tube.normals.push_back(vn);
		}
	}
	tube.amount = (int)tube.positions.size();
	numPoints = tube.amount / tubeSubDivision;

	for (int pointID = 0; pointID < numPoints - 1; pointID++) {

		int offset = 0;
		if (true) {
			float minDist = FLT_MAX;
			int nowStart = pointID * tubeSubDivision;
			int nextStart = (pointID + 1)*tubeSubDivision;
			for (int nowID = 0; nowID < tubeSubDivision; nowID++) {
				for (int nextID = 0; nextID < tubeSubDivision; nextID++) {
					float dist = (tube.positions[nowStart + nowID] - tube.positions[nextStart + nextID]).length_squared();
					if (dist < minDist) {
						minDist = dist;
						offset = (nextID - nowID + tubeSubDivision) % tubeSubDivision;
					}
				}
			}
		}
		for (int divID = 0; divID < tubeSubDivision; divID++) {
			int idx00 = pointID*tubeSubDivision + divID;
			int idx01 = pointID*tubeSubDivision + (divID + 1) % tubeSubDivision;
			int idx10 = (pointID + 1)*tubeSubDivision + (divID + offset) % tubeSubDivision;
			int idx11 = (pointID + 1)*tubeSubDivision + (divID + offset + 1) % tubeSubDivision;
			tube.indices.push_back(vec3i(idx00, idx01, idx10));
			tube.indices.push_back(vec3i(idx10, idx01, idx11));
		}
	}

	return true;
}

bool CurveUtil::cleanCurve(vector<vec3> &inCurve, vector<vec3> &outCurve) {

	vector<vec3> newCurve(0);
	for (vec3 p : inCurve) {
		if (newCurve.empty() || p != newCurve.back()) newCurve.push_back(p);
	}
	outCurve.swap(newCurve);

	return true;
}


bool CurveUtil::buildMatrix(vector<vec3> &inCurve, Eigen::Matrix3Xd &outMatrix) {

	outMatrix.resize(3, (int)inCurve.size());

#pragma omp parallel for
	for (int j = 0; j < (int)inCurve.size(); j++) {
		outMatrix.col(j) = Eigen::Vector3d(vec3d(inCurve[j]).data());
	}

	return true;
}

bool CurveUtil::computeCurveDirections(vector<vec3> &inCurve, vector<vec3> &outDirections) {

	int numPoints = (int)inCurve.size();
	outDirections.resize(numPoints);

#pragma omp parallel for
	for (int pointID = 0; pointID < numPoints; pointID++) {
		if (pointID == 0) {
			outDirections[pointID] = inCurve[pointID + 1] - inCurve[pointID];
		}
		else if (pointID == numPoints - 1) {
			outDirections[pointID] = inCurve[pointID] - inCurve[pointID - 1];
		}
		else {
			outDirections[pointID] = inCurve[pointID + 1] - inCurve[pointID - 1];
		}
		float len = outDirections[pointID].length();
		if (len) outDirections[pointID] /= len;
	}

	return true;
}

bool CurveUtil::computeCurveLength(vector<vec3> &inCurve, double &outLength) {

	int numSegments = (int)inCurve.size();
	if (numSegments < 2) {
		outLength = 0;
		return true;
	}

	outLength = 0;
	for (int segID = 0; segID < numSegments - 1; segID++) {
		outLength += (vec3d(inCurve[segID]) - vec3d(inCurve[segID + 1])).length();
	}

	return true;
}

bool CurveUtil::checkStraightCurve(vector<vec3> &inCurve, bool &outFlag) {

	double curvEps = 0.01;
	double angleEps = 1e-6;
	double prctEps = 0.9;

	double currentSegLen = (inCurve[0] - inCurve[1]).length();
	double totalLen = currentSegLen;
	double totalCurvature = 0;
	double totalWeights = 0;
	double maxSegLen = currentSegLen;
	double maxCurvature = 0;
	for (int segID = 0; segID < (int)inCurve.size() - 2; segID++) {
		double curv = 1 - cml::dot(cml::normalize(inCurve[segID] - inCurve[segID + 1]),
			cml::normalize(inCurve[segID + 1] - inCurve[segID + 2])); // not really "curvature" but somehow like that...

		double segLen0 = (inCurve[segID] - inCurve[segID + 1]).length();
		double segLen1 = (inCurve[segID + 1] - inCurve[segID + 2]).length();

		if (curv > angleEps) currentSegLen = 0;
		currentSegLen += segLen1;
		totalLen += segLen1;

		double weight = segLen0 + segLen1;
		totalCurvature += curv * weight;
		totalWeights += weight;

		maxSegLen = max(maxSegLen, currentSegLen);
		maxCurvature = max(maxCurvature, curv);
	}

	outFlag = true;
	if (maxSegLen < totalLen * prctEps) {
		if (totalWeights) totalCurvature /= totalWeights;
		if (totalCurvature > curvEps) {
			outFlag = false;
		}
	}

	return true;
}

bool CurveUtil::chainLines(vector<vec2i> &edges, vector<vec3> &vertices, vector<vector<int>> &indices, double maxAngle) {

	// chain edge lines with compatible directions into curves

	int numPoints = (int)vertices.size();
	int numEdges = (int)edges.size();

	// build graph

	vector<set<int>> neighborSet(numPoints, set<int>());
	for (vec2i &idx : edges) {
		neighborSet[idx[0]].insert(idx[1]);
		neighborSet[idx[1]].insert(idx[0]);
	}

	// build edge ID map

	map<vec2i, int> edgeMap;
	for (int edgeID = 0; edgeID < numEdges; edgeID++) {
		vec2i key = edges[edgeID];
		if (key[0] > key[1]) swap(key[0], key[1]);
		edgeMap[key] = edgeID;
	}

	// chain edges

	vector<bool> edgeFlags(numEdges, false); // visited flag
	for (int edgeID = 0; edgeID < numEdges; edgeID++) {
		if (edgeFlags[edgeID]) continue;
		edgeFlags[edgeID] = true;
		vec2i edge = edges[edgeID];

		// chain from two directions

		vector<int> currentLine(0);
		set<int> currentVertexSet;
		currentVertexSet.insert(edge[0]);
		currentVertexSet.insert(edge[1]);
		bool innerLoop = false;
		for (int alterID = 0; alterID < 2; alterID++) {
			int lastPointID = edge[1 - alterID];
			int currentPointID = edge[alterID];
			vec3 currentDirection = cml::normalize(vertices[currentPointID] - vertices[lastPointID]);
			vector<int> chain(0);
			chain.push_back(lastPointID);
			chain.push_back(currentPointID);
			while (true) {
				if (currentPointID == chain[0]) {
					break; // loop detected
				}
				float maxCosine = (float)cos(cml::rad(maxAngle));
				vec3 maxDirection;
				int maxID = -1;
				for (int neighborPointID : neighborSet[currentPointID]) {
					if (neighborPointID == lastPointID) continue;
					vec3 neighborDirection = cml::normalize(vertices[neighborPointID] - vertices[currentPointID]);
					float neighborCosine = cml::dot(currentDirection, neighborDirection);
					if (neighborCosine > maxCosine) {
						maxCosine = neighborCosine;
						maxDirection = neighborDirection;
						maxID = neighborPointID;
					}
				}
				if (maxID >= 0) {
					vec2i key(currentPointID, maxID);
					if (key[0] > key[1]) swap(key[0], key[1]);
					int nextEdgeID = edgeMap[key];
					if (edgeFlags[nextEdgeID]) {
						break; // loop detected
					}
					edgeFlags[nextEdgeID] = true;
					lastPointID = currentPointID;
					currentPointID = maxID;
					currentDirection = maxDirection;
					chain.push_back(maxID);
					if (currentVertexSet.find(maxID) != currentVertexSet.end() && maxID != chain[0]) {
						innerLoop = true;
						break;
					}
					currentVertexSet.insert(maxID);
				} else {
					break; // no next edge
				}
			}
			if (chain.front() == chain.back()) {
				// real loop -- use this
				currentLine = chain;
				break;
			}
			if (alterID == 0) {
				currentLine.assign(chain.rbegin(), chain.rend());
			} else {
				if ((int)chain.size() > 2) {
					currentLine.insert(currentLine.end(), chain.begin() + 2, chain.end());
				}
			}
			if (innerLoop) break;
		}
		if (innerLoop) {
			// skip chaining from this edge
			int numEdgeSeg = (int)currentLine.size() - 1;
			for (int id = 0; id < numEdgeSeg; id++) {
				vec2i key(currentLine[id], currentLine[id + 1]);
				if (key[0] > key[1]) swap(key[0], key[1]);
				int edgeID = edgeMap[key];
				edgeFlags[edgeID] = false; // reset edge flags
			}
			continue;
		}

		if (true) {
			// sanity check
			currentVertexSet.clear();
			for (int id = 0; id < (int)currentLine.size() - 1; id++) {
				int vertID = currentLine[id];
				if (currentVertexSet.find(vertID) != currentVertexSet.end()) {
					cout << "Error: duplicated vertex found" << endl;
					return false;
				}
				currentVertexSet.insert(vertID);
			}
		}

		indices.push_back(currentLine);
	}

	if (!chainCurves(vertices, indices, maxAngle)) return false;

	return true;
}

bool CurveUtil::chainCurves(vector<vec3> &vertices, vector<vector<int>> &indices, double maxAngle) {

	// chain curves with compatible directions into longer curves

	int numCurves = (int)indices.size();
	int numVertices = (int)vertices.size();

	vector<vector<int>> newCurves(0);

	// build graph

	vector<set<int>> graph(numVertices, set<int>());
	vector<bool> curveFlags(numCurves, false); // visited flag
	for (int curveID = 0; curveID < numCurves; curveID++) {
		auto &curve = indices[curveID];
		int pointID1 = curve.front();
		int pointID2 = curve.back();
		if (pointID1 != pointID2) {
			graph[pointID1].insert(curveID);
			graph[pointID2].insert(curveID);
		} else {
			curveFlags[curveID] = true; // skip chaining looping curve
			// shift loop start point to most "curvy" point
			int numPoints = (int)curve.size();
			double maxAngle = 0;
			double maxAngleSum = -DBL_MAX;
			int maxID = 0;
			for (int pointID = 0; pointID < numPoints; pointID++) {
				vec3d midPoint = vertices[curve[pointID]];
				vec3d prevPoint = vertices[curve[(pointID + numPoints - 1) % numPoints]];
				vec3d nextPoint = vertices[curve[(pointID + 1) % numPoints]];
				double angle = cml::unsigned_angle(midPoint - prevPoint, nextPoint - midPoint);
				double angleSum = cml::dot(midPoint, vec3d(1.0, 1.0, 1.0));
				if (angle > maxAngle || (fabs(angle-maxAngle) < 1e-6 && angleSum > maxAngleSum)) {
					maxAngle = angle;
					maxAngleSum = angleSum;
					maxID = pointID;
				}
			}
			vector<int> shiftedCurve(curve.begin() + maxID, curve.end());
			shiftedCurve.insert(shiftedCurve.end(), curve.begin() + 1, curve.begin() + maxID + 1);
			newCurves.push_back(shiftedCurve);
		}
	}

	for (int curveID = 0; curveID < numCurves; curveID++) {
		if (curveFlags[curveID]) continue;
		curveFlags[curveID] = true;

		vec2i endPoints(indices[curveID].front(), indices[curveID].back());

		// chain from two directions

		vector<int> currentChain(0);
		for (int alterID = 0; alterID < 2; alterID++) {

			vector<int> chain(0);
			if (alterID == 0) chain = indices[curveID];

			int lastPointID = endPoints[alterID];
			int currentPointID = endPoints[1 - alterID];
			int currentCurveID = curveID;
			vec3 currentDirection = cml::normalize(vertices[currentPointID] - vertices[lastPointID]);

			while (true) {
				float maxCosine = (float)cos(cml::rad(maxAngle));
				vec3 maxDirection;
				int maxNextCurveID = -1;
				int maxNextPointID = -1;
				bool maxNextCurveFlip = false;
				for (int nextCurveID : graph[currentPointID]) {
					if (nextCurveID == currentCurveID) continue;
					int nextPointID = indices[nextCurveID].back();
					bool nextCurveFlip = false;
					if (nextPointID == currentPointID) {
						nextPointID = indices[nextCurveID].front();
						nextCurveFlip = true;
					}
					vec3 neighborDirection = cml::normalize(vertices[nextPointID] - vertices[currentPointID]);
					float neighborCosine = cml::dot(currentDirection, neighborDirection);
					if (neighborCosine > maxCosine) {
						maxCosine = neighborCosine;
						maxDirection = neighborDirection;
						maxNextPointID = nextPointID;
						maxNextCurveID = nextCurveID;
						maxNextCurveFlip = nextCurveFlip;
					}
				}
				if (maxNextCurveID >= 0) {
					if (curveFlags[maxNextCurveID]) {
						break; // loop detected
					}
					curveFlags[maxNextCurveID] = true;
					lastPointID = currentPointID;
					currentPointID = maxNextPointID;
					currentDirection = maxDirection;
					currentCurveID = maxNextCurveID;
					auto &maxNextCurve = indices[maxNextCurveID];
					if (maxNextCurveFlip) {
						chain.insert(chain.end(), maxNextCurve.rbegin() + 1, maxNextCurve.rend());
					} else {
						chain.insert(chain.end(), maxNextCurve.begin() + 1, maxNextCurve.end());
					}
				} else {
					break; // no next curve
				}
			}
			if (alterID == 0) {
				currentChain.assign(chain.rbegin(), chain.rend());
			} else {
				if ((int)chain.size() > 0) {
					currentChain.insert(currentChain.end(), chain.begin(), chain.end());
				}
			}
		}

		newCurves.push_back(currentChain);
	}

	indices.swap(newCurves);

	return true;
}

bool CurveUtil::extractLines(vector<vector<int>> &indices, vector<vec3> &vertices, vector<vector<vec3>> &lines) {

	int numLines = (int)indices.size();
	lines.resize(numLines);

	for (int lineID = 0; lineID < numLines; lineID++) {
		auto &inLine = indices[lineID];
		auto &outLine = lines[lineID];
		int numPoints = (int)inLine.size();
		outLine.resize(numPoints);
		for (int pointID = 0; pointID < numPoints; pointID++) {
			outLine[pointID] = vertices[inLine[pointID]];
		}
	}

	return true;
}

bool CurveUtil::removeDuplicatePoints(vector<vec2i> &edges, vector<vec3> &vertices) {

	vector<int> vertexIndices;
	if (!MeshUtil::removeDuplicateVertices(vertices, vertices, vertexIndices)) return false;

	for (vec2i &edge : edges) {
		edge[0] = vertexIndices[edge[0]];
		edge[1] = vertexIndices[edge[1]];
	}

	return true;
}

bool CurveUtil::removeDuplicateLines(vector<vector<vec3>> &curves) {

	// only used for ground truth lines from mesh edges

	vec3 bbMin(FLT_MAX, FLT_MAX, FLT_MAX);
	vec3 bbMax(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	for (auto &curve : curves) {
		for (vec3 &v : curve) {
			bbMax.maximize(v);
			bbMin.minimize(v);
		}
	}
	float bbLen = (bbMax - bbMin).length();
	double maxDist = bbLen * 0.001; // UNDONE: param maximum distance for duplicated curves

	int numCurves = (int)curves.size();

	vector<bool> curveFlags(numCurves, true);

	for (int id1 = 0; id1 < numCurves - 1; id1++) {
		if (!curveFlags[id1]) continue;
		auto &curve1 = curves[id1];
		for (int id2 = id1 + 1; id2 < numCurves; id2++) {
			if (!curveFlags[id1]) continue;
			auto &curve2 = curves[id2];
			if (curve1.size() != curve2.size()) continue;
			int n = (int)curve1.size();
			double sumDist1 = 0;
			double sumDist2 = 0;
			for (int k = 0; k < n; k++) {
				sumDist1 += (curve1[k] - curve2[k]).length();
				sumDist2 += (curve1[n - k - 1] - curve2[k]).length();
			}
			if (sumDist1 < maxDist*n || sumDist2 < maxDist*n) {
				curveFlags[id2] = false;
			}
		}
	}

	vector<vector<vec3>> newCurves(0);
	int removedCount = 0;
	for (int curveID = 0; curveID < numCurves; curveID++) {
		if (curveFlags[curveID]) newCurves.push_back(curves[curveID]);
		else removedCount++;
	}
	//cout << "Removed " << removedCount << " out of " << numCurves << " curves" << endl;

	curves.swap(newCurves);

	return true;
}

bool CurveUtil::sampleLine(vector<vec3> &curve, vector<vec3> &sample, float radius) {

	sample.clear();
	int numSegments = (int)curve.size() - 1;
	for (int segID = 0; segID < numSegments; segID++) {
		vec3 p1 = curve[segID];
		vec3 p2 = curve[segID + 1];
		float len = (p2 - p1).length();
		int n = (int)ceil(len / radius);
		if (n <= 0) continue;
		vec3 delta = cml::normalize(p2 - p1) * (len / n);
		if (segID == 0) sample.push_back(p1);
		for (int k = 1; k <= n; k++) {
			sample.push_back(p1 + delta * k);
		}
	}

	return true;
}

bool CurveUtil::sampleLines(vector<vector<vec3>> &curves, vector<vector<vec3>> &samples, float radius) {

	int numCurves = (int)curves.size();
	samples.resize(numCurves);
	for (int curveID = 0; curveID < numCurves; curveID++) {
		if (!sampleLine(curves[curveID], samples[curveID], radius)) return false;
	}

	return true;
}

bool CurveUtil::filterShortLines(vector<vector<vec3>> &curves, float length) {

	vector<vector<vec3>> filteredCurves(0);
	for (auto &curve : curves) {
		double curveLen;
		if (!computeCurveLength(curve, curveLen)) return false;
		if (curveLen > (double)length) filteredCurves.push_back(curve);
	}

	curves.swap(filteredCurves);

	return true;
}

bool CurveUtil::saveCurves(string fileName, vector<vector<vec3>> &curves) {

	ofstream file(fileName, ios::binary);
	if (!file.is_open()) {
		cout << "Error: cannot write to file " << fileName << endl;
		return false;
	}
	if (!saveCurves(file, curves)) return false;
	file.close();

	return true;
}

bool CurveUtil::saveCurves(ostream &fileStream, vector<vector<vec3>> &curves) {

	int numCurves = (int)curves.size();
	fileStream.write((char*)(&numCurves), sizeof(numCurves));

	for (auto &curve : curves) {

		int numPoints = (int)curve.size();
		fileStream.write((char*)(&numPoints), sizeof(numPoints));

		for (vec3 &p : curve) {
			fileStream.write((char*)p.data(), sizeof(p));
		}
	}

	return true;
}

bool CurveUtil::loadCurves(string fileName, vector<vector<vec3>> &curves) {

	ifstream file(fileName, ios::binary);
	if (!file.is_open()) {
		cout << "Error: cannot open file " << fileName << endl;
		return false;
	}
	if (!loadCurves(file, curves)) return false;
	file.close();

	return true;
}

bool CurveUtil::loadCurves(istream &fileStream, vector<vector<vec3>> &curves) {

	int numCurves;
	fileStream.read((char*)(&numCurves), sizeof(numCurves));
	curves.resize(numCurves);

	for (auto &curve : curves) {

		int numPoints;
		fileStream.read((char*)(&numPoints), sizeof(numPoints));
		curve.resize(numPoints);

		for (int pointID = 0; pointID < numPoints; pointID++) {
			vec3 &point = curve[pointID];
			fileStream.read((char*)point.data(), sizeof(point));
		}
	}

	return true;
}

bool CurveUtil::saveCurveGroups(string fileName, vector<vector<vector<vec3>>> &curves) {

	ofstream file(fileName, ios::binary);
	if (!file.is_open()) {
		cout << "Error: cannot write to file " << fileName << endl;
		return false;
	}

	int numGroups = (int)curves.size();
	file.write((char*)&numGroups, sizeof(numGroups));
	for (auto &group : curves) {
		if (!saveCurves(file, group)) return false;
	}

	file.close();

	return true;
}

bool CurveUtil::loadCurveGroups(string fileName, vector<vector<vector<vec3>>> &curves) {

	ifstream file(fileName, ios::binary);
	if (!file.is_open()) {
		cout << "Error: cannot open file " << fileName << endl;
		return false;
	}

	int numGroups;
	file.read((char*)&numGroups, sizeof(numGroups));
	curves.resize(numGroups);
	for (auto &group : curves) {
		if (!loadCurves(file, group)) return false;
	}

	file.close();

	return true;
}

bool CurveUtil::visualizeCurves(string fileName, vector<vector<vec3>> &curves, int style) {

	// style 0: colored tube mesh
	// style 1: lines without colors
	// style 2: samples on lines without colors
	// style 3: colored samples on lines

	if (style == 0) {

		// colored tube mesh

		int tubeSubDivision = 10;

		vec3 bbMin(FLT_MAX, FLT_MAX, FLT_MAX);
		vec3 bbMax(-FLT_MAX, -FLT_MAX, -FLT_MAX);
		for (auto &curve : curves) {
			for (vec3 &v : curve) {
				bbMin.minimize(v);
				bbMax.maximize(v);
			}
		}
		float radius = (bbMax - bbMin).length() * 0.002f;

		PlyExporter pe;
		for (int curveID = 0; curveID < (int)curves.size(); curveID++) {

			TTriangleMesh curveMesh;
			if (!CurveUtil::makeTube(curves[curveID], curveMesh, radius)) return false;

			vec3i color = SegmentUtil::colorMapping(curveID);
			vector<vec3i> tubeC(curveMesh.indices.size(), color);
			if (!pe.addMesh(&curveMesh.positions, &curveMesh.normals, &curveMesh.indices, &tubeC)) return false;
		}
		if (!pe.output(fileName)) return false;

	} else if (style == 1) {

		// lines without colors
				
		vector<vec3> curvePoints;
		for (int curveID = 0; curveID < (int)curves.size(); curveID++) {
			auto &curve = curves[curveID];
			int numPoints = (int)curve.size();
			for (int segID = 0; segID < numPoints - 1; segID++) {
				curvePoints.push_back(curve[segID]);
				curvePoints.push_back(curve[segID + 1]);
			}
		}

		PlyExporter pe;
		if (!pe.addLine(&curvePoints)) return false;
		if (!pe.output(fileName)) return false;

	} else if (style == 2) {

		// samples on lines without colors

		vec3 bbMin(FLT_MAX, FLT_MAX, FLT_MAX);
		vec3 bbMax(-FLT_MAX, -FLT_MAX, -FLT_MAX);
		for (auto &curve : curves) {
			for (vec3 &v : curve) {
				bbMin.minimize(v);
				bbMax.maximize(v);
			}
		}
		float radius = (bbMax - bbMin).length() * 0.001f;

		vector<vector<vec3>> samples(0);
		if (!sampleLines(curves, samples, radius)) return false;

		PlyExporter pe;
		for (auto &sample : samples) {
			if (!pe.addPoint(&sample)) return false;
		}
		if (!pe.output(fileName)) return false;

	} else if (style == 3) {

		// colored samples on lines

		vec3 bbMin(FLT_MAX, FLT_MAX, FLT_MAX);
		vec3 bbMax(-FLT_MAX, -FLT_MAX, -FLT_MAX);
		for (auto &curve : curves) {
			for (vec3 &v : curve) {
				bbMin.minimize(v);
				bbMax.maximize(v);
			}
		}
		float radius = (bbMax - bbMin).length() * 0.001f;

		vector<vector<vec3>> curveSamples(0);
		if (!sampleLines(curves, curveSamples, radius)) return false;

		PlyExporter pe;
		for (int curveID = 0; curveID < (int)curves.size(); curveID++) {
			vec3i color = SegmentUtil::colorMapping(curveID);
			if (!pe.addPoint(&curveSamples[curveID], 0, cml::identity_4x4(), color)) return false;
		}
		if (!pe.output(fileName)) return false;

	} else {
		cout << "Error: unknown style " << style << endl;
		return false;
	}

	return true;
}

bool CurveUtil::generateViewPoints(vector<vec3> &viewPoints) {

	viewPoints.clear();

	float viewPointRadius = 2.0f;
	for (int vRot = 0; vRot < 3; vRot++) {
		float vAngle = cml::rad(30.0f) * (1 - vRot);
		for (int hRot = 0; hRot < 12; hRot++) {
			float hAngle = cml::rad(30.0f) * hRot;
			vec3 viewPoint = vec3(sin(hAngle) * cos(vAngle), sin(vAngle), cos(hAngle) * cos(vAngle)) * viewPointRadius;
			viewPoints.push_back(viewPoint);
		}
	}

	return true;
}

bool CurveUtil::saveViewPoints(string fileName, vector<vec3> &viewPoints) {

	ofstream file(fileName);
	if (!file.is_open()) {
		cout << "Error: cannot write to file " << fileName << endl;
		return false;
	}

	for (vec3 viewPoint : viewPoints) {
		file << viewPoint << endl;
	}

	file.close();

	return true;
}

bool CurveUtil::loadViewPoints(string fileName, vector<vec3> &viewPoints) {

	viewPoints.clear();

	ifstream file(fileName);
	if (!file.is_open()) {
		cout << "Error: cannot open file " << fileName << endl;
		return false;
	}

	while (!file.eof()) {
		vec3 viewPoint;
		file >> viewPoint[0] >> viewPoint[1] >> viewPoint[2];
		if (file.fail()) break;
		if(viewPoint.length_squared()) viewPoint.normalize();
		viewPoints.push_back(viewPoint);
	}

	file.close();

	return true;
}