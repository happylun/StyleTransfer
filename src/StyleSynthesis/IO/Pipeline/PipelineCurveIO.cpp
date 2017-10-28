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

#include "PipelineCurveIO.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>

#include "Mesh/MeshUtil.h"
#include "Sample/SampleUtil.h"
#include "Segment/SegmentUtil.h"
#include "Curve/CurveUtil.h"

#include "Sample/SamplePoissonDisk.h"

#include "Curve/CurveSupportLines.h"
#include "Curve/CurveRidgeValley.h"
#include "Curve/CurveBoundary.h"
#include "Curve/CurveContour.h"

#include "Curve/CurveLineSnapping.h"
#include "Curve/CurveLineVoting.h"

#include "Utility/PlyExporter.h"

#include "Data/StyleSynthesisConfig.h"

using namespace std;
using namespace StyleSynthesis;

//#define OUTPUT_PROGRESS

bool PipelineCurveIO::process() {

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string meshListFileName = datasetPrefix + "mesh/mesh-list-all.txt";
	vector<string> meshNameList;
	ifstream meshListFile(meshListFileName);
	while (!meshListFile.eof()) {
		string line;
		getline(meshListFile, line);
		if (line.empty()) break;
		meshNameList.push_back(StringUtil::trim(line));
	}
	meshListFile.close();

	int numMesh = (int)meshNameList.size();

	for (int meshID = 0; meshID < numMesh; meshID++) {

		string meshName = meshNameList[meshID];
		cout << "Processing curve " << meshName << endl;
		if (!runSingleModelSupportLines(meshName)) return error("single model SL error");
		if (!runSingleModelRidgesValleys(meshName)) return error("single model RV error");
		if (!runSingleModelBoundaries(meshName)) return error("single model B error");
		if (!runSingleModelContours(meshName)) return error("single model C error");
		if (!runSingleModelSnapping(meshName)) return error("single model snap error");

		//system("pause");
	}

	return true;
}

bool PipelineCurveIO::runSingleModelSupportLines(string meshName) {

#ifdef OUTPUT_PROGRESS
	cout << "Running support lines..." << endl;
#endif

	// names...

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string scaledPartFolder = datasetPrefix + "scaled-part/" + meshName + "/";

	string curvePrefix = datasetPrefix + "curve/" + meshName + "/";
	if (!FileUtil::makedir(curvePrefix)) return false;

	string slName = curvePrefix + "data-sl.txt";
	string slVName = curvePrefix + "vis-sl.ply";

	if (FileUtil::existsfile(slName)) return true; // early quit

	vector<vector<vec3>> allCurves(0);

	// process each part separately

	int numParts = 0;
	while (true) {
		string partName = scaledPartFolder + to_string(numParts) + ".ply";
		if (!FileUtil::existsfile(partName)) break;

		// data

		TTriangleMesh scaledMesh;
		TTriangleMesh weldMesh;
		vector<vector<vec3>> curves;

		// mesh

		if (!MeshUtil::loadMesh(partName, scaledMesh)) return false;
		if (!MeshUtil::recomputeNormals(scaledMesh)) return false;

		// weld mesh

		vector<int> weldFaceIndices;
		if (!MeshUtil::weldMeshFaces(scaledMesh, weldMesh, weldFaceIndices)) return false;

		// algorithm

		if (true) {
			CurveSupportLines csl(weldMesh);
			if (!csl.extractLines()) return false;
			if (!csl.output(curves)) return false;
		}

		allCurves.insert(allCurves.end(), curves.begin(), curves.end());
		
		numParts++;
	}

	if (!CurveUtil::visualizeCurves(slVName, allCurves, 0)) return false;
	if (!CurveUtil::saveCurves(slName, allCurves)) return false;

	return true;
}

bool PipelineCurveIO::runSingleModelRidgesValleys(string meshName) {

#ifdef OUTPUT_PROGRESS
	cout << "Running ridges & valleys..." << endl;
#endif

	// names...

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string scaledMeshName = datasetPrefix + "scaled-mesh/" + meshName + ".ply";

	string curvePrefix = datasetPrefix + "curve/" + meshName + "/";
	if (!FileUtil::makedir(curvePrefix)) return false;

	string rvName = curvePrefix + "data-rv.txt";
	string rvVName = curvePrefix + "vis-rv.ply";

	if (FileUtil::existsfile(rvName)) return true; // early quit

	// data...

	TTriangleMesh scaledMesh;
	vector<vector<vec3>> curves;

	// mesh

	if (!MeshUtil::loadMesh(scaledMeshName, scaledMesh)) return false;
	vector<int> degenFaceIndices;
	if (!MeshUtil::removeDegeneratedFaces(scaledMesh, scaledMesh, degenFaceIndices)) return false;
	if (!MeshUtil::recomputeNormals(scaledMesh)) return false;

	// algorithm

	if (true) {
		CurveRidgeValley crv(scaledMesh);
		if (!crv.extractCurve()) return false;
		if (!crv.output(curves)) return false;
		if (!crv.visualize(rvVName)) return false;
		if (!CurveUtil::saveCurves(rvName, curves)) return false;
	}

	return true;
}

bool PipelineCurveIO::runSingleModelBoundaries(string meshName) {

#ifdef OUTPUT_PROGRESS
	cout << "Running boundaries..." << endl;
#endif

	// names...

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string scaledMeshName = datasetPrefix + "scaled-mesh/" + meshName + ".ply";

	string curvePrefix = datasetPrefix + "curve/" + meshName + "/";
	if (!FileUtil::makedir(curvePrefix)) return false;

	string bName = curvePrefix + "data-b.txt";
	string bVName = curvePrefix + "vis-b.ply";

	if (FileUtil::existsfile(bName)) return true; // early quit

	// data...

	TTriangleMesh scaledMesh;
	TTriangleMesh weldMesh;
	vector<vector<vec3>> curves;

	// mesh

	if (!MeshUtil::loadMesh(scaledMeshName, scaledMesh)) return false;
	vector<int> degenFaceIndices;
	if (!MeshUtil::removeDegeneratedFaces(scaledMesh, scaledMesh, degenFaceIndices)) return false;

	vector<int> weldFaceIndices;
	if (!MeshUtil::weldMeshFaces(scaledMesh, weldMesh, weldFaceIndices)) return false;

	// algorithm
	
	if (true) {
		CurveBoundary cb(weldMesh);
		if (!cb.extractCurve()) return false;
		if (!cb.output(curves)) return false;
		if (!cb.visualize(bVName)) return false;
		if (!CurveUtil::saveCurves(bName, curves)) return false;
	}
	
	return true;
}

bool PipelineCurveIO::runSingleModelContours(string meshName) {

#ifdef OUTPUT_PROGRESS
	cout << "Running contours..." << endl;
#endif

	// names...

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string scaledMeshName = datasetPrefix + "scaled-mesh/" + meshName + ".ply";

	string curvePrefix = datasetPrefix + "curve/" + meshName + "/";
	if (!FileUtil::makedir(curvePrefix)) return false;

	string cName = curvePrefix + "data-c.txt";
	string cVName = curvePrefix + "vis-c.ply";

	string viewName = datasetPrefix + "curve/viewpoint.txt";

	if (FileUtil::existsfile(cName)) return true; // early quit

	// data...

	TTriangleMesh scaledMesh;
	vector<vector<vector<vec3>>> curves;

	// mesh

	if (!MeshUtil::loadMesh(scaledMeshName, scaledMesh)) return false;
	vector<int> degenFaceIndices;
	if (!MeshUtil::removeDegeneratedFaces(scaledMesh, scaledMesh, degenFaceIndices)) return false;

	// load view point

	vector<vec3> viewPointList(0);
	if (FileUtil::existsfile(viewName)) {
		if (!CurveUtil::loadViewPoints(viewName, viewPointList)) return false;
	} else {
		if (!CurveUtil::generateViewPoints(viewPointList)) return false;
		if (!CurveUtil::saveViewPoints(viewName, viewPointList)) return false;
	}
	int numViewPoints = (int)viewPointList.size();

	// algorithm

	if (true) {

		curves.assign(numViewPoints, vector<vector<vec3>>(0));

		ofstream curveFile(cName, ios::binary);
		curveFile.write((char*)&numViewPoints, sizeof(numViewPoints));
		for (int viewID = 0; viewID < numViewPoints; viewID++) {
#ifdef OUTPUT_PROGRESS
			cout << "\rRunning view " << (viewID+1) << " / " << numViewPoints << "           ";
#endif
			vec3 viewPoint = viewPointList[viewID];

			vector<vector<vec3>> &viewCurves = curves[viewID];
			CurveContour cc(scaledMesh);
			if (!cc.extractCurve(viewPoint)) return false;
			if (!cc.output(viewCurves)) return false;
			if (!CurveUtil::saveCurves(curveFile, viewCurves)) return false;
		}
		curveFile.close();
#ifdef OUTPUT_PROGRESS
		cout << endl;
#endif
	}

	return true;
}

bool PipelineCurveIO::runSingleModelSnapping(string meshName) {

#ifdef OUTPUT_PROGRESS
	cout << "Running snapping..." << endl;
#endif

	// names...

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string scaledMeshName = datasetPrefix + "scaled-mesh/" + meshName + ".ply";

	string curvePrefix = datasetPrefix + "curve/" + meshName + "/";
	if (!FileUtil::makedir(curvePrefix)) return false;

	string slName = curvePrefix + "data-sl.txt";
	string rvName = curvePrefix + "data-rv.txt";
	string bName = curvePrefix + "data-b.txt";
	string cName = curvePrefix + "data-c.txt";
	string allVPrefix = curvePrefix + "vis-all";

	string snapRVName = curvePrefix + "data-snap-rv.txt";
	string snapBName = curvePrefix + "data-snap-b.txt";
	string snapCName = curvePrefix + "data-snap-c.txt";
	string snapVisPrefix = curvePrefix + "vis-snap";

	string viewName = datasetPrefix + "curve/viewpoint.txt";

	if (FileUtil::existsfile(snapCName)) return true; // early quit

	// data...

	TTriangleMesh scaledMesh;

	// mesh

	if (!MeshUtil::loadMesh(scaledMeshName, scaledMesh)) return false;
	vector<int> degenFaceIndices;
	if (!MeshUtil::removeDegeneratedFaces(scaledMesh, scaledMesh, degenFaceIndices)) return false;

	// thresholds

	vec3 bbMax, bbMin;
	if (!MeshUtil::computeAABB(scaledMesh, bbMin, bbMax)) return false;
	float bbLen = (bbMax - bbMin).length();
	float visualSampleRadius = bbLen * 0.001f;
	float filterCurveLength = bbLen * 0.1f;

	// load view point

	vector<vec3> viewPointList(0);
	if (FileUtil::existsfile(viewName)) {
		if (!CurveUtil::loadViewPoints(viewName, viewPointList)) return false;
	} else {
		if (!CurveUtil::generateViewPoints(viewPointList)) return false;
		if (!CurveUtil::saveViewPoints(viewName, viewPointList)) return false;
	}
	int numViewPoints = (int)viewPointList.size();

	vector<int> visualViewList;
	visualViewList.push_back(10);
	visualViewList.push_back(11);
	visualViewList.push_back(0);
	visualViewList.push_back(1);
	visualViewList.push_back(2);

	if (false) {
		// HACK: skip non-visual views
		vector<vec3> newViewPointList;
		for (int visID : visualViewList) {
			newViewPointList.push_back(viewPointList[visID]);
		}
		numViewPoints = (int)visualViewList.size();
		for (int k = 0; k < numViewPoints; k++) visualViewList[k] = k;
		viewPointList.swap(newViewPointList);
	}

	// load support lines

	vector<vector<vec3>> curveSupportLines;
	if (!FileUtil::existsfile(slName)) {
		return error("run single model SL first");
	} else {
		if (!CurveUtil::loadCurves(slName, curveSupportLines)) return false;
	}
	//if (!CurveUtil::filterShortLines(curveSupportLines, filterCurveLength)) return false;

	// load ridges / valleys

	vector<vector<vec3>> curveRidgeValleys;
	if (!FileUtil::existsfile(rvName)) {
		return error("run single model RV first");
	} else {
		if (!CurveUtil::loadCurves(rvName, curveRidgeValleys)) return false;
	}

	vector<vector<vec3>> ridgeValleySamples(0);
	if (!CurveUtil::sampleLines(curveRidgeValleys, ridgeValleySamples, visualSampleRadius)) return false;

	// load boundaries

	vector<vector<vec3>> curveBoundaries;
	if (!FileUtil::existsfile(bName)) {
		return error("run single model B first");
	} else {
		if (!CurveUtil::loadCurves(bName, curveBoundaries)) return false;
	}

	vector<vector<vec3>> boundarySamples(0);
	if (!CurveUtil::sampleLines(curveBoundaries, boundarySamples, visualSampleRadius)) return false;

	// load contours from multiple views

	vector<vector<vector<vec3>>> curveContours(numViewPoints, vector<vector<vec3>>(0));
	if (!FileUtil::existsfile(cName)) {
		return error("run single model BC first");
	} else {
		ifstream curveFile(cName, ios::binary);
		int numViews;
		curveFile.read((char*)&numViews, sizeof(numViews));
		if (numViews != numViewPoints) return error("incompatible C file");
		for (int viewID = 0; viewID < numViewPoints; viewID++) {
			vector<vector<vec3>> &curveViewContours = curveContours[viewID];
			if (!CurveUtil::loadCurves(curveFile, curveViewContours)) return false;
		}
		curveFile.close();
	}

	// snap ridges/valleys/boundaries/contours to support lines

	vector<vector<vec3>> lineRidgeValleys;
	if (true) {
#ifdef OUTPUT_PROGRESS
		cout << "Snapping ridges & valleys" << endl;
#endif
		CurveLineSnapping cls(curveRidgeValleys, curveSupportLines);
		if (!cls.process()) return false;
		if (!cls.output(lineRidgeValleys)) return false;
		if (!CurveUtil::saveCurves(snapRVName, lineRidgeValleys)) return false;
	}

	vector<vector<vec3>> lineBoundaries;
	if (true) {
#ifdef OUTPUT_PROGRESS
		cout << "Snapping boundaries" << endl;
#endif
		CurveLineVoting clv(curveBoundaries, curveSupportLines);
		if (!clv.process()) return false;
		if (!clv.output(lineBoundaries)) return false;
		if (!CurveUtil::saveCurves(snapBName, lineBoundaries)) return false;
	}

	vector<vector<vector<vec3>>> lineContours(numViewPoints);
	ofstream snapCFile(snapCName, ios::binary);
	snapCFile.write((char*)&numViewPoints, sizeof(numViewPoints));
	for (int viewID = 0; viewID < numViewPoints; viewID++) {
#ifdef OUTPUT_PROGRESS
		cout << "\rSnapping contours " << (viewID+1) << " / " << numViewPoints << "           ";
#endif
		vector<vector<vec3>> &lineViewContours = lineContours[viewID];

		CurveLineVoting clv(curveContours[viewID], curveSupportLines);
		if (!clv.process()) return false;
		if (!clv.output(lineViewContours)) return false;
		if (!CurveUtil::saveCurves(snapCFile, lineViewContours)) return false;
	}
#ifdef OUTPUT_PROGRESS
	cout << endl;
#endif
	snapCFile.close();

	// visualize ridges/valleys/boundaries/contours

	vec3i ridgeValleyColor(0, 255, 0);
	vec3i boundaryColor(255, 0, 0);
	vec3i contourColor(0, 0, 255);

	for (int id = 0; id < (int)visualViewList.size(); id++) {
#ifdef OUTPUT_PROGRESS
		cout << "\rVisualizing view " << (id + 1) << " / " << (int)visualViewList.size() << "           ";
#endif

		int visID = visualViewList[id];
		string allVName = allVPrefix + "-" + to_string(id) + ".ply";
		string snapVisName = snapVisPrefix + "-" + to_string(id) + ".ply";

		vector<vector<vec3>> contourSamples(0);
		if (!CurveUtil::sampleLines(curveContours[visID], contourSamples, visualSampleRadius)) return false;

		PlyExporter peAll;
		for (auto &sample : boundarySamples) {
			if (!peAll.addPoint(&sample, 0, cml::identity_4x4(), boundaryColor)) return false;
		}
		for (auto &sample : ridgeValleySamples) {
			if (!peAll.addPoint(&sample, 0, cml::identity_4x4(), ridgeValleyColor)) return false;
		}
		for (auto &sample : contourSamples) {
			if (!peAll.addPoint(&sample, 0, cml::identity_4x4(), contourColor)) return false;
		}		
		if (!peAll.output(allVName)) return false;
		
		vector<vector<vec3>> lineAllSnap(0);
		lineAllSnap.insert(lineAllSnap.end(), lineBoundaries.begin(), lineBoundaries.end());
		lineAllSnap.insert(lineAllSnap.end(), lineRidgeValleys.begin(), lineRidgeValleys.end());
		lineAllSnap.insert(lineAllSnap.end(), lineContours[visID].begin(), lineContours[visID].end());
		if (!CurveUtil::removeDuplicateLines(lineAllSnap)) return false;
		if (!CurveUtil::visualizeCurves(snapVisName, lineAllSnap, 0)) return false;
	}
#ifdef OUTPUT_PROGRESS
	cout << endl;
#endif

	return true;
}

bool PipelineCurveIO::error(string s) {
	cout << "Error: " << s << endl;
	return false;
}