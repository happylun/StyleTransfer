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

#include "PipelineMeshIO.h"

#include <iostream>
#include <fstream>

#include "Mesh/MeshUtil.h"

#include "Data/StyleSynthesisConfig.h"

using namespace std;
using namespace StyleSynthesis;

bool PipelineMeshIO::process() {

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
		cout << "Processing mesh " << meshName << endl;
		if (!runSingleModel(meshName)) return error("single model error");

		//system("pause");
	}

	return true;
}

bool PipelineMeshIO::runSingleModel(string meshName) {

	// names...

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string rawMeshFileName = datasetPrefix + "raw-mesh/" + meshName + ".ply";
	string outMeshFileName = datasetPrefix + "scaled-mesh/" + meshName + ".ply";
	string rawPartFolderName = datasetPrefix + "raw-part/" + meshName + "/";
	string outPartFolderName = datasetPrefix + "scaled-part/" + meshName + "/";

	if (!FileUtil::makedir(outMeshFileName)) return false;
	if (!FileUtil::makedir(outPartFolderName)) return false;

	if (FileUtil::existsfile(outMeshFileName)) return true; // early quit

	// data

	TTriangleMesh mesh;
	vector<TTriangleMesh> parts;

	// load mesh and parts
	
	if (!MeshUtil::loadMesh(rawMeshFileName, mesh)) return false;
	
	int numParts = 0;
	while (true) {
		string rawPartName = rawPartFolderName + to_string(numParts) + ".ply";
		if (!FileUtil::existsfile(rawPartName)) break;
		parts.push_back(TTriangleMesh());
		if (!MeshUtil::loadMesh(rawPartName, parts.back())) return false;
		numParts++;
	}

	// make sure mesh and parts are on the same scale

	if (true) {
		// normalize mesh
		if (true) {
			vec3 bbMin, bbMax;
			if (!MeshUtil::computeAABB(mesh, bbMin, bbMax)) return false;
			vec3 bbSize = bbMax - bbMin;
			vec3 center = (bbMin + bbMax) / 2;
			float scale = 1.0f / max(bbSize[0], max(bbSize[1], bbSize[2]));
			for (vec3 &v : mesh.positions) {
				v = (v - center) * scale;
			}
		}
		// normalize parts
		if (true) {
			vec3 bbMin, bbMax;
			if (!MeshUtil::computeAABB(parts[0], bbMin, bbMax)) return false;
			for (int partID = 1; partID < numParts; partID++) {
				vec3 partBBMin, partBBMax;
				if (!MeshUtil::computeAABB(parts[partID], partBBMin, partBBMax)) return false;
				bbMin.minimize(partBBMin);
				bbMax.maximize(partBBMax);
			}
			vec3 bbSize = bbMax - bbMin;
			vec3 center = (bbMin + bbMax) / 2;
			float scale = 1.0f / max(bbSize[0], max(bbSize[1], bbSize[2]));
			for (int partID = 0; partID < numParts; partID++) {
				for (vec3 &v : parts[partID].positions) {
					v = (v - center) * scale;
				}
			}
		}
	}

	// compute transformation

	matrix3f meshRotate;
	meshRotate.identity();
	vec3 meshTranslate(0.0f, 0.0f, 0.0f);
	float meshScale = 1.0f;

	if (meshName[0] == 'F') {
		// furniture
		// drop to ground; scale shortest bounding box axis length to one

		vec3 bbMin, bbMax;
		if (!MeshUtil::computeAABB(mesh, bbMin, bbMax)) return false;
		vec3 bbLen = bbMax - bbMin;
		float minLen = bbLen[0];
		for (int k = 1; k < 3; k++) {
			minLen = min(minLen, bbLen[k]);
		}
		meshTranslate = vec3(0.0f, -bbMin[1], 0.0f);
		meshScale = 1.0f / minLen;
	} else if (meshName[0] == 'L') {
		// lamp
		// scale bounding box diagonal length of largest part to one; move center to origin

		float maxMinLen = 0;
		for (int partID = 0; partID < numParts; partID++) {
			vec3 bbMin, bbMax;
			if (!MeshUtil::computeAABB(parts[partID], bbMin, bbMax)) return false;
			vec3 bbLen = bbMax - bbMin;
			//float volume = (bbLen[0] * bbLen[1] * bbLen[2]);
			float minLen = min(bbLen[0], min(bbLen[1], bbLen[2]));
			if (minLen > maxMinLen) {
				maxMinLen = minLen;
				meshScale = 1.0f / bbLen.length();
			}
		}
		vec3 meshBBMin, meshBBMax;
		if (!MeshUtil::computeAABB(mesh, meshBBMin, meshBBMax)) return false;
		meshTranslate = -(meshBBMin + meshBBMax) / 2;
	} else if (meshName[0] == 'S') {
		// cutlery
		// rotate up; scale longest bounding box axis length to one; move center to origin

		vec3 bbMin, bbMax;
		if (!MeshUtil::computeAABB(mesh, bbMin, bbMax)) return false;
		vec3 bbLen = bbMax - bbMin;
		float maxLen = bbLen[0];
		for (int k = 1; k < 3; k++) {
			maxLen = max(maxLen, bbLen[k]);
		}
		meshTranslate = -(bbMin + bbMax) / 2;
		meshScale = 3.0f / maxLen; // HACK: scale a bit so that we have enough sample points for features...

		cml::matrix_rotation_world_x(meshRotate, cml::rad(90.0f));
	}

	// transform mesh

	for (vec3 &v : mesh.positions) {
		v = (meshRotate * v + meshTranslate) * meshScale;
	}
	if (!MeshUtil::saveMesh(outMeshFileName, mesh)) return false;

	// transform part

	for (int partID = 0; partID < numParts; partID++) {
		string outPartName = outPartFolderName + to_string(partID) + ".ply";
		TTriangleMesh &part = parts[partID];
		for (vec3 &v : part.positions) {
			v = (meshRotate * v + meshTranslate) * meshScale;
		}
		if (!MeshUtil::saveMesh(outPartName, part)) return false;
	}

	return true;
}

bool PipelineMeshIO::error(string s) {
	cout << "Error: " << s << endl;
	return false;
}