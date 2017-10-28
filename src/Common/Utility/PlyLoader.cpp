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

#include "PlyLoader.h"

#include <fstream>
#include <sstream>

bool PlyLoader::loadMesh(string fileName, vector<vec3i> *indices, vector<vec3> *vertices, vector<vec3> *normals, vector<string> *comments) {

	if(!indices) return false;
	if(!vertices) return false;

	int numVertices = 0;
	int numFaces = 0;

	indices->clear();
	vertices->clear();
	if(normals) normals->clear();

	ifstream file;

	file.open(fileName, ios::binary);
	if(!file) {
		cout << "File " << fileName << " does not exist!" << endl;
		return false;
	}

	string line;
	while(true) {
		getline(file, line, '\n');
		if( comments && line.find("comment") != string::npos) {
			comments->push_back(line);
		}
		if( line.find("element vertex") != string::npos) {
			stringstream ss(line);
			string s;
			ss >> s; ss >> s;
			ss >> numVertices;
		}
		if( line.find("element face") != string::npos) {
			stringstream ss(line);
			string s;
			ss >> s; ss >> s;
			ss >> numFaces;
		}
		if( line.find("end_header") != string::npos) {
			break;
		}
	}

	for(int i=0; i<numVertices; i++) {
		vec3 vPosition;
		vec3 vNormal;
		file.read((char*)vPosition.data(), sizeof(vPosition));
		file.read((char*)vNormal.data(), sizeof(vNormal));
		vertices->push_back(vPosition);
		if(normals) normals->push_back(vNormal);
	}
	for(int i=0; i<numFaces; i++) {
		unsigned char n;
		file.read((char*)(&n), sizeof(n));
		if(n == 3) {
			vec3i vIndex;
			file.read((char*)vIndex.data(), sizeof(vIndex));
			indices->push_back(vIndex);
		} else {
			file.seekg(n*sizeof(int), ios::cur);
		}
	}
	file.close();

	return true;
}

bool PlyLoader::loadPoint(string fileName, vector<vec3> *vertices, vector<vec3> *normals, vector<string> *comments) {

	if (!vertices) return false;

	int numVertices = 0;

	vertices->clear();
	if (normals) normals->clear();
	ifstream file;

	file.open(fileName, ios::binary);
	if(!file) {
		cout << "File " << fileName << " does not exist!" << endl;
		return false;
	}

	string line;
	while(true) {
		getline(file, line, '\n');
		if( comments && line.find("comment") != string::npos) {
			comments->push_back(line);
		}
		if( line.find("element vertex") != string::npos) {
			stringstream ss(line);
			string s; ss >> s; ss >> s;
			ss >> numVertices;
		}
		if( line.find("end_header") != string::npos) {
			break;
		}
	}

	for(int i=0; i<numVertices; i++) {
		vec3 vPosition;
		vec3 vNormal;
		file.read((char*)vPosition.data(), sizeof(vPosition));
		file.read((char*)vNormal.data(), sizeof(vNormal));
		vertices->push_back(vPosition);
		if(normals) normals->push_back(vNormal);
	}
	file.close();

	return true;
}
