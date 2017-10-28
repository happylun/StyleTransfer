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

#pragma once

// definition of template member function

// DO NOT INCLUDE THIS FILE EXCEPT IN FILE "PlyLoader.h"

#ifdef PLY_LOADER_TEMPLATE_DEFINITION

template<class T> bool PlyLoader::loadSample(string fileName, vector<T> *samples, vector<string> *comments) {

	if(!samples) return false;

	int numSamples = 0;

	samples->clear();

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
			ss >> numSamples;
		}
		if( line.find("end_header") != string::npos) {
			break;
		}
	}

	for(int i=0; i<numSamples; i++) {
		T sample;
		file.read((char*)&sample, sizeof(sample));
		samples->push_back(sample);
	}
	file.close();

	return true;
}

#endif