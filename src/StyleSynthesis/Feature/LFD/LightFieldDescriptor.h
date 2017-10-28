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

#include <vector>
#include <string>
#include <iostream>

using namespace std;

namespace LFD {
	class LightFieldDescriptor {

	public:

		LightFieldDescriptor();
		~LightFieldDescriptor();

	public:

		static bool init();
		static bool finish();
		static bool compare(vector<double> &inDescriptor1, vector<double> &inDescriptor2, vector<double> &outDistance);

	public:

		bool renderMesh(vector<float> &inVertices, vector<int> &inFaces);
		bool renderPointCloud(vector<float> &inPoints, float sampleRadius);
		bool calculate(vector<double> &outDescriptors);
		bool visualize(string fileName);

	private:

		static bool normalize(vector<float> &vertices, vector<int> &faces);
		static bool normalize(vector<float> &points, float &scale);

		static bool validateImage(vector<unsigned char> &image);
		static bool findCenter(vector<unsigned char> &inImage, pair<double,double> &outCenter);

		inline static bool error(string s) { cout << "Error: " << s << endl; return false; }

	private:

		vector<unsigned char> **mpImage;
		pair<double, double> **mpCenter;
	};
}