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

#include "FeatureLightFieldDescriptors.h"

#include <fstream>
#include <iostream>

#include "LFD/LightFieldDescriptor.h"

using namespace StyleSynthesis;

bool FeatureLightFieldDescriptors::mInitialized = false;

bool FeatureLightFieldDescriptors::init() {

	if (!mInitialized) {
		if (!LFD::LightFieldDescriptor::init()) {
			cout << "Error: cannot initialize LFD" << endl;
			return false;
		}
		mInitialized = true;
	}
	return mInitialized;
}

bool FeatureLightFieldDescriptors::finish() {

	if (mInitialized) {
		if (!LFD::LightFieldDescriptor::finish()) {
			cout << "Error: cannot shutdown LFD" << endl;
			return false;
		}
		mInitialized = false;
	}
	return !mInitialized;
}

FeatureLightFieldDescriptors::FeatureLightFieldDescriptors(TTriangleMesh *mesh, vector<double> *features) {

	mpMesh = mesh;
	mpSample = 0;
	mpFeatures = features;
	mpLFD = new LFD::LightFieldDescriptor;

	if (!mInitialized) init();
}

FeatureLightFieldDescriptors::FeatureLightFieldDescriptors(TSampleSet *sample, vector<double> *features) {

	mpMesh = 0;
	mpSample = sample;
	mpFeatures = features;
	mpLFD = new LFD::LightFieldDescriptor;

	if (!mInitialized) init();
}

FeatureLightFieldDescriptors::~FeatureLightFieldDescriptors() {

	if (mpLFD) delete mpLFD;
}

bool FeatureLightFieldDescriptors::calculate() {

	if (mpMesh) {
		vector<float> vecVertices;
		vector<int> vecFaces;

		vecVertices.reserve(mpMesh->amount * 3);
		for (vec3 v : mpMesh->positions) {
			for (int j = 0; j < 3; j++) vecVertices.push_back(v[j]);
		}

		vecFaces.reserve(mpMesh->indices.size() * 3);
		for (vec3i vi : mpMesh->indices) {
			for (int j = 0; j < 3; j++) vecFaces.push_back(vi[j]);
		}

		if (!mpLFD->renderMesh(vecVertices, vecFaces)) return false;
	}

	if (mpSample) {
		vector<float> vecVertices;

		vecVertices.reserve(mpSample->amount * 3);
		for (vec3 v : mpSample->positions) {
			for (int j = 0; j < 3; j++) vecVertices.push_back(v[j]);
		}

		if (!mpLFD->renderPointCloud(vecVertices, mpSample->radius)) return false;
	}

	if (!mpLFD->calculate(*mpFeatures)) return false;

	return true;
}

bool FeatureLightFieldDescriptors::visualize(string fileName) {

	return mpLFD->visualize(fileName);
}

bool FeatureLightFieldDescriptors::compareFeatures(vector<double> &feature1, vector<double> &feature2, vector<double> &distances) {

	return LFD::LightFieldDescriptor::compare(feature1, feature2, distances);
}