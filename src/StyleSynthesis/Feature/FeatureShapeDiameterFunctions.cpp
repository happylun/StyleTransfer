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

#include "FeatureShapeDiameterFunctions.h"

#include <fstream>
#include <algorithm>
#include <set>

#include "Utility/PlyExporter.h"

#include "Data/StyleSynthesisConfig.h"

#include "Feature/FeatureUtil.h"

using namespace StyleSynthesis;

FeatureShapeDiameterFunctions::FeatureShapeDiameterFunctions(TSampleSet *samples, TTriangleMesh *mesh, vector<double> *features) {

	mpSamples = samples;
	mpFeatures = features;

	// build KD tree

	mTreeData.resize(mesh->indices.size());

#pragma omp parallel for
	for (int faceID = 0; faceID<(int)mesh->indices.size(); faceID++) {
		vec3i idx = mesh->indices[faceID];
		G3D::Vector3 v0(mesh->positions[idx[0]].data());
		G3D::Vector3 v1(mesh->positions[idx[1]].data());
		G3D::Vector3 v2(mesh->positions[idx[2]].data());
		mTreeData[faceID].set(TKDT::NamedTriangle(v0, v1, v2, faceID));
	}

	mpMeshTree = new TKDTree(mTreeData.begin(), mTreeData.end());
}

FeatureShapeDiameterFunctions::FeatureShapeDiameterFunctions(TSampleSet *samples, TKDTree *meshTree, vector<double> *features) {

	mpSamples = samples;
	mpMeshTree = meshTree;
	mpFeatures = features;
	mTreeData.clear();
}

FeatureShapeDiameterFunctions::~FeatureShapeDiameterFunctions() {

	if (!mTreeData.empty() && mpMeshTree) delete mpMeshTree;
}

bool FeatureShapeDiameterFunctions::calculate() {
	
	if (!calculateSDF()) return false;
	if (!normalizeSDF()) return false;

	return true;
}

bool FeatureShapeDiameterFunctions::calculateSDF() {

	float coneRange = cos(cml::rad(20.0f)); // UNDONE: param ray cone angle range
	int numRays = 30;
	float eps = mpSamples->radius * 0.001f;
	double pruneStd = 1.0;

	mpFeatures->resize(mpSamples->amount);
#pragma omp parallel for
	for (int sampleID = 0; sampleID < mpSamples->amount; sampleID++) {

		vec3 p = mpSamples->positions[sampleID];
		vec3 n = mpSamples->normals[sampleID];

		// calculate diameter by sampling rays
		vector<double> rayLengths(0);
		vector<double> rayWeights(0);
		for (int rayID = 0; rayID < numRays; rayID++) {

			// random sampling on ray direction within cone range
			vec3 rayDir;
			while (true) {
				double r1 = cml::random_unit();
				double r2 = cml::random_unit();
				rayDir = vec3d(
					2.0 * cos(cml::constantsd::two_pi()*r1) * cml::sqrt_safe(r2*(1 - r2)),
					2.0 * sin(cml::constantsd::two_pi()*r1) * cml::sqrt_safe(r2*(1 - r2)),
					1.0 - 2.0*r2); // random direction on unit sphere
				float cosDN = cml::dot((vec3)rayDir, n);
				if (cosDN > coneRange || cosDN < -coneRange) { // retain ray direction within cone
					if (cosDN > 0) rayDir = -rayDir; // flip ray direction
					break;
				}
			}

			// calculate ray intersection
			vec3 rayOrigin = p - n * eps; // add eps for offset
			Thea::Ray3 ray(G3D::Vector3(rayOrigin.data()), G3D::Vector3(rayDir.data()));
			auto hitResult = mpMeshTree->rayStructureIntersection(ray);
			double dist = hitResult.getTime();
			if (dist > 0) { // intersects with mesh

				//auto hitNormal = hitResult.getNormal();
				//vec3 rayNormal((float)hitNormal[0], (float)hitNormal[1], (float)hitNormal[2]);
				//if (cml::dot(rayDir, rayNormal) > 0) { // skip intersection with 'outside' of the mesh
				{ // no skip: in case of two sided mesh
					rayLengths.push_back(dist);
					rayWeights.push_back(cml::dot(-rayDir, n));
				}
			}
		}

		if (rayLengths.empty()) {
			// no intersection at all
			(*mpFeatures)[sampleID] = 0; // any special handling ?
			continue;
		}

		// get statistics
		int numLen = (int)rayLengths.size();
		double avgLen = 0;
		double stdLen = 0;
		double medianLen = 0;
		double maxLen = 0;
		{
			for (double len : rayLengths) avgLen += len;
			avgLen /= numLen;
			for (double len : rayLengths) stdLen += cml::sqr(len - avgLen);
			stdLen = cml::sqrt_safe(stdLen / numLen);
			vector<double> sortedLen(rayLengths.begin(), rayLengths.end());
			nth_element(sortedLen.begin(), sortedLen.begin() + numLen / 2, sortedLen.end());
			medianLen = sortedLen[numLen / 2];
			for (double len : rayLengths) maxLen = max(maxLen, len);
		}

		// prune rays
		vector<bool> flagRays(numLen, true);
		for (int rayID = 0; rayID < numLen; rayID++) {
			double len = rayLengths[rayID];
			if (len < medianLen - stdLen*pruneStd || len > medianLen + stdLen*pruneStd) {
				flagRays[rayID] = false;
			}
		}

		// export SDF value as weighted sum
		double totalValue = 0;
		double totalWeight = 0;
		for (int rayID = 0; rayID < numLen; rayID++) {
			if (flagRays[rayID]) {
				totalValue += rayLengths[rayID] * rayWeights[rayID];
				totalWeight += rayWeights[rayID];
			}
		}
		if (totalWeight == 0) {
			cout << "Error: zero weights for SDF" << endl;
		}
		(*mpFeatures)[sampleID] = totalValue / totalWeight;
	}

	return true;
}

bool FeatureShapeDiameterFunctions::smoothSDF() {

	// build neighboring graph

	SKDTreeData treeData(mpSamples->amount);
#pragma omp parallel for
	for (int sampleID = 0; sampleID<mpSamples->amount; sampleID++) {
		vec3 v = mpSamples->positions[sampleID];
		treeData[sampleID] = SKDT::NamedPoint(v[0], v[1], v[2], (size_t)sampleID);
	}
	SKDTree tree(treeData.begin(), treeData.end());

	vector<set<int>> graph;
	graph.clear();
	graph.resize(mpSamples->amount);
	for (auto &it : graph) it.clear();

	double nbRange = mpSamples->radius * 3.5;
#pragma omp parallel for
	for (int sampleID = 0; sampleID < mpSamples->amount; sampleID++) {
		vec3 v = mpSamples->positions[sampleID];
		SKDT::NamedPoint queryPoint(v[0], v[1], v[2]);
		Thea::BoundedSortedArray<SKDTree::Neighbor> queryResult(10);
		tree.kClosestElements<Thea::MetricL2>(queryPoint, queryResult, nbRange);
		for (int j = 0; j < queryResult.size(); j++) {
			int neighborID = (int)tree.getElements()[queryResult[j].getIndex()].id;
#pragma omp critical
			{
				graph[sampleID].insert(neighborID);
				graph[neighborID].insert(sampleID);
			}
		}
	}
	
	// smoothing

	double gaussianDeviation = mpSamples->radius * 1.0; // UNDONE: param neighbor contribution Gaussian sigma
	double clampDifference = mpSamples->radius * 10.0; // UNDONE: param maximum neighbor distance
	int numIterations = 10;
	for (int iterID = 0; iterID < numIterations; iterID++) {
		vector<double> prevFeatures(mpFeatures->begin(), mpFeatures->end());
#pragma omp parallel for
		for (int sampleID = 0; sampleID < mpSamples->amount; sampleID++) {
			vec3 sampleP = mpSamples->positions[sampleID];
			double sampleV = prevFeatures[sampleID];

			double totalValue = sampleV;
			double totalWeight = 1.0;
			for (int neighborID : graph[sampleID]) {
				vec3 neighborP = mpSamples->positions[neighborID];
				double neighborV = prevFeatures[neighborID];

				double weight = exp(-(sampleP - neighborP).length_squared() / cml::sqr(gaussianDeviation));
				if (fabs(neighborV - sampleV) > clampDifference) weight = 0;

				totalValue += neighborV * weight;
				totalWeight += weight;
			}

			(*mpFeatures)[sampleID] = totalValue / totalWeight;
		}
	}

	return true;
}

bool FeatureShapeDiameterFunctions::normalizeSDF() {

	if (mpFeatures->empty()) {
		cout << "Error: no computed features" << endl;
		return false;
	}

	double minValue = DBL_MAX;
	double maxValue = -DBL_MAX;

	for (double v : (*mpFeatures)) {
		minValue = min(minValue, v);
		maxValue = max(maxValue, v);
	}

	if (minValue >= maxValue) { // no need to normalize
		for (double &v : (*mpFeatures)) v = 1.0;
		return true;
	}

	double alpha = 4.0; // UNDONE: param normalization power parameter

#pragma omp parallel for
	for (int id = 0; id < (int)mpFeatures->size(); id++) {
		double &v = (*mpFeatures)[id];
		if(v) v = log((v - minValue) / (maxValue - minValue) * alpha + 1) / log(alpha + 1);
	}

	return true;
}

bool FeatureShapeDiameterFunctions::visualize(string fileName) {

	PlyExporter pe;

	vector<vec3i> vColors;
	for (int sampleID = 0; sampleID < mpSamples->amount; sampleID++) {
		double v = (*mpFeatures)[sampleID];
		vColors.push_back(FeatureUtil::colorMapping(v));
	}

	if (!pe.addPoint(&mpSamples->positions, &mpSamples->normals, &vColors)) return false;
	if (!pe.output(fileName)) return false;

	return true;
}

bool FeatureShapeDiameterFunctions::compareFeatures(vector<double> &feature1, vector<double> &feature2, vector<double> &distance) {

	distance.clear();

	int histBins[] = { 16, 32, 64, 128 };
	for (int k = 0; k < 4; k++) {
		vector<double> hist1, hist2;
		if (!FeatureUtil::computeHistogram(feature1, hist1, histBins[k])) return false;
		if (!FeatureUtil::computeHistogram(feature2, hist2, histBins[k])) return false;
		double dist = FeatureUtil::computeEMD(hist1, hist2);
		distance.push_back(dist);
	}

	return true;
}
