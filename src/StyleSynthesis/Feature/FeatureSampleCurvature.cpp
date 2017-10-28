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

#include "FeatureSampleCurvature.h"

#include <fstream>
#include <cmath>
#include <set>

#include <Eigen/Eigen>

#include "Utility/PlyExporter.h"

#include "Data/StyleSynthesisConfig.h"

#include "Sample/SampleUtil.h"
#include "Segment/SegmentUtil.h"
#include "Feature/FeatureUtil.h"

using namespace StyleSynthesis;

FeatureSampleCurvature::FeatureSampleCurvature(TSampleSet *samples, vector<TCurvature> *curvatures) {

	mpSamples = samples;
	mpCurvatures = curvatures;
}

FeatureSampleCurvature::~FeatureSampleCurvature() {
}

bool FeatureSampleCurvature::calculate() {

	if (!buildGeodesicGraph()) return false;
	if (!getPatchNeighbors()) return false;
	if (!runPatchFitting()) return false;

	return true;
}

bool FeatureSampleCurvature::buildGeodesicGraph() {

	cout << "Building graph..." << endl;

	vector<bool> flags;
	if (!SegmentUtil::buildKNNGraph(*mpSamples, mGraph, flags)) return false;

	return true;
}

bool FeatureSampleCurvature::getPatchNeighbors() {

	cout << "Finding neighbors..." << endl;

	// find points within geodesic radius

	float patchRadius = (float)(mpSamples->radius * 2.0); // UNDONE: param geodesic radius

	mPatch.resize(mpSamples->amount);
#pragma omp parallel for
	for (int sampleID = 0; sampleID < mpSamples->amount; sampleID++) {

		vec3 sampleP = mpSamples->positions[sampleID];
		vec3 sampleN = mpSamples->normals[sampleID];

		// not optimal in complexity, but correct and fast enough
		set<int> pointSet;
		map<int, float> pointDists; // sample ID : min distance
		pointDists[sampleID] = 0;
		vector<pair<int, float>> queue; // (sample ID, distance) : queue length
		queue.push_back(make_pair(sampleID, 0.0f));
		int head = 0;
		while (head < (int)queue.size()) {
			int nodeID = queue[head].first;
			vec3 nodeP = mpSamples->positions[nodeID];
			vec3 nodeN = mpSamples->normals[nodeID];
			float nodeDist = queue[head].second;
			for (int neighborID : mGraph[nodeID]) {
				vec3 neighborP = mpSamples->positions[neighborID];
				vec3 neighborN = mpSamples->normals[neighborID];
				float edgeLen = (nodeP-neighborP).length();
				float neighborDist = nodeDist + edgeLen;
				if (neighborDist > patchRadius || cml::dot(nodeN, neighborN) < -0.01f) continue;
				auto &found = pointDists.find(neighborID);
				if (found == pointDists.end() || found->second > neighborDist) {
					pointDists[neighborID] = neighborDist;
					queue.push_back(make_pair(neighborID, neighborDist));
					pointSet.insert(neighborID);
				}
			}
			head++;
		}
		mPatch[sampleID].assign(pointSet.begin(), pointSet.end());
	}	

	return true;
}

bool FeatureSampleCurvature::runPatchFitting() {

	cout << "Fitting patch..." << endl;

	mpCurvatures->resize(mpSamples->amount);

	const double sigma = cml::sqr(mpSamples->radius * 2.0); // UNDONE: param neighbor gaussian radius
	const double maxKMagnitude = 1 / (mpSamples->radius * 2.0); // UNDONE: param maximum magnitude for curvature value

	// estimate curvature tensor from neighbors

#pragma omp parallel for
	for (int sampleID = 0; sampleID < mpSamples->amount; sampleID++) {

		vec3d position = mpSamples->positions[sampleID];
		vec3d normal = mpSamples->normals[sampleID];

		// establish local frame

		vec3d localCenter = position;
		vec3d localZ = normal;
		vec3d localX = cml::normalize(cml::cross(normal, fabs(normal[0])<fabs(normal[1]) ? vec3d(1.0, 0.0, 0.0) : vec3d(0.0, 1.0, 0.0)));
		vec3d localY = cml::normalize(cml::cross(localZ, localX));

		if (mPatch[sampleID].size() < 3) { // cannot reliably estimate curvature
			TCurvature &curvature = (*mpCurvatures)[sampleID];
			curvature.k1 = 0;
			curvature.k2 = 0;
			curvature.d1 = (vec3)localX;
			curvature.d2 = (vec3)localY;
			continue;
		}

		matrix3d localFrame; // world CS to local CS
		cml::matrix_set_basis_vectors(localFrame, localX, localY, localZ);
		localFrame.transpose();

		// build linear system

		int kNeighbors = (int)mPatch[sampleID].size();
		Eigen::MatrixXd matU(3 * kNeighbors, 7);
		Eigen::VectorXd matZ(3 * kNeighbors, 1);

		for (int j = 0; j < kNeighbors; j++) {

			int neighborID = mPatch[sampleID][j];

			vec3d neighborP = mpSamples->positions[neighborID];
			vec3d neighborN = mpSamples->normals[neighborID];
			vec3d localP = localFrame * (neighborP - localCenter);
			vec3d localN = cml::normalize(localFrame * neighborN);

			double x = localP[0];
			double y = localP[1];
			double z = localP[2];
			double a = localN[0];
			double b = localN[1];
			double c = localN[2];
			double w = exp(-(x*x + y*y) / sigma); // Gaussian weight

			if (!(w>0 && w <= 1)) {
				cout << "Error: incorrect weight" << endl;
			}
			if (!(c>=-1 && c <= 1)) {
				cout << "Error: incorrect c" << endl;
			}
			if (fabs(c) < 1e-7) {
				//cout << "Error: normal inconsistent with local frame" << endl;
				c = sqrt(a*a + b*b) * (c < 0 ? -1 : 1);
			}

			matU.row(j * 3) << 0.5f*x*x, x*y, 0.5f*y*y, x*x*x, x*x*y, x*y*y, y*y*y;
			matZ.row(j * 3) << z;

			matU.row(j * 3 + 1) << x, y, 0, 3*x*x, 2*x*y, y*y, 0;
			matZ.row(j * 3 + 1) << -a/c;

			matU.row(j * 3 + 2) << 0, x, y, 0, x*x, 2*x*y, 3*y*y;
			matZ.row(j * 3 + 2) << -b/c;

			matU.middleRows(j*3, 3) *= w;
			matZ.middleRows(j*3, 3) *= w;
		}

		if (!matU.allFinite() || !matZ.allFinite()) {
			cout << "Error: matrix is invalid!" << endl;
			cout << matU.allFinite() << ", " << matZ.allFinite() << endl;
		}

		// solve linear system: U * X = Z
		Eigen::VectorXd matX = matU.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(matZ);

		// curvature tensor
		Eigen::Matrix2d matK;
		matK << -matX(0), -matX(1), -matX(1), -matX(2);

		// extract principal curvature and principal curvature direction

		Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigenSolver(matK);
		if (eigenSolver.info() != Eigen::Success) {
			cout << "Error: fail to solve for principal curvatures" << endl;
		}

		TCurvature &curvature = (*mpCurvatures)[sampleID];
		curvature.k1 = (float)(eigenSolver.eigenvalues()(1));
		curvature.k2 = (float)(eigenSolver.eigenvalues()(0));
		curvature.d1 = (vec3)(cml::normalize(localX * eigenSolver.eigenvectors()(0, 1) + localY * eigenSolver.eigenvectors()(1, 1)));
		curvature.d2 = (vec3)(cml::normalize(localX * eigenSolver.eigenvectors()(0, 0) + localY * eigenSolver.eigenvectors()(1, 0)));

		if (curvature.k1 != curvature.k1 || curvature.k2 != curvature.k2) {
			cout << "Error: NaN data in curvature estimation" << endl;
		}

		// clamp curvature value to prevent outlier
		curvature.k1 = cml::clamp(curvature.k1, (float)-maxKMagnitude, (float)maxKMagnitude);
		curvature.k2 = cml::clamp(curvature.k2, (float)-maxKMagnitude, (float)maxKMagnitude);

		// re-orient curvature direction to make {n, d1, d2} a right-hand-side coordinate system
		if (cml::dot(cml::cross(curvature.d1, curvature.d2), (vec3)normal) < 0) curvature.d2 = -curvature.d2;
	}

	return true;
}

bool FeatureSampleCurvature::visualize(string fileName) {

	PlyExporter pe;

	// export curvature direction as thin rects

	vector<vec3i> vI; // indices of a rect
	vI.push_back(vec3i(0, 1, 2));
	vI.push_back(vec3i(1, 3, 2));

	int numPoints = mpSamples->amount;
	float radius = mpSamples->radius;

	float maxAbsK1 = 0;
	float maxAbsK2 = 0;
	vector<float> allAbsK1;
	vector<float> allAbsK2;
	for (TCurvature &curvature : (*mpCurvatures)) {
		allAbsK1.push_back(fabs(curvature.k1));
		allAbsK2.push_back(fabs(curvature.k2));
	}
	int n = (int)(numPoints * 0.8);
	nth_element(allAbsK1.begin(), allAbsK1.begin() + n, allAbsK1.end());
	nth_element(allAbsK2.begin(), allAbsK2.begin() + n, allAbsK2.end());
	maxAbsK1 = allAbsK1[n];
	maxAbsK2 = allAbsK2[n];

	for (int pointID = 0; pointID < numPoints; pointID++) {

		vec3 position = mpSamples->positions[pointID];
		TCurvature &curvature = (*mpCurvatures)[pointID];

		vec3 p0 = position + ( curvature.d1.normalize() + curvature.d2.normalize() * 0.2f) * radius * 0.5f;
		vec3 p1 = position + (-curvature.d1.normalize() + curvature.d2.normalize() * 0.2f) * radius * 0.5f;
		vec3 p2 = position + ( curvature.d1.normalize() - curvature.d2.normalize() * 0.2f) * radius * 0.5f;
		vec3 p3 = position + (-curvature.d1.normalize() - curvature.d2.normalize() * 0.2f) * radius * 0.5f;

		vector<vec3> vP;
		vP.push_back(p0);
		vP.push_back(p1);
		vP.push_back(p2);
		vP.push_back(p3);

		float c = cml::clamp(curvature.k1/maxAbsK1, -1.0f, 1.0f);
		int c1 = max(0, (int)(c * 255));
		int c2 = max(0, (int)(-c * 255));
		vec3i color = vec3i(255-c2, 255-c1-c2, 255-c1);

		if (!pe.addMesh(&vI, &vP, 0, cml::identity_4x4(), color)) return false;
	}

	if (!pe.output(fileName)) return false;

	return true;
}

bool FeatureSampleCurvature::getMetrics(TCurvature inCurvature, vector<double> &outMetrics) {

	double k1 = inCurvature.k1;
	double k2 = inCurvature.k2;
	double k1p = fabs(k1) > fabs(k2) ? k1 : k2;
	double k2p = fabs(k1) > fabs(k2) ? k2 : k1;

	// first 6 metrics can be either positive or negative
	// the remaining 7 metrics should be positive

	outMetrics.clear();
	outMetrics.push_back(k1);                        // max by value
	outMetrics.push_back(k2);                        // min by value
	outMetrics.push_back(k1p);                       // max by magnitude
	outMetrics.push_back(k2p);                       // min by magnitude
	outMetrics.push_back((k1 + k2) / 2);             // mean
	outMetrics.push_back(k1*k2);                     // Gaussian
	outMetrics.push_back(fabs(k1));                  // abs max by value
	outMetrics.push_back(fabs(k2));                  // abs min by value
	outMetrics.push_back(fabs(k1p));                 // abs max by magnitude
	outMetrics.push_back(fabs(k2p));                 // abs min by magnitude
	outMetrics.push_back(fabs(k1 + k2) / 2);         // abs mean
	outMetrics.push_back(fabs(k1*k2));               // abs Gaussian
	outMetrics.push_back((fabs(k1) + fabs(k2)) / 2); // mean magnitude

	return true;
}

bool FeatureSampleCurvature::compareFeatures(vector<TCurvature> &curvature1, vector<TCurvature> &curvature2, vector<double> &distance) {

	int n1 = (int)curvature1.size();
	int n2 = (int)curvature2.size();
	int dim = 0;

	Eigen::MatrixXd feature1, feature2; // # points X dim

	for (int id = 0; id < n1; id++) {
		vector<double> metrics;
		if (!getMetrics(curvature1[id], metrics)) return false;
		if (dim == 0) {
			dim = (int)metrics.size();
			feature1.resize(n1, dim);
			feature2.resize(n2, dim);
		}
		for (int d = 0; d < dim; d++) feature1(id, d) = metrics[d];
	}

	for (int id = 0; id < n2; id++) {
		vector<double> metrics;
		if (!getMetrics(curvature2[id], metrics)) return false;
		if (dim != (int)metrics.size()) return false;
		for (int d = 0; d < dim; d++) feature2(id, d) = metrics[d];
	}

	int histBins[] = { 16, 32, 64, 128 };

	distance.clear();
	for (int d = 0; d < dim; d++) {
		double *ptr1 = feature1.col(d).data();
		double *ptr2 = feature2.col(d).data();
		vector<double> metric1(ptr1, ptr1 + n1);
		vector<double> metric2(ptr2, ptr2 + n2);

		// find value range (min: 10-th percentile; max: 90-th percentile)
		vector<double> allValues(0);
		allValues.reserve(n1 + n2);
		allValues.insert(allValues.end(), metric1.begin(), metric1.end());
		allValues.insert(allValues.end(), metric2.begin(), metric2.end());
		int minN = (int)((n1 + n2)*0.1);
		int maxN = (int)((n1 + n2)*0.9);
		nth_element(allValues.begin(), allValues.begin() + minN, allValues.end());
		double minValue = allValues[minN];
		nth_element(allValues.begin() + minN, allValues.begin() + maxN, allValues.end());
		double maxValue = allValues[maxN];

		for (int k = 0; k < 4; k++) {
			vector<double> hist1, hist2;
			if (!FeatureUtil::computeHistogram(metric1, hist1, histBins[k], minValue, maxValue)) return false;
			if (!FeatureUtil::computeHistogram(metric2, hist2, histBins[k], minValue, maxValue)) return false;
			double dist = FeatureUtil::computeEMD(hist1, hist2);
			distance.push_back(dist);
		}
	}

	return true;
}

bool FeatureSampleCurvature::saveFeature(string fileName, vector<TCurvature> &curvatures) {

	ofstream outFile(fileName, ios::binary);
	if (!outFile.is_open()) {
		cout << "Error: cannot write to file " << fileName << endl;
		return false;
	}

	int numPoints = (int)curvatures.size();
	outFile.write((char*)&numPoints, sizeof(numPoints));

	for (TCurvature &curvature : curvatures) {
		outFile.write((char*)&curvature, sizeof(curvature));
	}

	return true;
}

bool FeatureSampleCurvature::loadFeature(string fileName, vector<TCurvature> &curvatures) {

	ifstream inFile(fileName, ios::binary);
	if (!inFile.is_open()) {
		cout << "Error: cannot open file " << fileName << endl;
		return false;
	}

	int numPoints;
	inFile.read((char*)&numPoints, sizeof(numPoints));
	curvatures.resize(numPoints);
	for (TCurvature &curvature : curvatures) {
		inFile.read((char*)&curvature, sizeof(curvature));
	}
	inFile.close();

	return true;
}
