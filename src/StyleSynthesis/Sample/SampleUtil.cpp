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

#include "SampleUtil.h"

#include <fstream>

#include "Utility/PlyExporter.h"
#include "Utility/PlyLoader.h"
#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

bool SampleUtil::buildKdTree(vector<vec3> &samples, SKDTree &tree, SKDTreeData &treeData) {

	treeData.resize(samples.size());
	for(size_t i=0; i<samples.size(); i++) {
		vec3 &v = samples[i];
		treeData[i] = SKDT::NamedPoint(v[0], v[1], v[2], i);
	}
	tree.init(treeData.begin(), treeData.end());
	
	return true;
}

bool SampleUtil::buildKdTree(Eigen::Matrix3Xd &samples, SKDTree &tree, SKDTreeData &treeData) {

	treeData.resize(samples.cols());
#pragma omp parallel for
	for (int i = 0; i<samples.cols(); i++) {
		treeData[i] = SKDT::NamedPoint(
			(float)samples(0, i),
			(float)samples(1, i),
			(float)samples(2, i),
			(size_t)i);
	}
	tree.init(treeData.begin(), treeData.end());

	return true;
}

bool SampleUtil::findNearestNeighbors(SKDTree &tree, Eigen::Matrix3Xd &inPoints, Eigen::VectorXi &outIndices) {

	outIndices.resize(inPoints.cols());
#pragma omp parallel for
	for (int j = 0; j < inPoints.cols(); j++) {
		SKDT::NamedPoint queryPoint((float)inPoints(0, j), (float)inPoints(1, j), (float)inPoints(2, j));
		Thea::BoundedSortedArray<SKDTree::Neighbor> queryResult(1);
		tree.kClosestElements<Thea::MetricL2>(queryPoint, queryResult);
		if (queryResult.isEmpty()) outIndices[j] = 0; // too far away, whatever...
		else outIndices[j] = (int)tree.getElements()[queryResult[0].getIndex()].id;
	}

	return true;
}

bool SampleUtil::findNearestNeighbors(SKDTree &tree, Eigen::Matrix3Xd &inPoints, set<int> &outIndices, double maxDist) {

	outIndices.clear();
	for (int j = 0; j < inPoints.cols(); j++) {
		SKDT::NamedPoint queryPoint((float)inPoints(0, j), (float)inPoints(1, j), (float)inPoints(2, j));
		Thea::BoundedSortedArray<SKDTree::Neighbor> queryResult(1);
		tree.kClosestElements<Thea::MetricL2>(queryPoint, queryResult, maxDist);
		if (!queryResult.isEmpty()) {
			int pointID = (int)tree.getElements()[queryResult[0].getIndex()].id;
			outIndices.insert(pointID);
		}
	}

	return true;
}

bool SampleUtil::sliceMatrices(Eigen::Matrix3Xd &inMatrix, Eigen::VectorXi &inIndices, Eigen::Matrix3Xd &outMatrix) {

	outMatrix.resize(inMatrix.rows(), inIndices.rows());
	for (int j = 0; j < inIndices.rows(); j++) {
		outMatrix.col(j) = inMatrix.col(inIndices[j]);
	}

	return true;
}

bool SampleUtil::computeNormalPCA(TPointSet &inPoints, vec3 &outCenter, vec3 &outNormal) {

	int numPoints = (int)inPoints.positions.size();
	if (numPoints < 3) {
		cout << "Error: not enough points for PCA" << endl;
		return false;
	}

	// build matrix
	Eigen::MatrixXf matM(numPoints, 3);
	for (int row = 0; row < numPoints; row++) {
		for (int col = 0; col < 3; col++) {
			matM(row, col) = inPoints.positions[row][col];
		}
	}

	// subtract mean
	Eigen::MatrixXf matO = matM;
	Eigen::RowVectorXf vecO = matM.colwise().mean();
	matO.rowwise() -= vecO;
	outCenter = vec3(vecO(0), vecO(1), vecO(2));

	// PCA (perform SVD rather than eigen solver for better numerical precision)
	Eigen::JacobiSVD< Eigen::MatrixXf > svd(matO, Eigen::ComputeThinV);
	Eigen::MatrixXf matV = svd.matrixV();
	outNormal = cml::normalize(vec3(matV(0, 2), matV(1, 2), matV(2, 2)));

	return true;
}

bool SampleUtil::computeAABB(TPointSet &inPoints, vec3 &outBBMin, vec3 &outBBMax) {

	if (inPoints.amount <= 0) {
		cout << "Error: null point set" << endl;
		return false;
	}

	outBBMin = inPoints.positions[0];
	outBBMax = outBBMin;

	for (vec3 &v : inPoints.positions) {
		outBBMin.minimize(v);
		outBBMax.maximize(v);
	}

	return true;
}

bool SampleUtil::computeVolume(TSampleSet &inPoints, double &outVolume) {

	Eigen::Matrix3Xd matP, matN;
	if (!buildMatrices(inPoints, matP, matN)) return false;

	if (matP.cols() < 3) {
		cout << "Error: too few points; cannot estimate OBB" << endl;
		return false;
	}

	// PCA (perform SVD rather than eigen solver for better numerical precision)
	Eigen::Matrix3Xd matO = matP.colwise() - matP.rowwise().mean(); // zero-mean points
	Eigen::JacobiSVD< Eigen::Matrix3Xd > svd(matO, Eigen::ComputeThinU);
	Eigen::Matrix3d matU = svd.matrixU();

	// compute OBB extent
	Eigen::Matrix3Xd matT = matU.transpose() * matO; // transform points to OBB local CS
	Eigen::Vector3d vecOBBMax = matT.rowwise().maxCoeff();
	Eigen::Vector3d vecOBBMin = matT.rowwise().minCoeff();
	Eigen::Vector3d vecOBBExtent = vecOBBMax - vecOBBMin;
	vecOBBExtent = vecOBBExtent.cwiseMax(inPoints.radius);
	double volumeOBB = vecOBBExtent[0] * vecOBBExtent[1] * vecOBBExtent[2];

	// compute AABB extent
	Eigen::Vector3d vecAABBMax = matP.rowwise().maxCoeff();
	Eigen::Vector3d vecAABBMin = matP.rowwise().minCoeff();
	Eigen::Vector3d vecAABBExtent = vecAABBMax - vecAABBMin;
	vecAABBExtent = vecAABBExtent.cwiseMax(inPoints.radius);
	double volumeAABB = vecAABBExtent[0] * vecAABBExtent[1] * vecAABBExtent[2];

	// use smaller one
	outVolume = min(volumeOBB, volumeAABB);

	return true;
}

bool SampleUtil::recomputeNormals(TSampleSet &samples) {

	SKDTree tree;
	SKDTreeData treeData;
	if (!buildKdTree(samples.positions, tree, treeData)) return false;

	double queryRange = samples.radius * 2.5; // UNDONE: param

#pragma omp parallel for
	for (int sampleID = 0; sampleID < samples.amount; sampleID++) {

		vec3 position = samples.positions[sampleID];
		vec3 normal = samples.normals[sampleID];

		// get neighbors
		TPointSet neighbors;
		SKDT::NamedPoint queryPoint(position[0], position[1], position[2]);
		Thea::BoundedSortedArray<SKDTree::Neighbor> queryResult(20);
		tree.kClosestElements<Thea::MetricL2>(queryPoint, queryResult, queryRange);

		// check compatibility of neighbors
		neighbors.positions.clear();
		for (int j = 0; j < queryResult.size(); j++) {
			int neighborID = (int)tree.getElements()[queryResult[j].getIndex()].id;
			vec3 nbPosition = samples.positions[neighborID];
			vec3 nbNormal = samples.normals[neighborID];
			if (cml::dot(normal, nbNormal) >= 0) {
				neighbors.positions.push_back(nbPosition - position);
			}
		}
		int numNeighbors = (int)neighbors.positions.size();
		neighbors.amount = numNeighbors;
		if (numNeighbors < 6) continue; // unstable

		vec3 newCenter, newNormal;
		computeNormalPCA(neighbors, newCenter, newNormal);
		if (cml::dot(normal, newNormal) < 0) newNormal = -newNormal;

		samples.normals[sampleID] = newNormal;
	}

	return true;
}

bool SampleUtil::buildMatrices(vector<vec3> &inPoints, Eigen::Matrix3Xd &outMat) {

	outMat.resize(3, (int)inPoints.size());

#pragma omp parallel for
	for (int j = 0; j < (int)inPoints.size(); j++) {
		outMat.col(j) = Eigen::Vector3d(vec3d(inPoints[j]).data());
	}

	return true;
}

bool SampleUtil::buildMatrices(TPointSet &inPoints, Eigen::Matrix3Xd &outMat) {

	outMat.resize(3, inPoints.amount);

#pragma omp parallel for
	for (int j = 0; j < inPoints.amount; j++) {
		outMat.col(j) = Eigen::Vector3d(vec3d(inPoints.positions[j]).data());
	}

	return true;
}

bool SampleUtil::buildMatrices(TPointSet &inPoints, Eigen::Matrix3Xd &outMatP, Eigen::Matrix3Xd &outMatN) {	

	outMatP.resize(3, inPoints.amount);
	outMatN.resize(3, inPoints.amount);

#pragma omp parallel for
	for (int j = 0; j < inPoints.amount; j++) {
		outMatP.col(j) = Eigen::Vector3d(vec3d(inPoints.positions[j]).data());
		outMatN.col(j) = Eigen::Vector3d(vec3d(inPoints.normals[j]).data());
	}

	return true;
}

bool SampleUtil::subSampleMatrices(Eigen::Matrix3Xd &inMat, Eigen::VectorXi &outIdx, int subsample) {

	int cols = (int)inMat.cols();
	if (cols <= subsample) {
		outIdx.resize(cols);
#pragma omp parallel for
		for (int j = 0; j < cols; j++) {
			outIdx[j] = j;
		}
	} else {
		outIdx.resize(subsample);
		vector<int> randIdx(cols);
		for (int j = 0; j < cols; j++) randIdx[j] = j;
		random_shuffle(randIdx.begin(), randIdx.end());
#pragma omp parallel for
		for (int j = 0; j < subsample; j++) {
			outIdx[j] = randIdx[j];
		}
	}

	return true;
}

bool SampleUtil::saveSample(string fileName, TSampleSet &samples) {

	vector<TPlySample> outSamples(samples.amount);
	for (int sampleID = 0; sampleID < samples.amount; sampleID++) {
		TPlySample &sample = outSamples[sampleID];
		sample.p = samples.positions[sampleID];
		sample.n = samples.normals[sampleID];
		sample.f = samples.indices[sampleID];
	}

	vector<string> outProperties;
	outProperties.push_back("float x");
	outProperties.push_back("float y");
	outProperties.push_back("float z");
	outProperties.push_back("float nx");
	outProperties.push_back("float ny");
	outProperties.push_back("float nz");
	outProperties.push_back("uint faceID");

	stringstream ss; ss << "SAMPLE_RADIUS " << samples.radius;
	vector<string> outComments;
	outComments.push_back(ss.str());

	if (!PlyExporter::exportSample(fileName, &outSamples, &outProperties, &outComments)) return false;

	return true;
}

bool SampleUtil::loadSample(string fileName, TSampleSet &samples) {

	vector<TPlySample> inSamples;
	vector<string> comments;

	if (!PlyLoader::loadSample(fileName, &inSamples, &comments)) return false;
	samples.positions.clear();
	samples.normals.clear();
	samples.indices.clear();
	for (TPlySample &sample : inSamples) {
		samples.positions.push_back(sample.p);
		samples.normals.push_back(sample.n);
		samples.indices.push_back(sample.f);
	}
	samples.amount = (int)samples.positions.size();

	samples.radius = 0;
	for (auto line : comments) {
		string::size_type pos = line.find("SAMPLE_RADIUS");
		if (pos != string::npos) {
			stringstream ss(line.substr(pos));
			string s; ss >> s;
			ss >> samples.radius;
			break;
		}
	}
	if (samples.radius == 0) {
		cout << "Error: no sample radius information" << endl;
		return false;
	}

	return true;
}

double SampleUtil::getPreciseRandomNumber() {

	long long rng = rand()*(RAND_MAX + 1) + rand();
	return double(rng) / ((RAND_MAX + 1)*(RAND_MAX + 1));
}

int SampleUtil::samplingCDF(vector<double> &cdf) {

	int n = (int)cdf.size();

	double rng = SampleUtil::getPreciseRandomNumber();
	// binary search
	int p1 = -1;
	int p2 = n - 1;
	while (p1 + 1 < p2) {
		int pm = (p1 + p2) / 2;
		if (cdf[pm] < rng) {
			p1 = pm;
		} else {
			p2 = pm;
		}
	}

	int sample = p2;
	if (rng > 0) {
		// choose first valid sample (leftmost among samples with the same CDF)
		while (sample > 0 && cdf[sample] == cdf[sample-1]) sample--;
	} else {
		// choose first valid sample (first with non-zero CDF)
		while (sample < n-1 && cdf[sample] == 0) sample++;
	}

	return sample;
}