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

#include "FeatureMeshCurvature.h"

#include <unordered_set>

#include "Mesh/MeshUtil.h"

#include <Eigen/Eigen>

using namespace StyleSynthesis;

bool FeatureMeshCurvature::computeCurvature(TTriangleMesh &mesh, vector<TCurvature> &curvature) {

	int numVertices = mesh.amount;
	int numFaces = (int)mesh.indices.size();

	vector<vec3d> curvatureTensor(0);
	curvatureTensor.resize(numVertices, vec3d(0.0, 0.0, 0.0));

	vector<double> pointArea(numVertices, 0.0);

	vector<vec3d> vertexDir1(numVertices), vertexDir2(numVertices);
	vector<bool> vertexFlag(numVertices, false);
	for (int faceID = 0; faceID < numFaces; faceID++) {
		vec3i faceIdx = mesh.indices[faceID];
		for (int j = 0; j < 3; j++) {
			vec3d dir = vec3d(mesh.positions[faceIdx[(j + 1) % 3]] - mesh.positions[faceIdx[j]]); // any direction
			vec3d n = vec3d(mesh.normals[faceIdx[j]]);
			if (cml::cross(dir, n).length_squared()) {
				vertexDir1[faceIdx[j]] = dir;
				vertexFlag[faceIdx[j]] = true;
			}
		}
	}

#pragma omp parallel for
	for (int vertID = 0; vertID < numVertices; vertID++) {
		if (vertexFlag[vertID]) {
			vec3d vertexN = vec3d(mesh.normals[vertID]);
			vertexDir1[vertID] = cml::normalize(cml::cross(vertexDir1[vertID], vertexN));
			vertexDir2[vertID] = cml::normalize(cml::cross(vertexN, vertexDir1[vertID]));
		}
	}

#pragma omp parallel for
	for (int faceID = 0; faceID < numFaces; faceID++) {

		vec3i faceIdx = mesh.indices[faceID];

		vec3d vertexP[3], vertexN[3];
		for (int j = 0; j < 3; j++) {
			vertexP[j] = vec3d(mesh.positions[faceIdx[j]]);
			vertexN[j] = cml::normalize(vec3d(mesh.normals[faceIdx[j]]));
		}
		vec3d faceEdge[] = { vertexP[2] - vertexP[1], vertexP[0] - vertexP[2], vertexP[1] - vertexP[0] };

		vec3d cornerArea;
		if (!computeCornerArea(faceEdge, cornerArea)) continue; // degenerated face

		// N-T-B coordinate system
		vec3d faceCS[3];
		faceCS[1] = cml::normalize(faceEdge[0]);
		faceCS[0] = cml::normalize(cml::cross(faceCS[1], cml::normalize(faceEdge[1])));
		faceCS[2] = cml::normalize(cml::cross(faceCS[0], faceCS[1]));

		// build linear system
		Eigen::MatrixXd matA = Eigen::MatrixXd::Zero(6, 3);
		Eigen::VectorXd matB = Eigen::MatrixXd::Zero(6, 1);
		for (int j = 0; j < 3; j++) {
			double u = cml::dot(faceEdge[j], faceCS[1]);
			double v = cml::dot(faceEdge[j], faceCS[2]);
			double nu = cml::dot(vertexN[(j + 2) % 3] - vertexN[(j + 1) % 3], faceCS[1]);
			double nv = cml::dot(vertexN[(j + 2) % 3] - vertexN[(j + 1) % 3], faceCS[2]);

			matA.row(j * 2 + 0) << u, v, 0;
			matA.row(j * 2 + 1) << 0, u, v;

			matB.row(j * 2 + 0) << nu;
			matB.row(j * 2 + 1) << nv;
		}

		// solve linear system: A * X = B
		if (!matA.allFinite() || !matB.allFinite()) {
			cout << "Error: matrix is invalid!" << endl;
		}
		Eigen::VectorXd matX = matA.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(matB);
		vec3d faceTensor = vec3d(matX(0), matX(1), matX(2));

		// project back to vertex CS
		vec3d vertexTensor[3];
		for (int j = 0; j < 3; j++) {
			vec3d vertexCS[] = { vertexN[j], vertexDir1[faceIdx[j]], vertexDir2[faceIdx[j]] };
			if (!projectCurvatureTensor(faceCS, faceTensor, vertexCS, vertexTensor[j])) {
				cout << "Error: projectCurvatureTensor" << endl;
			}
		}

		// add weighted tensor
#pragma omp critical
		{
			for (int j = 0; j < 3; j++) {
				curvatureTensor[faceIdx[j]] += vertexTensor[j] * cornerArea[j];
				pointArea[faceIdx[j]] += cornerArea[j];
			}
		}
	}

	curvature.clear();
	curvature.resize(numVertices);

	// extract principal curvature and direction from tensor
#pragma omp parallel for
	for (int vertexID = 0; vertexID < numVertices; vertexID++) {

		if (pointArea[vertexID] > 0) {
			curvatureTensor[vertexID] /= pointArea[vertexID]; // normalize

			// curvature tensor
			vec3d tensor = curvatureTensor[vertexID];
			Eigen::Matrix2d matK;
			matK << tensor[0], tensor[1], tensor[1], tensor[2];

			// solve for eigen values and eigen vectors
			Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigenSolver(matK);
			if (eigenSolver.info() != Eigen::Success) {
				cout << "Error: fail to solve for principal curvatures" << endl;
			}

			vec3d localX = vertexDir1[vertexID];
			vec3d localY = vertexDir2[vertexID];
			vec3d localZ = cml::normalize(vec3d(mesh.normals[vertexID]));

			TCurvature &curv = curvature[vertexID];
			curv.k1 = (float)eigenSolver.eigenvalues()(1);
			curv.k2 = (float)eigenSolver.eigenvalues()(0);
			curv.d1 = vec3(cml::normalize(localX * eigenSolver.eigenvectors()(0, 1) + localY * eigenSolver.eigenvectors()(1, 1)));
			curv.d2 = vec3(cml::normalize(localX * eigenSolver.eigenvectors()(0, 0) + localY * eigenSolver.eigenvectors()(1, 0)));

			// re-orient curvature direction to make {n, d1, d2} a right-hand-side coordinate system
			if (cml::dot(cml::cross(vec3d(curv.d1), vec3d(curv.d2)), localZ) < 0) curv.d2 = -curv.d2;
		} else {
			// zero weight - point Voronoi area underflows...

			TCurvature &curv = curvature[vertexID];
			vec3 n = mesh.normals[vertexID];
			curv.d1 = cml::normalize(cml::cross(n, n[0] < 0.5f ? vec3(1.0f, 0.0f, 0.0f) : vec3(0.0f, 1.0f, 0.0f)));
			curv.d2 = cml::normalize(cml::cross(n, curv.d1));
			curv.k1 = 0;
			curv.k2 = 0;
		}
	}

	return true;
}

bool FeatureMeshCurvature::computeCurvatureDerivative(TTriangleMesh &mesh, vector<TCurvature> &curvature, vector<vec4> &derivative) {

	int numVertices = mesh.amount;
	int numFaces = (int)mesh.indices.size();

	vector<double> totalWeights(numVertices, 0.0);
	vector<vec4d> accumDerivative(numVertices, vec4d(0.0, 0.0, 0.0, 0.0));

#pragma omp parallel for
	for (int faceID = 0; faceID < numFaces; faceID++) {

		vec3i faceIdx = mesh.indices[faceID];

		vec3d vertexP[3], vertexN[3];
		for (int j = 0; j < 3; j++) {
			vertexP[j] = vec3d(mesh.positions[faceIdx[j]]);
			vertexN[j] = cml::normalize(vec3d(mesh.normals[faceIdx[j]]));
		}
		vec3d faceEdge[] = { vertexP[2] - vertexP[1], vertexP[0] - vertexP[2], vertexP[1] - vertexP[0] };

		vec3d cornerArea;
		if (!computeCornerArea(faceEdge, cornerArea)) continue; // degenerated face

		// curvature tensor
		TCurvature vertexCurvatures[3];
		for (int j = 0; j < 3; j++) {
			vertexCurvatures[j] = curvature[faceIdx[j]];
		}

		// N-T-B coordinate system
		vec3d faceCS[3];
		faceCS[1] = cml::normalize(faceEdge[0]);
		faceCS[0] = cml::normalize(cml::cross(faceCS[1], cml::normalize(faceEdge[1])));
		faceCS[2] = cml::normalize(cml::cross(faceCS[0], faceCS[1]));

		// project tensor to face's CS
		vec3d faceTensor[3];
		for (int j = 0; j < 3; j++) {
			vec3d vertexCS[] = { vertexN[j], vertexCurvatures[j].d1, vertexCurvatures[j].d2 };
			vec3d vertexTensor = vec3(vertexCurvatures[j].k1, 0.0f, vertexCurvatures[j].k2);
			if (!projectCurvatureTensor(vertexCS, vertexTensor, faceCS, faceTensor[j])) {
				cout << "Error: projectCurvatureTensor" << endl;
			}
		}

		// build linear system
		Eigen::MatrixXd matA = Eigen::MatrixXd::Zero(12, 4);
		Eigen::VectorXd matB = Eigen::MatrixXd::Zero(12, 1);
		for (int j = 0; j < 3; j++) {
			double u = cml::dot(faceEdge[j], faceCS[1]);
			double v = cml::dot(faceEdge[j], faceCS[2]);
			vec3d d = faceTensor[(j + 2) % 3] - faceTensor[(j + 1) % 3];

			matA.row(j * 4 + 0) << u, v, 0, 0;
			matA.row(j * 4 + 1) << 0, u, v, 0;
			matA.row(j * 4 + 2) << 0, u, v, 0;
			matA.row(j * 4 + 3) << 0, 0, u, v;

			matB.row(j * 4 + 0) << d[0];
			matB.row(j * 4 + 1) << d[1];
			matB.row(j * 4 + 2) << d[1];
			matB.row(j * 4 + 3) << d[2];
		}

		// solve linear system: A * X = B
		if (!matA.allFinite() || !matB.allFinite()) {
			cout << "Error: matrix is invalid!" << endl;
		}
		Eigen::VectorXd matX = matA.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(matB);
		vec4d faceDerv = vec4d(matX(0), matX(1), matX(2), matX(3));

		// project back to vertex CS
		vec4d vertexDerv[3];
		for (int j = 0; j < 3; j++) {
			vec3d vertexCS[] = { vertexN[j], vertexCurvatures[j].d1, vertexCurvatures[j].d2 };
			if (!projectCurvatureDerivative(faceCS, faceDerv, vertexCS, vertexDerv[j])) {
				cout << "Error: projectCurvatureDerivative" << endl;
			}
		}

		// add weighted derivation
#pragma omp critical
		{
			for (int j = 0; j < 3; j++) {
				accumDerivative[faceIdx[j]] += vertexDerv[j] * cornerArea[j];
				totalWeights[faceIdx[j]] += cornerArea[j];
			}
		}
	}

	// output derivative

	derivative.resize(numVertices);
#pragma omp parallel for
	for (int vertexID = 0; vertexID < numVertices; vertexID++) {
		if (totalWeights[vertexID] > 0) {
			derivative[vertexID] = vec4(accumDerivative[vertexID] / totalWeights[vertexID]);
		} else {
			derivative[vertexID] = vec4(0.0f, 0.0f, 0.0f, 0.0f);
		}
	}

	return true;
}

bool FeatureMeshCurvature::smoothNormal(TTriangleMesh &mesh, vector<TCurvature> &curvature) {

	int numVertices = mesh.amount;

	vector<vector<int>> vertexGraph;
	if (!MeshUtil::buildVertexGraph(mesh, vertexGraph)) return false;

	vector<float> pointArea;
	if (!computePointArea(mesh, pointArea)) return false;

	float featureSize;
	if (!computeFeatureSize(mesh, curvature, featureSize)) return false;
	float invSigmaSq = 1 / cml::sqr(0.5f * featureSize);

	vector<vec3> newNormal = mesh.normals;

#pragma omp parallel for
	for (int vertexID = 0; vertexID < numVertices; vertexID++) {

		vec3 vertexP = mesh.positions[vertexID];
		vec3 vertexN = mesh.normals[vertexID];

		// find neighbors and compute weights

		vector<int> neighborVertices;
		vector<double> neighborWeights;

		unordered_set<int> visitedVertices;
		neighborVertices.push_back(vertexID);
		neighborWeights.push_back((double)pointArea[vertexID]);

		int head = 0;
		while (head < (int)neighborVertices.size()) {
			int currentID = neighborVertices[head];
			for (int neighborID : vertexGraph[vertexID]) {
				if (visitedVertices.find(neighborID) != visitedVertices.end()) continue;
				visitedVertices.insert(neighborID);
				// compute weight
				float dist2 = (mesh.positions[neighborID] - vertexP).length_squared();
				if (dist2 * invSigmaSq >= 9) continue; // 3 standard deviation
				double weight = exp(-0.5 * dist2 * invSigmaSq);
				float cosN = cml::dot(mesh.normals[neighborID], vertexN);
				if (cosN <= 0) continue; // incompatible normal
				weight *= cosN;
				weight *= pointArea[neighborID];
				neighborVertices.push_back(neighborID);
				neighborWeights.push_back(weight);
			}
			head++;
		}

		// diffuse normal

		vec3d accumWeightedNormal(0.0, 0.0, 0.0);

		int numNeighbors = (int)neighborVertices.size();
		for (int id = 0; id < numNeighbors; id++) {
			int neighborID = neighborVertices[id];
			double neighborWeight = neighborWeights[id];
			vec3d neighborNormal = vec3d(mesh.normals[neighborID]);

			accumWeightedNormal += neighborNormal * neighborWeight;
		}
		if (accumWeightedNormal.length_squared()) accumWeightedNormal.normalize();

		newNormal[vertexID] = accumWeightedNormal;
	}

	mesh.normals.swap(newNormal);

	return true;
}

bool FeatureMeshCurvature::smoothCurvature(TTriangleMesh &mesh, vector<TCurvature> &curvature) {

	int numVertices = mesh.amount;

	vector<vector<int>> vertexGraph;
	if (!MeshUtil::buildVertexGraph(mesh, vertexGraph)) return false;

	vector<float> pointArea;
	if (!computePointArea(mesh, pointArea)) return false;

	float featureSize;
	if (!computeFeatureSize(mesh, curvature, featureSize)) return false;
	float invSigmaSq = 1 / cml::sqr(0.5f * featureSize);;

	vector<TCurvature> newCurvature = curvature;

#pragma omp parallel for
	for (int vertexID = 0; vertexID < numVertices; vertexID++) {

		vec3 vertexP = mesh.positions[vertexID];
		vec3 vertexN = mesh.normals[vertexID];
		TCurvature &vertexCurvature = curvature[vertexID];

		// find neighbors and compute weights

		vector<int> neighborVertices;
		vector<double> neighborWeights;

		unordered_set<int> visitedVertices;
		neighborVertices.push_back(vertexID);
		neighborWeights.push_back((double)pointArea[vertexID]);

		int head = 0;
		while (head < (int)neighborVertices.size()) {
			int currentID = neighborVertices[head];
			for (int neighborID : vertexGraph[vertexID]) {
				if (visitedVertices.find(neighborID) != visitedVertices.end()) continue;
				visitedVertices.insert(neighborID);
				// compute weight
				float dist2 = (mesh.positions[neighborID] - vertexP).length_squared();
				if (dist2 * invSigmaSq >= 9) continue; // 3 standard deviation
				double weight = exp(-0.5 * dist2 * invSigmaSq);
				float cosN = cml::dot(mesh.normals[neighborID], vertexN);
				if (cosN <= 0) continue; // incompatible normal
				weight *= cosN;
				weight *= pointArea[neighborID];
				neighborVertices.push_back(neighborID);
				neighborWeights.push_back(weight);
			}
			head++;
		}

		// diffuse tensor

		vec3d vertexCS[] = { vertexN, vertexCurvature.d1, vertexCurvature.d2 };
		vec3d accumWeightedTensor(0.0, 0.0, 0.0);
		double accumWeights = 0;

		int numNeighbors = (int)neighborVertices.size();
		for (int id = 0; id < numNeighbors; id++) {
			int neighborID = neighborVertices[id];
			double neighborWeight = neighborWeights[id];

			TCurvature &neighborCurvature = curvature[neighborID];
			vec3d neighborCS[] = { mesh.normals[neighborID], neighborCurvature.d1, neighborCurvature.d2 };
			vec3d neighborTensor((double)neighborCurvature.k1, 0.0, (double)neighborCurvature.k2);
			vec3d projTensor;
			if (!projectCurvatureTensor(neighborCS, neighborTensor, vertexCS, projTensor)) {
				cout << "Error: projCurvature tensor" << endl;
			}
			
			accumWeightedTensor += projTensor * neighborWeight;
			accumWeights += neighborWeight;
		}

		vec3d newTensor = accumWeightedTensor / accumWeights;

		// compute curvature

		Eigen::Matrix2d matK;
		matK << newTensor[0], newTensor[1], newTensor[1], newTensor[2];
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigenSolver(matK);
		if (eigenSolver.info() != Eigen::Success) {
			cout << "Error: fail to solve for principal curvatures" << endl;
		}
		auto &curv = newCurvature[vertexID];
		curv.k1 = (float)eigenSolver.eigenvalues()(1);
		curv.k2 = (float)eigenSolver.eigenvalues()(0);
	}

	curvature.swap(newCurvature);

	return true;
}

bool FeatureMeshCurvature::smoothCurvatureDerivative(TTriangleMesh &mesh, vector<TCurvature> &curvature, vector<vec4> &derivative) {

	int numVertices = mesh.amount;

	vector<vector<int>> vertexGraph;
	if (!MeshUtil::buildVertexGraph(mesh, vertexGraph)) return false;

	vector<float> pointArea;
	if (!computePointArea(mesh, pointArea)) return false;

	float featureSize;
	if (!computeFeatureSize(mesh, curvature, featureSize)) return false;
	float invSigmaSq = 1 / cml::sqr(0.5f * featureSize);

	vector<vec4> newDerivative = derivative;

#pragma omp parallel for
	for (int vertexID = 0; vertexID < numVertices; vertexID++) {

		vec3 vertexP = mesh.positions[vertexID];
		vec3 vertexN = mesh.normals[vertexID];
		TCurvature &vertexCurvature = curvature[vertexID];;

		// find neighbors and compute weights

		vector<int> neighborVertices;
		vector<double> neighborWeights;

		unordered_set<int> visitedVertices;
		neighborVertices.push_back(vertexID);
		neighborWeights.push_back((double)pointArea[vertexID]);

		int head = 0;
		while (head < (int)neighborVertices.size()) {
			int currentID = neighborVertices[head];
			for (int neighborID : vertexGraph[vertexID]) {
				if (visitedVertices.find(neighborID) != visitedVertices.end()) continue;
				visitedVertices.insert(neighborID);
				// compute weight
				float dist2 = (mesh.positions[neighborID] - vertexP).length_squared();
				if (dist2 * invSigmaSq >= 9) continue; // 3 standard deviation
				double weight = exp(-0.5 * dist2 * invSigmaSq);
				float cosN = cml::dot(mesh.normals[neighborID], vertexN);
				if (cosN <= 0) continue; // incompatible normal
				weight *= cosN;
				weight *= pointArea[neighborID];
				neighborVertices.push_back(neighborID);
				neighborWeights.push_back(weight);
			}
			head++;
		}

		// diffuse derivative

		vec3d vertexCS[] = { vertexN, vertexCurvature.d1, vertexCurvature.d2 };
		vec4d accumWeightedDerivative(0.0, 0.0, 0.0, 0.0);
		double accumWeights = 0;

		int numNeighbors = (int)neighborVertices.size();
		for (int id = 0; id < numNeighbors; id++) {
			int neighborID = neighborVertices[id];
			double neighborWeight = neighborWeights[id];

			TCurvature &neighborCurvature = curvature[neighborID];
			vec4d neighborDerivative = derivative[neighborID];
			vec3d neighborCS[] = { mesh.normals[neighborID], neighborCurvature.d1, neighborCurvature.d2 };
			vec4d projDerivative;
			if (!projectCurvatureDerivative(neighborCS, neighborDerivative, vertexCS, projDerivative)) {
				cout << "Error: projCurvatureDerivative" << endl;
			}

			accumWeightedDerivative += projDerivative * neighborWeight;
			accumWeights += neighborWeight;
		}

		newDerivative[vertexID] = accumWeightedDerivative / accumWeights;
	}

	derivative.swap(newDerivative);

	return true;
}

bool FeatureMeshCurvature::computeFeatureSize(TTriangleMesh &mesh, vector<TCurvature> &curvature, float &featureSize) {

	// ref: rtsc.cc

	// feature size = 1% of the reciprocal of the 10-th percentile curvature

	vec3 bbMin, bbMax;
	if (!MeshUtil::computeAABB(mesh, bbMin, bbMax)) return false;
	float maxFeatureSize = 0.025f * (bbMax - bbMin).length();

	int numVertices = mesh.amount;
	vector<float> samples(numVertices * 2);
	for (int vertID = 0; vertID < numVertices; vertID++) {
		samples[vertID * 2 + 0] = fabs(curvature[vertID].k1);
		samples[vertID * 2 + 1] = fabs(curvature[vertID].k2);
	}

	int n = (int)(samples.size() * 0.1f);
	nth_element(samples.begin(), samples.begin() + n, samples.end());
	float prctl = samples[n];

	featureSize = prctl ? min(0.01f / prctl, maxFeatureSize) : maxFeatureSize;

	return true;
}

bool FeatureMeshCurvature::computePointArea(TTriangleMesh &mesh, vector<float> &area) {

	int numVertices = mesh.amount;
	int numFaces = (int)mesh.indices.size();

	vector<vec3d> faceCornerArea(numFaces);
#pragma omp parallel for
	for (int faceID = 0; faceID < numFaces; faceID++) {
		vec3i faceIdx = mesh.indices[faceID];
		vec3d faceEdge[3];
		for (int k = 0; k < 3; k++) {
			faceEdge[k] = vec3d(mesh.positions[faceIdx[(k + 2) % 3]] - mesh.positions[faceIdx[(k + 1) % 3]]);
		}
		computeCornerArea(faceEdge, faceCornerArea[faceID]);
	}

	area.assign(numVertices, 0.0f);
	for (int faceID = 0; faceID < numFaces; faceID++) {
		vec3i faceIdx = mesh.indices[faceID];
		vec3d cornerArea = faceCornerArea[faceID];
		for (int k = 0; k < 3; k++) {
			area[faceIdx[k]] += (float)cornerArea[k];
		}
	}

	return true;
}

bool FeatureMeshCurvature::computeNormalCurvature(TCurvature &tensor, vec3 &direction, float &curvature) {

	float dk1 = cml::sqr(cml::dot(cml::normalize(tensor.d1), direction));
	float dk2 = cml::sqr(cml::dot(cml::normalize(tensor.d2), direction));

	if (dk1 == 0 && dk2 == 0) {
		curvature = tensor.k1; // should be +oo? whatever...
	} else {
		curvature = (tensor.k1 * dk1 + tensor.k2 * dk2) / (dk1 + dk2);
	}

	return true;
}

bool FeatureMeshCurvature::projectCurvatureTensor(
	const vec3d(&oldCS)[3], vec3d oldTensor,
	const vec3d(&newCS)[3], vec3d &newTensor)
{
	// ref: TriMesh_curvature.cc

	// rotate new CS to be perp to normal of old CS
	double normalCos = cml::dot(oldCS[0], newCS[0]);
	vec3d newU, newV;
	if (normalCos > -1.0) {
		vec3d perp = oldCS[0] - newCS[0] * normalCos;
		vec3d dperp = (oldCS[0] + newCS[0]) / (normalCos + 1);
		newU = newCS[1] - dperp * cml::dot(newCS[1], perp);
		newV = newCS[2] - dperp * cml::dot(newCS[2], perp);
	} else {
		newU = -newCS[1];
		newV = -newCS[2];
	}

	// reproject curvature tensor
	double u1 = cml::dot(newU, oldCS[1]);
	double v1 = cml::dot(newU, oldCS[2]);
	double u2 = cml::dot(newV, oldCS[1]);
	double v2 = cml::dot(newV, oldCS[2]);
	newTensor[0] = oldTensor[0] * u1*u1 + oldTensor[1] * (2.0*u1*v1) + oldTensor[2] * v1*v1;
	newTensor[1] = oldTensor[0] * u1*u2 + oldTensor[1] * (u1*v2 + u2*v1) + oldTensor[2] * v1*v2;
	newTensor[2] = oldTensor[0] * u2*u2 + oldTensor[1] * (2.0*u2*v2) + oldTensor[2] * v2*v2;

	if (newTensor != newTensor) {
#pragma omp critical
		{
			cout << "NaN in project tensor: " << endl;
			/*
			cout << newTensor << endl << oldTensor << endl;
			cout << "new CS: " << endl;
			for (int k = 0; k < 3; k++) cout << newCS[k] << endl;
			cout << "old CS: " << endl;
			for (int k = 0; k < 3; k++) cout << oldCS[k] << endl;
			*/
		}
	}

	return true;
}

bool FeatureMeshCurvature::projectCurvatureDerivative(
	const vec3d(&oldCS)[3], vec4d oldDerv,
	const vec3d(&newCS)[3], vec4d &newDerv)
{
	// ref: TriMesh_curvature.cc

	// rotate new CS to be perp to normal of old CS
	double normalCos = cml::dot(oldCS[0], newCS[0]);
	vec3d newU, newV;
	if (normalCos > -1.0) {
		vec3d perp = oldCS[0] - newCS[0] * normalCos;
		vec3d dperp = (oldCS[0] + newCS[0]) / (normalCos + 1);
		newU = newCS[1] - dperp * cml::dot(newCS[1], perp);
		newV = newCS[2] - dperp * cml::dot(newCS[2], perp);
	} else {
		newU = -newCS[1];
		newV = -newCS[2];
	}

	// reproject derivative tensor
	double u1 = cml::dot(newU, oldCS[1]);
	double v1 = cml::dot(newU, oldCS[2]);
	double u2 = cml::dot(newV, oldCS[1]);
	double v2 = cml::dot(newV, oldCS[2]);
	newDerv[0] = oldDerv[0] * u1*u1*u1
		+ oldDerv[1] * 3.0*u1*u1*v1
		+ oldDerv[2] * 3.0*u1*v1*v1
		+ oldDerv[3] * v1*v1*v1;
	newDerv[1] = oldDerv[0] * u1*u1*u2
		+ oldDerv[1] * (u1*u1*v2 + 2.0*u2*u1*v1)
		+ oldDerv[2] * (u2*v1*v1 + 2.0*u1*v1*v2)
		+ oldDerv[3] * v1*v1*v2;
	newDerv[2] = oldDerv[0] * u1*u2*u2
		+ oldDerv[1] * (u2*u2*v1 + 2.0*u1*u2*v2)
		+ oldDerv[2] * (u1*v2*v2 + 2.0*u2*v2*v1)
		+ oldDerv[3] * v1*v2*v2;
	newDerv[3] = oldDerv[0] * u2*u2*u2
		+ oldDerv[1] * 3.0*u2*u2*v2
		+ oldDerv[2] * 3.0*u2*v2*v2
		+ oldDerv[3] * v2*v2*v2;

	return true;
}

bool FeatureMeshCurvature::computeCornerArea(const vec3d(&faceEdge)[3], vec3d &cornerArea) {

	// ref: TriMesh_pointareas.cc

	double faceArea = cml::cross(faceEdge[0], faceEdge[1]).length() / 2;
	if (faceArea <= 0) return false; // degenerated face

	double lengthSq[3];
	for (int k = 0; k < 3; k++) {
		lengthSq[k] = faceEdge[k].length_squared();
	}
	double edgeWeight[3];
	for (int k = 0; k < 3; k++) {
		edgeWeight[k] = lengthSq[k] * (lengthSq[(k + 1) % 3] + lengthSq[(k + 2) % 3] - lengthSq[k]);
		if (edgeWeight[k] <= 0) {
			cornerArea[(k + 1) % 3] = -0.25 * lengthSq[(k + 2) % 3] * faceArea / cml::dot(faceEdge[k], faceEdge[(k + 2) % 3]);
			cornerArea[(k + 2) % 3] = -0.25 * lengthSq[(k + 1) % 3] * faceArea / cml::dot(faceEdge[k], faceEdge[(k + 1) % 3]);
			cornerArea[k] = faceArea - cornerArea[(k + 1) % 3] - cornerArea[(k + 2) % 3];
			return true;
		}
	}
	double edgeWeightScale = 0.5 * faceArea / (edgeWeight[0] + edgeWeight[1] + edgeWeight[2]);
	for (int k = 0; k < 3; k++) {
		cornerArea[k] = edgeWeightScale * (edgeWeight[(k + 1) % 3] + edgeWeight[(k + 2) % 3]);
	}

	return true;
}