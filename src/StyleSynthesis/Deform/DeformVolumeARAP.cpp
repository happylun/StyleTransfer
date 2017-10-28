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

#include "DeformVolumeARAP.h"

#include <unordered_set>
#include <set>

#pragma warning(push)
#pragma warning(disable:4018)
#pragma warning(disable:4244)
#include <Eigen/Sparse>
#pragma warning(pop)

#include "Library/libiglHelperBegin.h"
#include <igl/arap.h>
#include "Library/libiglHelperEnd.h"

#include "Mesh/MeshUtil.h"

#include "Utility/PlyExporter.h"

#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

DeformVolumeARAP::DeformVolumeARAP() {
}

DeformVolumeARAP::~DeformVolumeARAP() {
}

bool DeformVolumeARAP::loadMesh(TTetrahedralMesh &mesh) {

	mpMesh = &mesh;

	return true;
}

bool DeformVolumeARAP::loadConstraints(map<int, vec3d> &constraints) {

	mConstraints = constraints;

	return true;
}

bool DeformVolumeARAP::process() {

	if (!computeConstraints()) return false;
	if (!solveDeformation()) return false;

	return true;
}

bool DeformVolumeARAP::computeConstraints() {

	if (!mConstraints.empty()) return true;

	// HACK: manual deformation handle constraints

	cout << "Computing constraints..." << endl;

	mConstraints.clear();

	int numVerts = mpMesh->amount;

	for (int vertID = 0; vertID < numVerts; vertID++) {
		vec3 pos = mpMesh->positions[vertID];

		if (pos[0] > 0.4f) {
			mConstraints[vertID] = vec3d(0.5f, -1.0f, 0.0f);
			//if (pos[1] > 0.1) mConstraints[vertID] = vec3d(0.0f, 0.5f, 0.0f);
			//if (pos[1] <-0.1) mConstraints[vertID] = vec3d(0.0f, -0.5f, 0.0f);
		} else if (pos[0] < -0.4f) {
			mConstraints[vertID] = vec3d(0.0f, 0.0f, 0.0f);
		}
	}

	return true;
}

bool DeformVolumeARAP::solveDeformation() {

	int numVertices = mpMesh->amount;
	int numTets = (int)mpMesh->indices.size();
	int numConstraints = (int)mConstraints.size();

	cout << "Deforming volume..." << endl;

	mDeformedMesh = *mpMesh;
	for (auto &cons : mConstraints) {
		int vertID = cons.first;
		vec3d disp = cons.second;
		mDeformedMesh.positions[vertID] += (vec3)disp;
	}

	if (numConstraints == 0) return true; // no deformation
	if (numConstraints == numVertices) return true; // no unknown deformation

	cout << "Building matrices" << endl;

	Eigen::MatrixXd matV;
	Eigen::MatrixXi matT;
	Eigen::VectorXi vecB(numConstraints);
	Eigen::MatrixXd matBC(numConstraints, 3);

	if (!MeshUtil::tet2mat(*mpMesh, matV, matT)) return false;

	int consID = 0;
	for (auto &cons : mConstraints) {
		int vertID = cons.first;
		vec3d disp = cons.second;
		vec3d newPos = (vec3d)(mpMesh->positions[vertID]) + disp;
		vecB(consID) = vertID;
		matBC.row(consID) = Eigen::Vector3d(newPos[0], newPos[1], newPos[2]);
		consID++;
	}

	cout << "Precomputing ARAP" << endl;

	igl::ARAPData arapData;
	arapData.max_iter = 100;
#pragma omp parallel num_threads(1)
	igl::arap_precomputation(matV, matT, 3, vecB, arapData);

	cout << "Solving ARAP" << endl;

	Eigen::MatrixXd matX = matV; // initial guess and final result
#pragma omp parallel num_threads(1)
	igl::arap_solve(matBC, arapData, matX);

	cout << "Computing displacement" << endl;

	for (int vertID = 0; vertID < numVertices; vertID++) {
		mDeformedMesh.positions[vertID] = vec3d(matX(vertID, 0), matX(vertID, 1), matX(vertID, 2));
	}

	return true;
}

bool DeformVolumeARAP::output(TTetrahedralMesh &mesh) {

	mesh = mDeformedMesh;

	return true;
}

bool DeformVolumeARAP::visualize(string filePath) {

	// TODO:

	/*
	if (true) {

		unordered_set<int> fixedVerts;
		unordered_set<int> movedVerts;
		for (auto &it : mConstraints) {
			if (it.second == vec3d(0.0, 0.0, 0.0)) {
				fixedVerts.insert(it.first);
			} else {
				movedVerts.insert(it.first);
			}
		}
		vector<vec3i> faceColors(mpMesh->indices.size(), vec3i(127, 127, 127));
		for (int faceID = 0; faceID < (int)mpMesh->indices.size(); faceID++) {
			vec3i idx = mpMesh->indices[faceID];
			bool fixed = true;
			bool moved = true;
			for (int k = 0; k < 3; k++) {
				if (fixedVerts.find(idx[k]) == fixedVerts.end()) fixed = false;
				if (movedVerts.find(idx[k]) == movedVerts.end()) moved = false;
			}
			if (fixed) {
				faceColors[faceID] = vec3i(0, 0, 255);
			} else if (moved) {
				faceColors[faceID] = vec3i(255, 0, 0);
			}
		}

		vector<vec3i> vertexColors(mpMesh->amount, vec3i(127, 127, 127));
		for (int vertID : fixedVerts) vertexColors[vertID] = vec3i(0, 0, 255);
		for (int vertID : movedVerts) vertexColors[vertID] = vec3i(255, 0, 0);

		vec3 origBBMin, origBBMax;
		vec3 deformBBMin, deformBBMax;
		if (!MeshUtil::computeAABB(*mpMesh, origBBMin, origBBMax)) return false;
		if (!MeshUtil::computeAABB(mDeformedMesh, deformBBMin, deformBBMax)) return false;
		float spacing = (origBBMax[0] - origBBMin[0] + deformBBMax[0] - deformBBMin[0]) / 2 * 1.2f;

		PlyExporter pe;
		//if (!pe.addMesh(&mpMesh->positions, 0, &mpMesh->indices, &faceColors, vec3(-spacing / 2, 0.0f, 0.0f))) return false;
		//if (!pe.addMesh(&mDeformedMesh.positions, 0, &mDeformedMesh.indices, &faceColors, vec3(spacing / 2, 0.0f, 0.0f))) return false;
		if (!pe.addMesh(&mpMesh->indices, &mpMesh->positions, 0, &vertexColors, vec3(-spacing / 2, 0.0f, 0.0f))) return false;
		if (!pe.addMesh(&mDeformedMesh.indices, &mDeformedMesh.positions, 0, &vertexColors, vec3(spacing / 2, 0.0f, 0.0f))) return false;
		if (!pe.output(filePath + "visual-deform-surface.ply")) return false;
	}
	*/
	return true;
}