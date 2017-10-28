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

#include "CurveDeformVertex.h"

#include <fstream>

#include "Library/libiglHelperBegin.h"
#include <igl/tetgen/tetrahedralize.h>
#include "Library/libiglHelperEnd.h"

#include "Mesh/MeshUtil.h"
#include "Sample/SampleUtil.h"
#include "Curve/CurveUtil.h"

#include "Match/MatchCurveICP.h"

#include "Deform/DeformVolumeARAP.h"
#include "Deform/DeformVolumetricGraph.h"

#include "Utility/PlyExporter.h"
#include "Utility/FileUtil.h"
#include "Utility/Timer.h"

#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

#define OUTPUT_VISUALIZATION

CurveDeformVertex::CurveDeformVertex() {
}

CurveDeformVertex::~CurveDeformVertex() {
}

bool CurveDeformVertex::loadMesh(TTriangleMesh &mesh) {

	mpMesh = &mesh;
	return true;
}

bool CurveDeformVertex::loadCurves(
	vector<vector<vector<vec3>>> &exemplarCurvesGroup,
	vector<vector<vector<vec3>>> &candidateCurvesGroup)
{
	mExemplarCurvesGroup = exemplarCurvesGroup;
	mCandidateCurvesGroup = candidateCurvesGroup;
	return true;
}

bool CurveDeformVertex::loadNames(string visualizeFolder) {

	mVisualizeFolder = visualizeFolder;
	return true;
}

bool CurveDeformVertex::process() {

	mNoDeformation = (StyleSynthesisConfig::mDeform_Method == "None");

	if (!alignCurves()) return false;
	if (!extractHandles()) return false;
	if (!deformMesh()) return false;

	return true;
}

bool CurveDeformVertex::alignCurves() {

	cout << "Aligning matched curves..." << endl;

	mExemplarCurvePoints.clear();
	mCandidateCurvePoints.clear();
	mInterpolatedCurveFlags.clear();

	int numGroups = (int)mExemplarCurvesGroup.size();
	for (int groupID = 0; groupID < numGroups; groupID++) {

		cout << "Aligning curve group " << (groupID + 1) << " / " << numGroups << endl;

		vector<vector<vec3>> &exemGroup = mExemplarCurvesGroup[groupID];
		vector<vector<vec3>> &candGroup = mCandidateCurvesGroup[groupID];

		if (!alignGroup(exemGroup)) return false;
		if (!alignGroup(candGroup)) return false;

		int numMatchings = (int)exemGroup.size();

		vector<vector<vec3d>> groupExemPoints(0);
		vector<vector<vec3d>> groupCandPoints(0);
		groupExemPoints.reserve(numMatchings);
		groupCandPoints.reserve(numMatchings);

		for (int matchID = 0; matchID < numMatchings; matchID++) {
			vector<vec3> &exemCurve = exemGroup[matchID];
			vector<vec3> &candCurve = candGroup[matchID];

			vector<vec3d> alignedExemPoints;
			vector<vec3d> alignedCandPoints;
			if (!alignCurve(exemCurve, candCurve, alignedExemPoints, alignedCandPoints)) return false;
			if (alignedExemPoints.empty() || alignedCandPoints.empty()) continue;

			groupExemPoints.push_back(alignedExemPoints);
			groupCandPoints.push_back(alignedCandPoints);
		}

		int numOriginalCurves = (int)groupExemPoints.size();
		if (numOriginalCurves <= 0) continue;
		if (!interpolateGroup(groupExemPoints, groupCandPoints)) return false;
		int numInterpolatedCurves = (int)groupExemPoints.size() - numOriginalCurves;

		mExemplarCurvePoints.insert(mExemplarCurvePoints.end(), groupExemPoints.begin(), groupExemPoints.end());
		mCandidateCurvePoints.insert(mCandidateCurvePoints.end(), groupCandPoints.begin(), groupCandPoints.end());
		for (int k = 0; k < numOriginalCurves; k++) mInterpolatedCurveFlags.push_back(false);
		for (int k = 0; k < numInterpolatedCurves; k++) mInterpolatedCurveFlags.push_back(true);
	}

#ifdef OUTPUT_VISUALIZATION
	if (!mNoDeformation) {
		// visualize curve alignment

		cout << "Visualizing curve alignment" << endl;

		vec3i baseColor(127, 127, 127);
		vec3i exemplarColor(255, 0, 0);
		vec3i candidateColor(0, 255, 0);
		double tubeRadius = StyleSynthesisConfig::mCurve_VisualizationTubeRadius;

		PlyExporter peAll;
		PlyExporter peBase;
		PlyExporter peExemplar;
		PlyExporter peCandidate;
		PlyExporter peLink;

		int numMatchings = (int)mExemplarCurvePoints.size();
		for (int matchID = 0; matchID < numMatchings; matchID++) {

			vector<vec3> exemCurve;
			vector<vec3> candCurve;
			for (vec3d point : mExemplarCurvePoints[matchID]) {
				exemCurve.push_back(vec3(point));
			}
			for (vec3d point : mCandidateCurvePoints[matchID]) {
				candCurve.push_back(vec3(point));
			}

			TTriangleMesh exemTube, candTube;
			if (!CurveUtil::makeTube(exemCurve, exemTube, (float)tubeRadius)) return false;
			if (!CurveUtil::makeTube(candCurve, candTube, (float)tubeRadius)) return false;

			if (!peExemplar.addMesh(&exemTube.positions, &exemTube.normals, &exemTube.indices, exemplarColor)) return false;
			if (!peAll.addMesh(&exemTube.positions, &exemTube.normals, &exemTube.indices, exemplarColor)) return false;

			if (!peCandidate.addMesh(&candTube.positions, &candTube.normals, &candTube.indices, candidateColor)) return false;
			if (!peAll.addMesh(&candTube.positions, &candTube.normals, &candTube.indices, candidateColor)) return false;
			
			int numPoints = (int)exemCurve.size();
			vector<vec2i> indices(numPoints-1);
			for (int k = 0; k < numPoints - 1; k++) indices[k] = vec2i(k, k + 1);
			if (!peLink.addLine(&indices, &exemCurve)) return false;
			if (!peLink.addLine(&indices, &candCurve)) return false;
			if (!peLink.addLine(&exemCurve, &candCurve)) return false;
		}

		if (!peBase.addMesh(&mpMesh->positions, &mpMesh->normals, &mpMesh->indices, baseColor)) return false;
		if (!peAll.addMesh(&mpMesh->positions, &mpMesh->normals, &mpMesh->indices, baseColor)) return false;

		if (!peAll.output(mVisualizeFolder + "match-all.ply")) return false;
		if (!peExemplar.output(mVisualizeFolder + "match-exemplar.ply")) return false;
		if (!peCandidate.output(mVisualizeFolder + "match-candidate.ply")) return false;
		if (!peBase.output(mVisualizeFolder + "match-base.ply")) return false;
		if (!peLink.output(mVisualizeFolder + "match-link.ply")) return false;
	}
#endif

	return true;
}

bool CurveDeformVertex::alignGroup(vector<vector<vec3>> &group) {

	// align curves within a group to have consistent orientations

	if (group.empty()) return true;

	vec3 pivotDir = group[0].front() - group[0].back();
	if (cml::dot(pivotDir, vec3(1.0f, 1.0f, 1.0f)) < 0) pivotDir = -pivotDir;

	for (vector<vec3> &curve : group) {
		vec3 currentDir = curve.front() - curve.back();
		if (cml::dot(pivotDir, currentDir) < 0) {
			// flip it!
			vector<vec3> flippedCurve(curve.rbegin(), curve.rend());
			curve.swap(flippedCurve);
		}
	}

	return true;
}

bool CurveDeformVertex::alignCurve(
	vector<vec3> &exemCurve,
	vector<vec3> &candCurve,
	vector<vec3d> &alignedExemPoints,
	vector<vec3d> &alignedCandPoints)
{

	double sampleRadius = StyleSynthesisConfig::mCurve_SamplingRadius;

	alignedExemPoints.clear();
	alignedCandPoints.clear();

	// get sample points on curve

	vector<vec3> exemSamples, candSamples;
	if (true) {
		if (!CurveUtil::sampleLine(exemCurve, exemSamples, (float)sampleRadius)) return false;
		double minLen = DBL_MAX;
		for (int id = 1; id < (int)exemSamples.size(); id++) {
			double len = (double)(exemSamples[id] - exemSamples[id - 1]).length();
			if (len > 0 && len < minLen) minLen = len;
		}
		if (!CurveUtil::sampleLine(candCurve, candSamples, (float)(max(0.0001, minLen * 0.5)))) return false;
	}
	if ((int)exemSamples.size() < 2 || (int)candSamples.size() < 2) return true; // empty curve

	// align curve

	Eigen::Matrix3Xd matEP, matCP;
	if (true) {
		if (!CurveUtil::buildMatrix(exemSamples, matEP)) return false;
		if (!CurveUtil::buildMatrix(candSamples, matCP)) return false;
		Eigen::Affine3d transform;
		transform.setIdentity();
#ifdef OUTPUT_VISUALIZATION
		if (!MatchCurveICP::visualize(mVisualizeFolder + "align-before.ply", matEP, matCP, transform)) return false;
#endif
		if (!MatchCurveICP::prealign(matEP, matCP, transform)) return false;
#ifdef OUTPUT_VISUALIZATION
		if (!MatchCurveICP::visualize(mVisualizeFolder + "align-after.ply", matEP, matCP, transform)) return false;
#endif
		matEP = transform * matEP;
	}

	// align end point

	if (true) {
		Eigen::Vector3d vecEP1 = matEP.leftCols(1);
		Eigen::Vector3d vecEP2 = matEP.rightCols(1);
		Eigen::Vector3d vecCP1 = matCP.leftCols(1);
		Eigen::Vector3d vecCP2 = matCP.rightCols(1);
		double dist1 = min((vecEP1 - vecCP1).norm(), (vecEP2 - vecCP2).norm());
		double dist2 = min((vecEP1 - vecCP2).norm(), (vecEP2 - vecCP1).norm());
		if (dist2 < dist1) { // flip exemplar curve
			matEP = matEP.rowwise().reverse().eval();
		}
		if (StyleSynthesisConfig::mDeform_AlignCurveEndPoints) {
			Eigen::Affine3d transform;
			if (!MatchCurveICP::alignTwoPoints(vecEP1, vecEP2, vecCP1, vecCP2, transform)) return false;
#ifdef OUTPUT_VISUALIZATION
			if (!MatchCurveICP::visualize(mVisualizeFolder + "align-end.ply", matEP, matCP, transform)) return false;
#endif
			matEP = transform * matEP;
		}
	}

	// find curve point correspondences

	if (true) {
		Eigen::VectorXi matchIndices;
		if (!MatchCurveICP::findCorrespondences(matEP, matCP, matchIndices)) return false;
		Eigen::Matrix3Xd slicedEP;
		if (!SampleUtil::sliceMatrices(matEP, matchIndices, slicedEP)) return false;
		matEP.swap(slicedEP);
	}

	// output curve

	alignedExemPoints.resize((int)matEP.cols());
	alignedCandPoints.resize((int)matCP.cols());
	for (int id = 0; id < (int)matEP.cols(); id++) {
		alignedExemPoints[id] = vec3d(matEP(0, id), matEP(1, id), matEP(2, id));
	}
	for (int id = 0; id < (int)matCP.cols(); id++) {
		alignedCandPoints[id] = vec3d(matCP(0, id), matCP(1, id), matCP(2, id));
	}

	return true;
}

bool CurveDeformVertex::interpolateGroup(
	vector<vector<vec3d>> &exemGroup,
	vector<vector<vec3d>> &candGroup)
{

	// compute curve center

	int numMatchings = (int)candGroup.size();
	vector<vec3d> candCurveCenter(numMatchings);
	for (int matchID = 0; matchID < numMatchings; matchID++) {
		vector<vec3d> &curve = candGroup[matchID];
		vec3d center(0.0, 0.0, 0.0);
		for (vec3d point : curve) center += point;
		if (!curve.empty()) center *= 1.0 / curve.size();
		candCurveCenter[matchID] = center;
	}

	// find interpolation pairs

	vector<vec2i> interPairs(0);
	if (true) {
		// HACK: find interpolation loop within group

		// find loop by searching next nearest neighbor
		vector<int> interLoop(0);
		vector<bool> visitedFlag(numMatchings, false);
		interLoop.push_back(0);
		visitedFlag[0] = true;
		for (int iter = 0; iter < numMatchings - 1; iter++) {
			int currentID = interLoop.back();
			int nextID = -1;
			double nextDist = DBL_MAX;
			for (int otherID = 0; otherID < numMatchings; otherID++) {
				if (visitedFlag[otherID]) continue;
				double dist = (candCurveCenter[currentID] - candCurveCenter[otherID]).length();
				if (dist < nextDist) {
					nextDist = dist;
					nextID = otherID;
				}
			}
			interLoop.push_back(nextID);
			visitedFlag[nextID] = true;
		}

#ifdef OUTPUT_VISUALIZATION
		if (true) {
			cout << "Interpolation loop:";
			vector<vec3> loopPoints;
			for (int id = 0; id < (int)interLoop.size(); id++) {
				int curveID = interLoop[id];
				cout << " " << curveID;
				loopPoints.push_back(candCurveCenter[curveID]);
			}
			cout << endl;
			vector<vec2i> loopIndices(interLoop.size());
			for (int k = 0; k < (int)interLoop.size(); k++) loopIndices[k] = vec2i(k, (k + 1) % interLoop.size());

			PlyExporter pe;
			if (!pe.addLine(&loopIndices, &loopPoints)) return false;
			if (!pe.output(mVisualizeFolder + "inter-loop.ply")) return false;
		}
#endif

		// generate loop
		for (int id = 1; id < (int)interLoop.size(); id++) {
			interPairs.push_back(vec2i(interLoop[id - 1], interLoop[id]));
		}
		if ((int)interLoop.size() > 2) interPairs.push_back(vec2i(interLoop.back(), interLoop.front()));
	}

	// do interpolation

	double interMaxGap = StyleSynthesisConfig::mDeform_MeshSubdivisionRadius
		* StyleSynthesisConfig::mDeform_InterpolationInterval;

	int numInterPairs = (int)interPairs.size();
	for (int pairID = 0; pairID < numInterPairs; pairID++) {
		int sourceID = interPairs[pairID][0];
		int targetID = interPairs[pairID][1];
		cout << "Interpolating " << sourceID << ", " << targetID << endl;
		int numSteps = (int)((candCurveCenter[sourceID] - candCurveCenter[targetID]).length() / interMaxGap);

		for (int interID = 1; interID < numSteps; interID++) {
			double alpha = interID / (double)numSteps;
			vector<vec3d> interExemCurve;
			vector<vec3d> interCandCurve;
			if (!interpolateCurve(alpha, exemGroup[sourceID], exemGroup[targetID], interExemCurve)) return false;
			if (!interpolateCurve(alpha, candGroup[sourceID], candGroup[targetID], interCandCurve)) return false;
			exemGroup.push_back(interExemCurve);
			candGroup.push_back(interCandCurve);
		}
	}

	return true;
}

bool CurveDeformVertex::interpolateCurve(
	double alpha,
	vector<vec3d> &sourceCurve,
	vector<vec3d> &targetCurve,
	vector<vec3d> &interpolatedCurve)
{

	int numSourcePoints = (int)sourceCurve.size();
	int numTargetPoints = (int)targetCurve.size();

	// extract mapped curve (on target curve)
	vector<vec3d> mappedCurve(numSourcePoints);
	for (int id = 0; id < numSourcePoints; id++) {
		double factor = id * (numTargetPoints - 1.0) / (numSourcePoints - 1.0);
		int prevID = cml::clamp((int)(floor(factor)), 0, numTargetPoints - 1);
		int nextID = cml::clamp((int)(ceil(factor)), 0, numTargetPoints - 1);
		double beta = factor - prevID;
		mappedCurve[id] = targetCurve[prevID] * (1.0 - beta) + targetCurve[nextID] * beta;
	}

	// compute interpolated points
	interpolatedCurve.resize(numSourcePoints);
	for (int id = 0; id < numSourcePoints; id++) {
		 interpolatedCurve[id] = sourceCurve[id] * (1.0 - alpha) + mappedCurve[id] * alpha;
	}

	return true;
}

bool CurveDeformVertex::extractHandles() {

	cout << "Extracting deformation handles..." << endl;

	double impactRadius = StyleSynthesisConfig::mDeform_CurveImpactRadius;
	double minImpactRadius = StyleSynthesisConfig::mDeform_MeshSubdivisionRadius * 0.5;

	// get all points from all curves

	int numCurves = (int)mCandidateCurvePoints.size();
	vector<vec3> curvePoints; // sample points from all curves
	vector<vec2i> curveIndices; // (curve ID, point ID) : # of total curve sampel points
	for (int curveID = 0; curveID < numCurves; curveID++) {
		auto &curve = mCandidateCurvePoints[curveID];
		int numPoints = (int)curve.size();
		for (int pointID = 0; pointID < numPoints; pointID++) {
			curvePoints.push_back(curve[pointID]);
			curveIndices.push_back(vec2i(curveID, pointID));
		}
	}

	// build Kd tree for curve points

	SKDTree curveTree;
	SKDTreeData curveTreeData;
	if (!SampleUtil::buildKdTree(curvePoints, curveTree, curveTreeData)) return false;

	// compute curve handle displacement

	mDeformCurveHandle.assign(mpMesh->amount, false);
	mDeformHardHandle.assign(mpMesh->amount, false);
	mDeformDisplacement.assign(mpMesh->amount, vec3d(0.0, 0.0, 0.0)); // by default 0 displacement

#pragma omp parallel for
	for (int vertID = 0; vertID < mpMesh->amount; vertID++) {
		vec3 point = mpMesh->positions[vertID];
		SKDT::NamedPoint queryPoint(point[0], point[1], point[2]);
		int numResults = mNoDeformation ? 1 : 100;
		Thea::BoundedSortedArray<SKDTree::Neighbor> queryResult(numResults);
		curveTree.kClosestElements<Thea::MetricL2>(queryPoint, queryResult, impactRadius);
		if (!queryResult.isEmpty()) {
			double accumWeights = 0;
			vec3d accumDisplacement(0.0, 0.0, 0.0);
			for (int qID = 0; qID < queryResult.size(); qID++) {
				int nnID = (int)curveTree.getElements()[queryResult[qID].getIndex()].id;
				int curveID = curveIndices[nnID][0];
				int curvePointID = curveIndices[nnID][1];
				vec3d exemPoint = mExemplarCurvePoints[curveID][curvePointID];
				vec3d candPoint = mCandidateCurvePoints[curveID][curvePointID];
				vec3d displacement = exemPoint - candPoint;
				double weight = 1.0 / cml::sqr(max(minImpactRadius, (candPoint - vec3d(point)).length()));
				accumDisplacement += displacement * weight;
				accumWeights += weight;
				if (qID == 0 && mInterpolatedCurveFlags[curveID]) {
					mDeformHardHandle[vertID] = true;
					break;
				}
			}
			if (accumWeights) {
				mDeformDisplacement[vertID] = accumDisplacement / accumWeights;
				mDeformCurveHandle[vertID] = true;
			}
		}
	}

	// find fixed vertices

	mDeformFixedHandle.assign(mpMesh->amount, false);
	if (true) {
		// mark vertices lying on bounding box sides as fixed
		vec3 bbMin, bbMax;
		if (!MeshUtil::computeAABB(*mpMesh, bbMin, bbMax)) return false;
		vector<bool> bbMinFlag(3, true);
		vector<bool> bbMaxFlag(3, true);
		for (int vertID = 0; vertID < mpMesh->amount; vertID++) {
			if (!mDeformCurveHandle[vertID]) continue;
			vec3 v = mpMesh->positions[vertID];
			for (int k = 0; k < 3; k++) {
				if (v[k] == bbMin[k]) bbMinFlag[k] = false; // this side should not be fixed
				if (v[k] == bbMax[k]) bbMaxFlag[k] = false;
			}
		}
		/*
		// hack: fix ceiling and floor
		bbMinFlag[0] = bbMaxFlag[0] = false;
		bbMinFlag[1] = bbMaxFlag[1] = true;
		bbMinFlag[2] = bbMaxFlag[2] = false;
		*/
		for (int vertID = 0; vertID < mpMesh->amount; vertID++) {
			if (mDeformCurveHandle[vertID]) continue;
			vec3 v = mpMesh->positions[vertID];
			for (int k = 0; k < 3; k++) {
				if ((bbMinFlag[k] && v[k] == bbMin[k]) ||
					(bbMaxFlag[k] && v[k] == bbMax[k]))
				{
					mDeformFixedHandle[vertID] = true;
					break;
				}
			}
		}
	}

	return true;
}

bool CurveDeformVertex::deformMesh() {

	cout << "Deforming mesh..." << endl;

	int numVertices = mpMesh->amount;

	int numSteps = StyleSynthesisConfig::mDeform_Steps;

	// compute constraints for deformation handles

	map<int, vec3d> constraints;
	for (int vertexID = 0; vertexID < numVertices; vertexID++) {
		if (mDeformCurveHandle[vertexID]) {
			vec3d displacement = mDeformDisplacement[vertexID];
			displacement /= numSteps;
			constraints[vertexID] = displacement;
		}
		if (mDeformFixedHandle[vertexID]) {
			constraints[vertexID] = vec3d(0.0, 0.0, 0.0);
		}
	}

	// do deformation

	mDeformedMesh = *mpMesh;

	if (mNoDeformation) {
		// just move deformation handles
		for (auto &it : constraints) {
			int vertID = it.first;
			vec3d displacement = it.second;
			if (vertID < mDeformedMesh.amount) {
				mDeformedMesh.positions[vertID] += vec3(displacement);
			}
		}

		if (!smoothMesh(mDeformedMesh, mDeformHardHandle)) return false;

	} else {
		for (int step = 0; step < numSteps; step++) {

			Timer::tic();

			cout << "=== Step " << (step + 1) << " ===" << endl;

			string stepVisualizeFolder = mVisualizeFolder + "step-" + to_string(step + 1) + "/";
			if (!FileUtil::makedir(stepVisualizeFolder)) return false;

			TTetrahedralMesh tetMesh;
			DeformVolumetricGraph dvg;
			if (!dvg.loadMesh(mDeformedMesh)) return false;
			if (!dvg.loadNames(stepVisualizeFolder)) return false;
			if (!dvg.process()) return false;
			if (!dvg.output(tetMesh)) return false;

			TTetrahedralMesh deformedTetMesh;
			DeformVolumeARAP dva;
			if (!dva.loadMesh(tetMesh)) return false;
			if (!dva.loadConstraints(constraints)) return false;
			if (!dva.process()) return false;
			if (!dva.output(deformedTetMesh)) return false;
#ifdef OUTPUT_VISUALIZATION
			if (!dva.visualize(stepVisualizeFolder + "vertex-")) return false;
#endif

			for (int vertID = 0; vertID < mpMesh->amount; vertID++) {
				mDeformedMesh.positions[vertID] = deformedTetMesh.positions[vertID];
			}
			mDeformedMesh.amount = (int)mDeformedMesh.positions.size();

#ifdef OUTPUT_VISUALIZATION
			if (!MeshUtil::saveMesh(stepVisualizeFolder + "deformed.ply", mDeformedMesh)) return false;
#endif

			if (!smoothMesh(mDeformedMesh, mDeformHardHandle)) return false;

#ifdef OUTPUT_VISUALIZATION
			if (!MeshUtil::saveMesh(stepVisualizeFolder + "step-" + to_string(step + 1) + ".ply", mDeformedMesh)) return false;
#endif

			cout << "Step time: " << Timer::toString(Timer::toc());
		}
	}

	return true;
}

bool CurveDeformVertex::smoothMesh(TTriangleMesh &mesh, vector<bool> &constraintFlags) {

	// constrained Taubin smoothing

	cout << "Smoothing mesh" << endl;

	int numIteration = StyleSynthesisConfig::mDeform_PostSmoothingIterations;
	double lambda = 0.5;
	double mu = -0.53;

	// find neighbors and get edge info

	int numVertices = mesh.amount;
	vector<set<int>> vertexNeighbros(numVertices);
	map<vec2i, set<int>> edgeMap;

	for (vec3i idx : mesh.indices) {
		for (int k = 0; k < 3; k++) {
			int v1 = idx[k];
			int v2 = idx[(k + 1) % 3];
			int vt = idx[(k + 2) % 3];
			vertexNeighbros[v1].insert(v2);
			vertexNeighbros[v2].insert(v1);

			vec2i key = v1 < v2 ? vec2i(v1, v2) : vec2i(v2, v1);
			edgeMap[key].insert(vt);
		}
	}

	// mark fixed vertices

	vector<bool> fixedFlag(numVertices, false);
	for (int vertID = 0; vertID < numVertices; vertID++) {
		if (constraintFlags[vertID]) {
			// deformation handle
			fixedFlag[vertID] = true;
		}
	}
	for (auto &edge : edgeMap) {
		if (edge.second.size() == 1) {
			// boundary points
			fixedFlag[edge.first[0]] = true;
			fixedFlag[edge.first[1]] = true;
		}
	}

	// Taubin smoothing

	for (int iter = 0; iter < numIteration; iter++) {
		vector<vec3> newVertices = mesh.positions;
		for (int alterID = 0; alterID < 2; alterID++) {
#pragma omp parallel for
			for (int vertID = 0; vertID < numVertices; vertID++) {
				if (fixedFlag[vertID]) continue;
				vec3d laplacian(0.0, 0.0, 0.0);
				double totalWeight = 0;
				for (int nbID : vertexNeighbros[vertID]) {
					vec2i key = vertID < nbID ? vec2i(vertID, nbID) : vec2i(nbID, vertID);
					auto &edge = edgeMap.find(key);

					double weight = 0;
					//// cotangent weight
					//for (int oppoID : edge->second) {					
					//	vec3d v1 = mesh.positions[vertID] - mesh.positions[oppoID];
					//	vec3d v2 = mesh.positions[nbID] - mesh.positions[oppoID];
					//	double cotWeight = fabs(cml::dot(v1, v2) / cml::cross(v1, v2).length());
					//	weight += cotWeight;
					//}

					// uniform weight
					weight = 1.0;

					vec3d delta = mesh.positions[nbID] - mesh.positions[vertID];
					laplacian += delta * weight;
					totalWeight += weight;
				}
				if (totalWeight) laplacian *= 1.0 / totalWeight;
				vec3d displacement = laplacian * (alterID == 0 ? lambda : mu); // first shrink then inflate
				newVertices[vertID] = mesh.positions[vertID] + vec3(displacement);
			}
			mesh.positions = newVertices;
		}
	}

	if (!MeshUtil::recomputeNormals(mesh)) return false;

	//if (true) {
	//	vector<vec3i> vertexColor(numVertices, vec3i(127, 127, 127));
	//	for (int vertID = 0; vertID < numVertices; vertID++) {
	//		if (fixedFlag[vertID]) vertexColor[vertID] = vec3i(255, 0, 0);
	//	}
	//	PlyExporter pe;
	//	if (!pe.addMesh(&mesh.indices, &mesh.positions, &mesh.normals, &vertexColor)) return false;
	//	if (!pe.output("smooth-color.ply")) return false;
	//	system("pause");
	//}

	return true;
}

bool CurveDeformVertex::visualize() {

	// compute spacing

	vec3 originalBBMin, originalBBMax;
	vec3 deformedBBMin, deformedBBMax;
	if (!MeshUtil::computeAABB(*mpMesh, originalBBMin, originalBBMax)) return false;
	if (!MeshUtil::computeAABB(mDeformedMesh, deformedBBMin, deformedBBMax)) return false;
	float spacing = (originalBBMax[0] - originalBBMin[0] + deformedBBMax[0] - deformedBBMin[0]) / 2 * 1.2f;

	// compute vertex colors

	vector<vec3i> vertexColors;
	for (int vertID = 0; vertID < mpMesh->amount; vertID++) {
		if (mDeformCurveHandle[vertID]) {
			if (mDeformDisplacement[vertID].length_squared() == 0) {
				vertexColors.push_back(vec3i(0, 0, 255));
			} else {
				vertexColors.push_back(vec3i(255, 0, 0));
			}
		} else if (mDeformFixedHandle[vertID]) {
			vertexColors.push_back(vec3i(0, 0, 255));
		} else {
			vertexColors.push_back(vec3i(127, 127, 127));
		}
	}

	PlyExporter pe;
	if (!pe.addMesh(&mpMesh->indices, &mpMesh->positions, 0, &vertexColors, vec3(-spacing / 2, 0.0f, 0.0f))) return false;
	if (!pe.addMesh(&mDeformedMesh.indices, &mDeformedMesh.positions, 0, &vertexColors, vec3(spacing / 2, 0.0f, 0.0f))) return false;
	if (!pe.output(mVisualizeFolder + "visual-deform-vertex.ply")) return false;

	if (!MeshUtil::saveMesh(mVisualizeFolder + "visual-deform-shape.ply", mDeformedMesh)) return false;

	return true;
}

bool CurveDeformVertex::output(TTriangleMesh &mesh) {

	mesh = mDeformedMesh;

	return true;
}