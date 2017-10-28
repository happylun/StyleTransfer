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

#include "FeatureSaliency.h"

#include <fstream>
#include <iostream>
#include <set>

#include "Mesh/MeshUtil.h"

#include "Sample/SampleSimplePoissonDisk.h"
#include "Sample/SampleUtil.h"

#include "Feature/FeatureAsset.h"

#include "Context/ContextPartGraph.h"

#include "Similarity/SimilarityMetric.h"

#include "Utility/PlyExporter.h"

#include "Data/StyleSynthesisConfig.h"

using namespace StyleSynthesis;

//#define OUTPUT_PROGRESS

FeatureSaliency::FeatureSaliency() {
}

FeatureSaliency::~FeatureSaliency() {
}

bool FeatureSaliency::loadData(
	TTriangleMesh *mesh,
	TSampleSet *samples,
	vector<vector<vector<int>>> *segments,
	FeatureAsset *features,
	ContextPartGraph *graph)
{
	mpMesh = mesh;
	mpSamples = samples;
	mpSegments = segments;
	mpFeatures = features;
	mpGraph = graph;

	return true;
}

bool FeatureSaliency::process() {

	if (!computeNodeDescriptor()) return false;
	if (!computeNodeFeature()) return false;
	if (!enforceSymmetry()) return false;
	if (!computePointSaliency()) return false;

	return true;
}

bool FeatureSaliency::computeNodeDescriptor() {

#ifdef OUTPUT_PROGRESS
	cout << "Computing node descriptor..." << endl;
#endif

	vec3 meshBBMin, meshBBMax;
	if (!MeshUtil::computeAABB(*mpMesh, meshBBMin, meshBBMax)) return false;
	vec3 meshBBExtent(
		max(fabs(meshBBMin[0]), fabs(meshBBMax[0])),
		max(fabs(meshBBMin[1]), fabs(meshBBMax[1])),
		max(fabs(meshBBMin[2]), fabs(meshBBMax[2])));

	// gather nodes of interest
	for (ContextPartGraphNode *node : mpGraph->mRootNode->mChildren) {

		//// use group nodes
		//mNodes.push_back(node);

		// use largest sub-group nodes
		if (node->mChildren.empty()) {
			mNodes.push_back(node);
		} else {
			for (ContextPartGraphNode *child : node->mChildren) {
				mNodes.push_back(child);
			}
		}
	}

	int numNodes = (int)mNodes.size();
	int numFaces = (int)mpMesh->indices.size();

	vector<int> faceMapping(numFaces);
	mNodeDescriptors.resize(numNodes);

#pragma omp parallel for
	for (int nodeID = 0; nodeID < numNodes; nodeID++) {
		
		ContextPartGraphNode *node = mNodes[nodeID];
		vector<double> &nodeDescriptors = mNodeDescriptors[nodeID];
		vector<int> &nodeSegment = (*mpSegments)[node->mPartLevelID][node->mPartSegmentID];

		nodeDescriptors.clear();
		if (nodeSegment.empty()) continue;

		// face mapping

		for (int faceID : nodeSegment) faceMapping[faceID] = nodeID;

		// node mesh

		TTriangleMesh nodeMesh;
		if (!MeshUtil::extractSubMesh(*mpMesh, nodeSegment, nodeMesh)) error("node mesh");

		// sampling

		TSampleSet nodeSamples;
		Eigen::Matrix3Xd sampleMatP, sampleMatN;
		if (true) {
			int numSamples = 1000; // UNDONE: param desired number of samples
			SampleSimplePoissonDisk sspd(&nodeMesh);
			if (!sspd.runSampling(numSamples)) error("sampling node mesh");
			if (!sspd.exportSample(nodeSamples)) error("exporting node samples");
			if (!SampleUtil::buildMatrices(nodeSamples, sampleMatP, sampleMatN)) error("building sample matrices");
		}
		if (nodeSamples.amount == 0) continue; // empty mesh

		////////////////////////// node descriptor //////////////////////////		

		// those codes are duplicates of the REAL node descriptor inside ContextPartGraphNode
		// but this allows us to select only a subset of node descriptors for saliency

		if (true) {

			// position

			Eigen::Vector3d massCenter = sampleMatP.rowwise().mean();
			double lowHeight = sampleMatP.row(1).minCoeff();
			double highHeight = sampleMatP.row(1).maxCoeff();
			double radialDistance = cml::sqrt_safe(cml::sqr(massCenter[0]) + cml::sqr(massCenter[2]));

			nodeDescriptors.push_back(fabs(massCenter[0]) / meshBBExtent[0]);
			nodeDescriptors.push_back(fabs(massCenter[1]) / meshBBExtent[1]);
			nodeDescriptors.push_back(fabs(massCenter[2]) / meshBBExtent[2]);
			nodeDescriptors.push_back(lowHeight);
			nodeDescriptors.push_back(highHeight);
			nodeDescriptors.push_back(radialDistance);
		}

		if (true) {

			// size

			double sizeEpsilon = nodeSamples.radius;

			vec3 bbMin, bbMax;
			if (!MeshUtil::computeAABB(nodeMesh, bbMin, bbMax)) error("node mesh AABB");
			Eigen::AlignedBox3d boundingBox(
				Eigen::Vector3d(vec3d(bbMin).data()),
				Eigen::Vector3d(vec3d(bbMax).data()));
			Eigen::Vector3d boundingBoxSizes = boundingBox.sizes();

			Eigen::Matrix3Xd mat = sampleMatP;
			mat = mat.colwise() - mat.rowwise().mean();
			Eigen::Vector3d axisAlignedVariance = mat.array().square().rowwise().mean().max(sizeEpsilon*sizeEpsilon);

			Eigen::JacobiSVD< Eigen::Matrix3Xd > svd(mat, Eigen::ComputeThinU);
			Eigen::Vector3d vecS = svd.singularValues();
			Eigen::Vector3d principalVariance = vecS.cwiseMax(sizeEpsilon*sizeEpsilon);

			// orientation

			Eigen::Vector3d majorOrientation;
			Eigen::Matrix3d matU = svd.matrixU();
			Eigen::Vector3d aspect(vecS[0] - vecS[1], vecS[1] - vecS[2], vecS[2]);
			if (aspect[0] > aspect[1] && aspect[0] > aspect[2]) {
				// stick-like
				majorOrientation = matU.col(0);
			} else if (aspect[1] > aspect[0] && aspect[1] > aspect[2]) {
				// plane-like
				majorOrientation = matU.col(2);
			} else {
				// sphere-like
				majorOrientation = matU.col(0);
			}
			majorOrientation.normalize();


			nodeDescriptors.push_back(boundingBoxSizes[0]);
			nodeDescriptors.push_back(boundingBoxSizes[1]);
			nodeDescriptors.push_back(boundingBoxSizes[2]);
			nodeDescriptors.push_back(axisAlignedVariance[0]);
			nodeDescriptors.push_back(axisAlignedVariance[1]);
			nodeDescriptors.push_back(axisAlignedVariance[2]);
			nodeDescriptors.push_back(principalVariance[0]);
			nodeDescriptors.push_back(principalVariance[1]);
			nodeDescriptors.push_back(principalVariance[2]);

			nodeDescriptors.push_back(fabs(majorOrientation[0]));
			nodeDescriptors.push_back(fabs(majorOrientation[1]));
			nodeDescriptors.push_back(fabs(majorOrientation[2]));
		}
	}

	// gather samples for nodes

	mNodeSamples.assign(numNodes, vector<int>(0));
	mNodeSampleMapping.resize(mpSamples->amount);

	for (int sampleID = 0; sampleID < mpSamples->amount; sampleID++) {
		int faceID = mpSamples->indices[sampleID];
		int nodeID = faceMapping[faceID];

		mNodeSamples[nodeID].push_back(sampleID);
		mNodeSampleMapping[sampleID] = nodeID;
	}

	// extract symmetry cliques

	if (true) {
		map<int, int> nodeUIDMap; // UID => node ID
		for (int nodeID = 0; nodeID < numNodes; nodeID++) {
			nodeUIDMap[mNodes[nodeID]->mUID] = nodeID;
		}

		vector<bool> nodeFlags(numNodes, false); // visited flag
		for (int nodeID = 0; nodeID < numNodes; nodeID++) {
			if (nodeFlags[nodeID]) continue;
			ContextPartGraphNode *node = mNodes[nodeID];
			nodeFlags[nodeID] = true;

			bool emptyClique = false;
			if (mNodeSamples[nodeID].empty()) emptyClique = true;

			set<int> clique;
			clique.insert(nodeID);
			for (ContextPartGraphNode *otherNode : node->mSymmetry) {
				int uid = otherNode->mUID;
				auto it = nodeUIDMap.find(uid);
				if (it == nodeUIDMap.end()) continue;
				int otherNodeID = it->second;
				clique.insert(otherNodeID);
				nodeFlags[otherNodeID] = true;
				if (mNodeSamples[otherNodeID].empty()) emptyClique = true;
			}
			if (emptyClique) continue; // skip any clique having empty node
			vector<int> cliqueList(clique.begin(), clique.end());
			mNodeSymmetryCliques.push_back(cliqueList);
		}
	}

	return true;
}

bool FeatureSaliency::computeNodeFeature() {

#ifdef OUTPUT_PROGRESS
	cout << "Computing node feature..." << endl;
#endif

	int numSamples = mpSamples->amount;
	int numNodes = (int)mNodeSamples.size();

	// compute curvature features on samples

	vector<vector<double>> sampleCurvatureFeatures(mpSamples->amount);

#pragma omp parallel for
	for (int sampleID = 0; sampleID < numSamples; sampleID++) {
		if (!FeatureSampleCurvature::getMetrics(
			mpFeatures->mCurvature[sampleID],
			sampleCurvatureFeatures[sampleID])) error("curvature metrics");
	}

	// find samples near feature curves

	int numCurveTypes = (int)mpFeatures->mCurve.mPointClouds.size();
	vector<vector<bool>> sampleCurveFlags(numCurveTypes, vector<bool>(numSamples, false));

	if (true) {
		SKDTree sampleTree;
		SKDTreeData sampleTreeData;
		if (!SampleUtil::buildKdTree(mpSamples->positions, sampleTree, sampleTreeData)) return false;

		for (int typeID = 0; typeID < numCurveTypes; typeID++) {

			Eigen::Matrix3Xd curvePointMat;
			if (!SampleUtil::buildMatrices(mpFeatures->mCurve.mPointClouds[typeID], curvePointMat)) return false;

			Eigen::VectorXi curvePointNeighbors;
			if (!SampleUtil::findNearestNeighbors(sampleTree, curvePointMat, curvePointNeighbors)) return false;

			for (int id = 0; id < (int)curvePointNeighbors.size(); id++) {
				sampleCurveFlags[typeID][curvePointNeighbors[id]] = true;
			}
		}
	}

	// compute node features

	mNodeFeatures.resize(numNodes);

	for (int nodeID = 0; nodeID < numNodes; nodeID++) {

		vector<double> &features = mNodeFeatures[nodeID];
		vector<int> &samples = mNodeSamples[nodeID];

		features.clear();
		if (samples.empty()) continue;

		// curvature metrics (average over samples)

		vector<double> curvatureMetrics;
		if (true) {
			int numDimensions = (int)sampleCurvatureFeatures[0].size();
			curvatureMetrics.assign(numDimensions, 0);
			for (int sampleID : samples) {
				for (int dim = 0; dim < numDimensions; dim++) {
					curvatureMetrics[dim] += sampleCurvatureFeatures[sampleID][dim];
				}
			}
			double denom = 1.0 / samples.size();
			for (int dim = 0; dim < numDimensions; dim++) curvatureMetrics[dim] *= denom;
		}

		// percentage of curve points

		vector<double> curvePointPercentage(numCurveTypes);
		if (true) {
			int totalPointCount = (int)samples.size();
			for (int typeID = 0; typeID < numCurveTypes; typeID++) {
				int curvePointCount = 0;
				for (int sampleID : samples) {
					if (sampleCurveFlags[typeID][sampleID]) curvePointCount++;
				}
				curvePointPercentage[typeID] = curvePointCount / (double)totalPointCount;
			}			
		}


		features.insert(features.end(), curvatureMetrics.begin(), curvatureMetrics.end());
		features.insert(features.end(), curvePointPercentage.begin(), curvePointPercentage.end());
	}

	return true;
}

bool FeatureSaliency::enforceSymmetry() {

#ifdef OUTPUT_PROGRESS
	cout << "Enforcing symmetry..." << endl;
#endif

	int numCliques = (int)mNodeSymmetryCliques.size();

	int dimNodeDesc = (int)mNodeDescriptors[0].size();
	int dimNodeFeat = (int)mNodeFeatures[0].size();

	for (int cliqueID = 0; cliqueID < numCliques; cliqueID++) {
		int numNodes = (int)mNodeSymmetryCliques[cliqueID].size();
		Eigen::MatrixXd cliqueDescMat(dimNodeDesc, numNodes);
		Eigen::MatrixXd cliqueFeatMat(dimNodeFeat, numNodes);

		for (int id = 0; id < numNodes; id++) {
			int nodeID = mNodeSymmetryCliques[cliqueID][id];
			for (int dim = 0; dim < dimNodeDesc; dim++) {
				cliqueDescMat(dim, id) = mNodeDescriptors[nodeID][dim];
			}
			for (int dim = 0; dim < dimNodeFeat; dim++) {
				cliqueFeatMat(dim, id) = mNodeFeatures[nodeID][dim];
			}
		}

		Eigen::VectorXd cliqueMeanDesc = cliqueDescMat.rowwise().mean();
		Eigen::VectorXd cliqueMeanFeat = cliqueFeatMat.rowwise().mean();

		for (int id = 0; id < numNodes; id++) {
			int nodeID = mNodeSymmetryCliques[cliqueID][id];
			for (int dim = 0; dim < dimNodeDesc; dim++) {
				mNodeDescriptors[nodeID][dim] = cliqueMeanDesc[dim];
			}
			for (int dim = 0; dim < dimNodeFeat; dim++) {
				mNodeFeatures[nodeID][dim] = cliqueMeanFeat[dim];
			}
		}
	}

	return true;
}

bool FeatureSaliency::computePointSaliency() {

#ifdef OUTPUT_PROGRESS
	cout << "Computing point saliency..." << endl;
#endif

	mPointSaliency.resize(mpSamples->amount);

	for (int sampleID = 0; sampleID < mpSamples->amount; sampleID++) {

		auto &saliency = mPointSaliency[sampleID];
		saliency.clear();

		int nodeID = mNodeSampleMapping[sampleID];
		vector<double> &nodeDesc = mNodeDescriptors[nodeID];
		vector<double> &nodeFeat = mNodeFeatures[nodeID];

		if (StyleSynthesisConfig::mStyle_UseContextualSaliencyFeature) {
			saliency.insert(saliency.end(), nodeDesc.begin(), nodeDesc.end());
			//saliency.insert(saliency.end(), nodeFeat.begin(), nodeFeat.end());
		} else {
			saliency.insert(saliency.end(), nodeFeat.begin(), nodeFeat.end());
		}
	}

	return true;
}

bool FeatureSaliency::output(Eigen::MatrixXd &saliency) {

	int numPoints = (int)mPointSaliency.size();
	int numDimensions = (int)mPointSaliency[0].size();
	saliency.resize(numPoints, numDimensions);

	for (int pointID = 0; pointID < numPoints; pointID++) {
		for (int dim = 0; dim < numDimensions; dim++) {
			saliency(pointID, dim) = mPointSaliency[pointID][dim];
		}
	}

	return true;
}

bool FeatureSaliency::visualize(string fileName, TSampleSet &samples, Eigen::MatrixXd &saliency) {

	int numSamples = samples.amount;

	vector<double> saliencyValues(numSamples);

#pragma omp parallel for
	for (int sampleID = 0; sampleID < numSamples; sampleID++) {
		Eigen::VectorXd salVec = saliency.row(sampleID);
		saliencyValues[sampleID] = SimilarityMetric::evaluateSaliency(salVec);

		/*
#pragma omp critical
		if (true) {
			vec3 pos = samples.positions[sampleID];
			if (pos[0] > -0.401f && pos[0] < -0.4f &&
				pos[1] > 0.466f && pos[1] < 0.467f &&
				pos[2] > 0.304f && pos[2] < 0.305f)
			{
				Eigen::VectorXd normalizedVec = salVec.cwiseQuotient(SimilarityMetric::mScaleSaliency);
				normalizedVec = normalizedVec.cwiseMax(-1.0);
				normalizedVec = normalizedVec.cwiseMin(1.0);
				Eigen::VectorXd homoVec(normalizedVec.size() + 1);
				homoVec << normalizedVec, 1.0;
				Eigen::VectorXd salTerms = homoVec.array() * SimilarityMetric::mWeightsSaliency.array();
				Eigen::MatrixXd mat(salTerms.size(), 3);
				mat << homoVec, SimilarityMetric::mWeightsSaliency, salTerms;
				cout << "========== seat ===========" << endl;
				cout << "Saliency computation:" << endl;
				cout << mat << endl;
				cout << "Before sigmoid: " << homoVec.dot(SimilarityMetric::mWeightsSaliency) << endl;
				cout << "Saliency value: " << saliencyValues[sampleID] << endl;
			}
			if (pos[0] > -0.831f && pos[0] < -0.83f &&
				pos[1] > 0.054f && pos[1] < 0.055f &&
				pos[2] > 0.401f && pos[2] < 0.402f)
			{
				Eigen::VectorXd normalizedVec = salVec.cwiseQuotient(SimilarityMetric::mScaleSaliency);
				normalizedVec = normalizedVec.cwiseMax(-1.0);
				normalizedVec = normalizedVec.cwiseMin(1.0);
				Eigen::VectorXd homoVec(normalizedVec.size() + 1);
				homoVec << normalizedVec, 1.0;
				Eigen::VectorXd salTerms = homoVec.array() * SimilarityMetric::mWeightsSaliency.array();
				Eigen::MatrixXd mat(salTerms.size(), 3);
				mat << homoVec, SimilarityMetric::mWeightsSaliency, salTerms;
				cout << "========== leg ===========" << endl;
				cout << "Saliency computation:" << endl;
				cout << mat << endl;
				cout << "Before sigmoid: " << homoVec.dot(SimilarityMetric::mWeightsSaliency) << endl;
				cout << "Saliency value: " << saliencyValues[sampleID] << endl;
			}
		}
		*/
	}

	vector<vec3i> colors(numSamples);
	double colorScaleFactor = 50.0;
	for (int sampleID = 0; sampleID < numSamples; sampleID++) {
		double v = cml::clamp(saliencyValues[sampleID]*colorScaleFactor, 0.0, 1.0);
		//colors[sampleID] = FeatureUtil::colorMapping(v); // jet color
		colors[sampleID] = vec3i(255, (int)(255 * (1 - v)), (int)(255 * (1 - v))); // white to red
	}

	PlyExporter pe;
	if (!pe.addPoint(&samples.positions, &samples.normals, &colors)) return false;
	if (!pe.output(fileName)) return false;

	return true;
}