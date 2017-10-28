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

#include "PipelineFeatureIO.h"

#include <iostream>
#include <fstream>

#include "Mesh/MeshUtil.h"

#include "Segment/SegmentGroupApxCvx.h"
#include "Segment/SegmentUtil.h"

#include "Sample/SamplePoissonDisk.h"
#include "Sample/SampleSimplePoissonDisk.h"
#include "Sample/SampleUtil.h"

#include "Feature/FeatureAsset.h"
#include "Feature/FeatureSaliency.h"

#include "Context/ContextPartGraph.h"
#include "Context/ContextPartGraphNodeGenerator.h"

#include "Utility/Timer.h"

#include "Data/DataUtil.h"
#include "Data/StyleSynthesisConfig.h"

using namespace std;
using namespace StyleSynthesis;

#define OUTPUT_PROGRESS

bool PipelineFeatureIO::process() {

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string meshListFileName = datasetPrefix + "mesh/mesh-list-all.txt";
	vector<string> meshNameList;
	ifstream meshListFile(meshListFileName);
	while (!meshListFile.eof()) {
		string line;
		getline(meshListFile, line);
		if (line.empty()) break;
		meshNameList.push_back(StringUtil::trim(line));
	}
	meshListFile.close();

	int numMesh = (int)meshNameList.size();

	for (int meshID = 0; meshID < numMesh; meshID++) {
		string meshName = meshNameList[meshID];
		cout << "Processing feature " << meshName << endl;

		if (!runSingleModelSample(meshName)) return error("single model sample error");
		if (!runSingleModelFeature(meshName)) return error("single model feature error");
		if (!runSingleModelPartFeature(meshName)) return error("single model part feature error");
		if (!runSingleModelSaliency(meshName, meshID)) return error("single model saliency error");
	}
	return true;
}

bool PipelineFeatureIO::runSingleModelSample(string meshName) {

	// names...

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string wholeMeshName = datasetPrefix + "segment/" + meshName + ".ply";
	if (!FileUtil::existsfile(wholeMeshName)) {
		cout << "Error: run segmentation pipeline first" << endl;
		return false;
	}

	string sampleName = datasetPrefix + "sample/" + meshName + ".ply";
	if (!FileUtil::makedir(sampleName)) return false;

	if (FileUtil::existsfile(sampleName)) return true; // early quit

	// data

	TTriangleMesh mesh;
	TSampleSet sample;

	// algorithm

#ifdef OUTPUT_PROGRESS
	cout << "Sampling..." << endl;
#endif

	if (!MeshUtil::loadMesh(wholeMeshName, mesh)) return false;

	double sampleRadius = StyleSynthesisConfig::mSample_WholeMeshSampleRadius;

	if (!StyleSynthesisConfig::mSample_ApplyExtraFiltering) {
		// classical Poisson disk sampling
		SampleSimplePoissonDisk sspd(&mesh);
		if (!sspd.runSampling(sampleRadius)) return false;
		if (!sspd.exportSample(sample)) return false;
		if (!SampleUtil::saveSample(sampleName, sample)) return false;
	} else {
		// extra filtering:
		// -- interior face pruning
		// -- virtual ground
		// -- fail rate checking
		SamplePoissonDisk spd(&mesh);
		if (!spd.runSampling(sampleRadius)) return false;
		if (!spd.exportSample(sample)) return false;
		if (!SampleUtil::saveSample(sampleName, sample)) return false;
	}

	return true;
}

bool PipelineFeatureIO::runSingleModelFeature(string meshName) {

	// names...

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string wholeMeshName = datasetPrefix + "segment/" + meshName + ".ply";

	string sampleName = datasetPrefix + "sample/" + meshName + ".ply";

	string curveFolder = datasetPrefix + "curve/" + meshName + "/";

	string featureFolder = datasetPrefix + "feature/" + meshName + "/";
	if (!FileUtil::makedir(featureFolder)) return false;

	string featureCurvatureName = featureFolder + "curvature.txt";
	string featureCurvatureVName = featureFolder + "curvature.ply";
	string featureSDFName = featureFolder + "SDF.txt";
	string featureSDFVName = featureFolder + "SDF.ply";
	string featureLFDName = featureFolder + "LFD.txt";
	string featureLFDVName = featureFolder + "LFD.ppm";
	string featureSDName = featureFolder + "SD.txt";
	string featureSDVName = featureFolder + "SD.ply";
	string featureCurveName = featureFolder + "curve.txt";
	string featureCurveVName = featureFolder + "curve.ply";

	if (FileUtil::existsfile(featureCurveName)) return true; // early quit

	// data

	TTriangleMesh mesh;
	TSampleSet sample;

	if (!MeshUtil::loadMesh(wholeMeshName, mesh)) return false;
	if (!SampleUtil::loadSample(sampleName, sample)) return false;

	// algorithm

	if (!FileUtil::existsfile(featureCurvatureName)) {
#ifdef OUTPUT_PROGRESS
		cout << "Feature: curvature..." << endl;
#endif
		vector<FeatureSampleCurvature::TCurvature> feature;
		FeatureSampleCurvature fsc(&sample, &feature);
		if (!fsc.calculate()) return false;
		if (!fsc.visualize(featureCurvatureVName)) return false;
		if (!FeatureSampleCurvature::saveFeature(featureCurvatureName, feature)) return false;
	}

	if (!FileUtil::existsfile(featureSDFName)) {
#ifdef OUTPUT_PROGRESS
		cout << "Feature: Shape Diameter Feature..." << endl;
#endif
		vector<double> feature;
		FeatureShapeDiameterFunctions fsdf(&sample, &mesh, &feature);
		if (!fsdf.calculate()) return false;
		if (!fsdf.visualize(featureSDFVName)) return false;
		if (!FeatureAsset::saveFeature(featureSDFName, feature)) return false;
	}

	if (!FileUtil::existsfile(featureLFDName)) {
#ifdef OUTPUT_PROGRESS
		cout << "Feature: Light Field Descriptors..." << endl;
#endif
		vector<double> feature;
		FeatureLightFieldDescriptors flfd(&mesh, &feature);
		if (!flfd.calculate()) return false;
		if (!flfd.visualize(featureLFDVName)) return false;
		if (!FeatureAsset::saveFeature(featureLFDName, feature)) return false;
	}

	if (!FileUtil::existsfile(featureSDName)) {
#ifdef OUTPUT_PROGRESS
		cout << "Feature: Shape Distributions..." << endl;
#endif
		vector<double> feature;
		FeatureShapeDistributions fsd(&sample, &feature);
		if (!fsd.calculate()) return false;
		if (!fsd.visualize(featureSDVName)) return false;
		if (!FeatureAsset::saveFeature(featureSDName, feature)) return false;
	}

	if (!FileUtil::existsfile(featureCurveName)) {
#ifdef OUTPUT_PROGRESS
		cout << "Feature: Curves..." << endl;
#endif
		FeatureCurve fc;
		if (!fc.loadCurves(curveFolder)) return false;
		if (!fc.extractPointClouds(sample.radius)) return false;
		if (!fc.visualize(featureCurveVName)) return false;
		if (!fc.saveFeature(featureCurveName)) return false;
	}

	return true;
}


bool PipelineFeatureIO::runSingleModelPartFeature(string meshName) {

	// names...

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string wholeMeshName = datasetPrefix + "segment/" + meshName + ".ply";

	string segmentName = datasetPrefix + "segment/" + meshName + "-element.txt";

	string featureFolder = datasetPrefix + "feature/" + meshName + "/";
	if (!FileUtil::makedir(featureFolder)) return false;

	string featurePartLFDName = featureFolder + "part-LFD.txt";
	string featurePartSDName = featureFolder + "part-SD.txt";

	if (FileUtil::existsfile(featurePartSDName)) return true; // early quit

	// data

	TTriangleMesh mesh;
	vector<vector<int>> segments;

	vector<TTriangleMesh> partMeshes;
	vector<TSampleSet> partSamples;

	if (!MeshUtil::loadMesh(wholeMeshName, mesh)) return false;
	if (!SegmentUtil::loadSegmentationData(segmentName, segments)) return false;
	int numParts = (int)segments.size();

	partMeshes.resize(numParts);
#pragma omp parallel for
	for (int partID = 0; partID < numParts; partID++) {
		TTriangleMesh &part = partMeshes[partID];
		part.positions = mesh.positions;
		part.normals = mesh.normals;
		part.amount = mesh.amount;
		part.indices.clear();
		for (int faceID : segments[partID]) part.indices.push_back(mesh.indices[faceID]);
		if (!MeshUtil::cleanUp(part)) error("part mesh clean up");
	}

	partSamples.resize(numParts);
#pragma omp parallel for
	for (int partID = 0; partID < numParts; partID++) {
		SampleSimplePoissonDisk sspd(&partMeshes[partID]);
		if (!sspd.runSampling(1024)) error("sampling part mesh");
		if (!sspd.exportSample(partSamples[partID])) error("exporting part sampels");
	}

	// algorithm
	
	if (!FileUtil::existsfile(featurePartLFDName)) {
#ifdef OUTPUT_PROGRESS
		cout << "Part Feature: Light Field Descriptors..." << endl;
#endif
		vector<vector<double>> partFeatures(numParts);
		for (int partID = 0; partID < numParts; partID++) {
			FeatureLightFieldDescriptors flfd(&partMeshes[partID], &partFeatures[partID]);
			if (!flfd.calculate()) return false;			
		}
		if (!FeatureAsset::savePartFeature(featurePartLFDName, partFeatures)) return false;
	}

	if (!FileUtil::existsfile(featurePartSDName)) {
#ifdef OUTPUT_PROGRESS
		cout << "Part Feature: Shape Distributions..." << endl;
#endif
		vector<vector<double>> partFeatures(numParts);
#pragma omp parallel for
		for (int partID = 0; partID < numParts; partID++) {
			FeatureShapeDistributions fsd(&partSamples[partID], &partFeatures[partID]);
			if (!fsd.calculate()) error("part feature shape distributions");
		}
		if (!FeatureAsset::savePartFeature(featurePartSDName, partFeatures)) return false;
	}

	return true;
}

bool PipelineFeatureIO::runSingleModelSaliency(string meshName, int meshID) {

	// names...

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string wholeMeshName = datasetPrefix + "segment/" + meshName + ".ply";

	string segmentName = datasetPrefix + "segment/" + meshName + "-segment.txt";

	string sampleName = datasetPrefix + "sample/" + meshName + ".ply";

	string featureFolder = datasetPrefix + "feature/" + meshName + "/";

	string graphFolder = datasetPrefix + "graph/" + meshName + "/";
	string graphHierarchyName = graphFolder + "graph-hierarchy.txt";
	string graphContextName = graphFolder + "graph-context.txt";

	string saliencyName = featureFolder + "saliency.txt";

	string meshSaliencyName = datasetPrefix + "saliency/" + to_string(meshID + 1) + ".txt"; // 1-indexed
	if (!FileUtil::makedir(meshSaliencyName)) return false;

	if (FileUtil::existsfile(meshSaliencyName)) return true; // early quit

	// algorithm

	//if (!FileUtil::existsfile(saliencyName)) {
	if (true) {

		TTriangleMesh mesh;
		TSampleSet samples;
		vector<vector<vector<int>>> segments;
		FeatureAsset features;
		Eigen::MatrixXd saliency;
		ContextPartGraph graph;
		ContextPartGraphNodeGenerator graphNodeGenerator;

		if (!MeshUtil::loadMesh(wholeMeshName, mesh)) return false;
		if (!SampleUtil::loadSample(sampleName, samples)) return false;
		if (!SegmentGroupApxCvx::loadSegments(segmentName, segments)) return false;
		if (!features.loadAllFeatures(featureFolder)) return false;

		if (!graph.loadGraphHierarchy(graphHierarchyName, graphNodeGenerator, &mesh, &segments)) return false;
		if (!graph.loadGraphContext(graphContextName)) return false;

		FeatureSaliency fs;
		if (!fs.loadData(&mesh, &samples, &segments, &features, &graph)) return false;
		if (!fs.process()) return false;
		if (!fs.output(saliency)) return false;

		if (!DataUtil::saveMatrixBinary(saliencyName, saliency)) return false;
	}

	if (!FileUtil::copyfile(saliencyName, meshSaliencyName)) return false;


	return true;
}

bool PipelineFeatureIO::error(string s) {
	cout << "Error: " << s << endl;
	return false;
}