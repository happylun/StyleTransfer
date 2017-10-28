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

#include "SimilarityData.h"

#include "Mesh/MeshUtil.h"
#include "Sample/SampleUtil.h"
#include "Segment/SegmentUtil.h"
#include "Data/DataUtil.h"

#include "Data/StyleSynthesisConfig.h"

using namespace std;
using namespace StyleSynthesis;

SimilarityData::SimilarityData() {
}

SimilarityData::~SimilarityData() {
}

bool SimilarityData::loadData(string sourceName, string targetName) {

	string datasetPrefix = StyleSynthesisConfig::mData_DataSetRootFolder;

	string sourceMeshName = datasetPrefix + "segment/" + sourceName + ".ply";
	string targetMeshName = datasetPrefix + "segment/" + targetName + ".ply";

	string sourceSampleName = datasetPrefix + "sample/" + sourceName + ".ply";
	string targetSampleName = datasetPrefix + "sample/" + targetName + ".ply";

	string sourceSegmentName = datasetPrefix + "segment/" + sourceName + "-element.txt";
	string targetSegmentName = datasetPrefix + "segment/" + targetName + "-element.txt";

	string sourceFeatureFolder = datasetPrefix + "feature/" + sourceName + "/";
	string targetFeatureFolder = datasetPrefix + "feature/" + targetName + "/";

	string sourceSaliencyName = sourceFeatureFolder + "saliency.txt";
	string targetSaliencyName = targetFeatureFolder + "saliency.txt";

	if (!MeshUtil::loadMesh(sourceMeshName, mSourceMesh)) return false;
	if (!MeshUtil::loadMesh(targetMeshName, mTargetMesh)) return false;

	if (!SampleUtil::loadSample(sourceSampleName, mSourceSamples)) return false;
	if (!SampleUtil::loadSample(targetSampleName, mTargetSamples)) return false;
	if (!SampleUtil::buildKdTree(mSourceSamples.positions, mSourceSamplesKdTree, mSourceSamplesKdTreeData)) return false;
	if (!SampleUtil::buildKdTree(mTargetSamples.positions, mTargetSamplesKdTree, mTargetSamplesKdTreeData)) return false;

	if (!SegmentUtil::loadSegmentationData(sourceSegmentName, mSourceSegments)) return false;
	if (!SegmentUtil::loadSegmentationData(targetSegmentName, mTargetSegments)) return false;

	if (!mSourceFeatures.loadAllFeatures(sourceFeatureFolder)) return false;
	if (!mTargetFeatures.loadAllFeatures(targetFeatureFolder)) return false;

	if (!DataUtil::loadMatrixBinary(sourceSaliencyName, mSourceSaliency)) return false;
	if (!DataUtil::loadMatrixBinary(targetSaliencyName, mTargetSaliency)) return false;

	return true;
}