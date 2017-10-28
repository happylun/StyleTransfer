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

#pragma once

#include <vector>

#include <Eigen/Dense>

#include "Feature/FeatureAsset.h"

#include "Data/StyleSynthesisTypes.h"

using namespace std;

namespace StyleSynthesis {

	class SimilarityData {

		friend class SimilarityDistance;
		friend class SimilarityMetric;

	public:

		SimilarityData();
		~SimilarityData();

	public:

		bool loadData(string sourceName, string targetName);

	protected:

		// mesh

		TTriangleMesh mSourceMesh;
		TTriangleMesh mTargetMesh;

		// sample

		TSampleSet mSourceSamples;
		TSampleSet mTargetSamples;

		SKDTree mSourceSamplesKdTree;
		SKDTreeData mSourceSamplesKdTreeData;

		SKDTree mTargetSamplesKdTree;
		SKDTreeData mTargetSamplesKdTreeData;

		// segment

		vector<vector<int>> mSourceSegments; // face ID : # of faces : # of segments
		vector<vector<int>> mTargetSegments; // face ID : # of faces : # of segments

		// feature

		FeatureAsset mSourceFeatures;
		FeatureAsset mTargetFeatures;

		// saliency

		Eigen::MatrixXd mSourceSaliency; // # of source samples X # of feature dimensions
		Eigen::MatrixXd mTargetSaliency; // # of target samples X # of feature dimensions
	};

}