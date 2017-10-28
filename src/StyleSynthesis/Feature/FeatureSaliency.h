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

#include "Data/StyleSynthesisTypes.h"

using namespace std;

namespace StyleSynthesis {

	class FeatureAsset;
	class ContextPartGraph;
	class ContextPartGraphNode;

	class FeatureSaliency {

		// feature saliency

	public:

		FeatureSaliency();
		~FeatureSaliency();
	
	public:

		bool loadData(
			TTriangleMesh *mesh,
			TSampleSet *samples,
			vector<vector<vector<int>>> *segments,
			FeatureAsset *features,
			ContextPartGraph *graph);

		bool process();
		bool output(Eigen::MatrixXd &saliency); // # of points X # of feature vector dimensions
		
		static bool visualize(string fileName, TSampleSet &samples, Eigen::MatrixXd &saliency);

	private:

		bool computeNodeDescriptor();
		bool computeNodeFeature();
		bool enforceSymmetry();
		bool computePointSaliency();

		inline static bool error(string s) { cout << "Error: " << s << endl; return false; }

	private:

		TTriangleMesh *mpMesh;
		TSampleSet *mpSamples;
		vector<vector<vector<int>>> *mpSegments;
		FeatureAsset *mpFeatures;
		ContextPartGraph *mpGraph;

		vector<ContextPartGraphNode*> mNodes; // graph node : # of nodes used for saliency
		vector<vector<int>> mNodeSamples; // sample ID : # of samples in node : # of nodes
		vector<int> mNodeSampleMapping; // node ID : # of samples
		vector<vector<int>> mNodeSymmetryCliques; // node ID : # of nodes in clique : # of symmetry cliques
		vector<vector<double>> mNodeDescriptors; // descriptor vector : # of nodes
		vector<vector<double>> mNodeFeatures; // feature vector : # of nodes
		vector<vector<double>> mPointSaliency; // feature value : # of features : # of points

	};
}