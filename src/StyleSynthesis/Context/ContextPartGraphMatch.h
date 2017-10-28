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
#include <string>

#include <Eigen/Dense>

#include "ContextPartGraph.h"

#include "Data/StyleSynthesisTypes.h"

using namespace std;

namespace StyleSynthesis {

	class ContextPartGraphMatch {

	public:

		ContextPartGraphMatch();
		~ContextPartGraphMatch();

		friend class ContextPartGraphTrain;
		friend class ContextPartGraphLearning;
		friend class ContextPartGraphLearningCurve;
		friend class ContextPartGraphTabuSearchPart;
		friend class ContextPartGraphCrossValidation;
		friend class DebugContextGraphIO;
		friend class DebugAssembleIO;
		friend class PipelineMatchPartIO;
		friend class PipelineTrainPartIO;
		friend class PipelineTrainValidationIO;
		friend class PipelineTrainLearningIO;

		typedef ContextPartGraphNode TNode;
		typedef ContextPartGraph TGraph;

	public:

		bool loadWeights(string weightsFolder);
		bool loadGraph(TGraph &source, TGraph &target);
		bool loadMatchingMode(int mode);
		bool process();
		bool visualize(string fileName);

		bool exportSimilarityMatrix(Eigen::MatrixXd &outMatrix);
		bool exportGraphSimilarity(double &similarity);
		bool exportMatchings(vector<vector<vec2i>> &matchings);
		bool exportMatchingsSimilarity(vector<double> &matchingsSimilarity);
		bool exportMatchingMode(int &mode);

		static bool saveMatchings(string fileName, vector<vector<vec2i>> &matchings);
		static bool loadMatchings(string fileName, vector<vector<vec2i>> &matchings);

	protected:
		
		bool computeNodeDistance();
		bool computeEdgeDistance();
		bool computeNodeSigma();
		bool computeEdgeSigma();
		bool evaluateNodeSimilarity();
		bool evaluateEdgeSimilarity();
		bool computeGraphSimilarity();

		bool normalizeGraphSimilarity(Eigen::MatrixXd &sourceSimilarity, Eigen::MatrixXd &targetSimilarity);

		bool matchGraphNodes();
		bool exportDerivatives(vector<Eigen::MatrixXd>&, vector<Eigen::MatrixXd>&, vector<Eigen::MatrixXd>&, vector<Eigen::MatrixXd>&, vector<Eigen::MatrixXd>&); // V addition

		bool computeEdgeSimilarity(
			int sourceID, int targetID,
			int sourceNeighborID, int targetNeighborID,
			double &similarity);

		bool sliceOutSimilarityMatrix(
			int sourceID, int targetID,
			Eigen::MatrixXd &similarityMatrix);

		inline static bool error(string s) { cout << "Error: " << s << endl; return false; }

	protected:

		TGraph *mpSourceGraph;
		TGraph *mpTargetGraph;

		int mNumSourceNodes;
		int mNumTargetNodes;
		int mNumIterations;

		vector<double> mNodeWeights;
		vector<double> mEdgeWeights;
		vector<double> mPathWeights;

		vector<vector<double>> mRawNodeDistance; // node distance vector : # of node pairs
		vector<vector<double>> mRawEdgeDistance; // edge distance vector : # of edge quads
		vector<double> mNodeSigma; // sigma : dimension of node distance vector
		vector<double> mEdgeSigma; // sigma : dimension of edge distance vector
		vector<double> mNodeSigmaMultipliers; // learned sigma multipliers for node distances
		vector<double> mEdgeSigmaMultipliers; // learned sigma multipliers for edge distances

		int mMatchingMode;
		
		vector<vec4i> mAllEdgeQuads; // (source ID, target ID, source NB ID, target NB ID) : # of edge quads
		vector<vector<int>> mQuadsForPairs; // edge quad ID : # of edge quads for this node pair : # of node pairs

		typedef Eigen::Matrix<long long, -1, -1, 0, -1, -1> EigenMatrixXll;
		
		Eigen::MatrixXd mNodeSimilarityMatrix; // # of source nodes X # of target nodes
		Eigen::VectorXd mEdgeSimilarityVector; // # of edge quads
		vector<EigenMatrixXll> mPathCountMatrix; // // # of path lengths : # of source nodes X # of target nodes
		vector<Eigen::MatrixXd> mGraphSimilarityMatrixPerPathLength; // # of path lengths : # of source nodes X # of target nodes
		Eigen::MatrixXd mGraphSimilarityMatrix; // node similarity : # of source nodes X # of target nodes

		vector<vector<vec2i>> mMatchings; // (source node ID, target node ID) : # of matched pairs in group : # of matching groups
		vector<double> mMatchingsSimilarity; // functionality similarity : # of matching groups
	};
}