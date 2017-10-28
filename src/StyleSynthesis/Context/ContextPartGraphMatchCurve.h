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

#include <Eigen/Dense>

#include "Data/StyleSynthesisTypes.h"

namespace StyleSynthesis {

	class ContextPartGraph;
	class ContextPartGraphNode;
	class ContextPartGraphAssembleResult;

	class ContextPartGraphMatchCurve {

	public:

		ContextPartGraphMatchCurve();
		~ContextPartGraphMatchCurve();

	public:

		bool loadWeights(string weightsFolder);
		bool loadGraph(ContextPartGraph *srcGraph, ContextPartGraph *tgtGraph, ContextPartGraphAssembleResult *solution, Eigen::MatrixXd &similarity);
		bool loadCurve(string srcCurveFolder, string tgtCurveFolder, vector<vec3> &viewPoints);
		bool loadCandidates(vector<vec2i> &curvePairs);

		bool process();
		bool initialize();
		bool gatherCandidates();
		bool computeDistance(bool doPruning = true);
		bool computeSigma();
		bool computeScore();
		bool computeScoreDerivative();
		bool extractPairs();

		bool visualize(string fileName);
		bool output(
			vector<vector<vector<vec3>>> &srcCurves,
			vector<vector<vector<vec3>>> &tgtCurves,
			vector<vector<int>> &srcNodes,
			vector<vector<int>> &tgtNodes);

		bool loadCurveDistances(string fileName);
		bool saveCurveDistances(string fileName);

	private:

		bool initContext();
		bool initCurves();
		bool findCliques();

		bool findCurveNodes(
			vector<vector<vec3>> &inCurves,
			ContextPartGraph &inGraph,
			vector<vector<int>> &outNodes); // node UID : # of nodes containing curve : # of curves

		bool detectSymmetry(
			Eigen::Matrix3Xd &curveMatSP,
			Eigen::Matrix3Xd &curveMatSN,
			Eigen::Matrix3Xd &curveMatTP,
			Eigen::Matrix3Xd &curveMatTN,
			bool &isSymmetric);

		bool computeCurveDistance(
			int sourceCurveID,
			int targetCurveID,
			vector<double> &outDistances,
			int &outSourceNodeID,
			int &outTargetNodeID,
			double &outNodeScore,
			bool doPruning);

		bool computeMatchingDistances(
			int sourceCurveID,
			int targetCurveID,
			vector<double> &outDistances);

		bool clampCurveToNode(
			vector<vec3> &inCurve,
			ContextPartGraphNode *inNode,
			vector<vec3> &outCurve);

		static inline bool error(string info) { cout << "\nError: " << info << endl; return false; }

	public:

		vector<double> mCurveWeights;
		vector<double> mCurveSigma;
		vector<double> mCurveSigmaMultipliers;

	private:

		ContextPartGraph *mpSourceGraph;
		ContextPartGraph *mpTargetGraph;
		ContextPartGraphAssembleResult *mpSolution;
		Eigen::MatrixXd mNodeSimilarity;

		vector<vec3> mCurveViewPoints;

		vector<vector<vec3>> mSourceAllCurves; // curve point : # of points : # of curves (mixed types)
		vector<vector<vec3>> mTargetAllCurves;

		vector<Eigen::Matrix3Xd> mCurveMatSP; // curve sample points position : # of curves (mixed types)
		vector<Eigen::Matrix3Xd> mCurveMatSN; // curve sample points direction : # of curves (mixed types)
		vector<Eigen::Matrix3Xd> mCurveMatTP;
		vector<Eigen::Matrix3Xd> mCurveMatTN;

		vector<vector<int>> mSourceAllCurvesNodes; // node ID : # of nodes containing curve : # of curves
		vector<vector<int>> mTargetAllCurvesNodes;

		vector<bool> mSourceAllCurvesValidFlags; // whether is valid for matching : # of curves
		vector<bool> mTargetAllCurvesValidFlags;
		vector<bool> mSourceAllCurvesStraightFlags; // whether is straight line : # of curves
		vector<bool> mTargetAllCurvesStraightFlags;

		vector<int> mSourceAllCurvesCliqueID; // clique ID : # of curves
		vector<int> mTargetAllCurvesCliqueID;
		vector<vector<int>> mSourceCliques; // curve ID : # of curves in clique : # of cliques
		vector<vector<int>> mTargetCliques;

		vector<Eigen::Vector3d> mSourceAllCurveRelPos; // relative position to mesh : # of curves
		vector<Eigen::Vector3d> mTargetAllCurveRelPos;

		// type: 0 - ridge/valley; 1 - boundary; 2+ - contours
		vector<vec2i> mSourceCurveIndex; // (start index, end index) : # of curve types
		vector<vec2i> mTargetCurveIndex;

		vector<vec2i> mAllPairIndices; // (source curve ID, target curve ID) : # of all valid matching pairs
		vector<double> mAllPairViewDistances; // view distance : # of all valid matching pairs

	protected:

		vector<vector<double>> mAllPairDistances; // distance entry : # of distance entries (empty if pruned) : # of all valid matching pairs
		vector<vec2i> mAllPairNodes; // (source node ID, target node ID) : # of all valid matching pairs
		vector<double> mAllPairNodeScores; // context score : # of all valid matching pairs

	public:

		Eigen::VectorXd mAllPairScores; // curve score : # of all valid matching pairs
		Eigen::MatrixXd mAllPairScoreDerivatives; // derivative wrt weights : # of all valid matching pairs X # of weights dimensions

		vector<vector<vector<vec3>>> mMatchingsSourceCurves; // curve point : # of points on source curve : # of matched curve pairs : # of matching groups
		vector<vector<vector<vec3>>> mMatchingsTargetCurves; // curve point : # of points on target curve : # of matched curve pairs : # of matching groups
		vector<vector<int>> mMatchingsSourceNodes; // ID of node associated with source curve : # of matched curve pairs : # of matching groups
		vector<vector<int>> mMatchingsTargetNodes; // ID of node associated with target curve : # of matched curve pairs : # of matching groups

	};
}