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

#include "Data/StyleSynthesisTypes.h"

using namespace std;

namespace StyleSynthesis {

	class ContextPartGraph;
	class ContextPartGraphNode;
	class ContextPartGraphNodeGenerator;

	class ContextPartGraphAssembleUtil {

	private:

		ContextPartGraphAssembleUtil();
		~ContextPartGraphAssembleUtil();

		friend class ContextPartGraphAssemble;
		friend class ContextPartGraphAssembleResult;

	public:

		typedef ContextPartGraph TGraph;
		typedef ContextPartGraphNode TNode;
		typedef ContextPartGraphNodeGenerator TNodeGen;

	protected:

		struct TSlot {
			Eigen::Matrix3Xd samples;
			Eigen::Vector3d center;
			vec2i selfIndex; // (node ID, slot ID)
			vec2i adjacentIndex; // (node ID, slot ID)
			bool isVirtualSlot;
			bool isReUsedSlot;
			TSlot() {
				selfIndex = vec2i(-1, -1);
				adjacentIndex = vec2i(-1, -1);
				isVirtualSlot = false;
				isReUsedSlot = false;
			}
		};

		enum TNodeType {
			TNodeType_ACTIVE,
			TNodeType_INACTIVE,
			TNodeType_REPLACED,
			TNodeType_REMOVED,
			TNodeType_MODIFIED
		};

		static bool markNodeTypes(
			TGraph *inGraph,
			vector<int> &inReplaceNodeList, // currently being replaced
			vector<int> &inRemoveNodeList, // currently being removed
			vector<int> &inPreviousNodeList, // previously replaced
			vector<TNodeType> &outTypes);

		static bool extractNodeMesh(
			TNode *inNode,
			TTriangleMesh &outMesh);

		static bool extractMeshSamples(
			TTriangleMesh &inMesh,
			TSampleSet &outSamples);

		static bool extractNodeSlots(
			TGraph *inGraph,
			vector<TSampleSet> &inSamples,
			vector<int> &inWorkingNodes,
			vector<vector<TSlot>> &outSlots);

		static bool computeSlotDistance(
			TSlot &sourceSlot,
			TSlot &targetSlot,
			double &nearestDistance,
			double &farthestDistance);

		static bool alignMatchedParts(
			int alignmentMode,
			TNode *sourceNode,
			TNode *targetNode,
			TSampleSet &sourceSamples,
			TSampleSet &targetSamples,
			Eigen::Affine3d &matchTransform,
			Eigen::Affine3d &fittingTransform);

		static bool alignMatchedSlots(
			Eigen::Matrix3Xd &inSourceFittingPoints,
			Eigen::Matrix3Xd &inSourceXformedPoints,
			vector<TSlot> &inSourceFittingSlots,
			vector<TSlot> &inSourceXformedSlots,
			vector<TSlot> &inTargetSlots,
			vector<TSlot> &outAlignedSlots);

		static bool alignGlobalSlots(
			int inAlignmentMode,
			vector<vector<int>> &inAlignmentGroups,
			vector<int> &inAlignmentPrimitive, // -1 for non-replaced nodes
			vector<Eigen::Vector3d> &inoutAlignmentOrientation, // transformed
			vector<vector<TSlot>> &inoutWorkingSlots,
			vector<Eigen::Affine3d> &outNodeTransformation);

		static bool alignGlobalSlots(
			vector<int> &inWorkingNodes,
			vector<vector<TSlot>> &inoutWorkingSlots,
			vector<Eigen::Affine3d> &outNodeTransformation); // only used for slots post alignment

		static bool alignExemplarSlots(
			vector<vec2i> &replaceMapping,
			vector<Eigen::Affine3d> &exemplarXform,
			vector<Eigen::Affine3d> &candidateXform,
			vector<vector<TSlot>> &exemplarSlots,
			vector<vector<TSlot>> &candidateSlots,
			vector<vector<TSlot>> &alignedSlots); // only used for slots post alignment

		static bool alignSlotPairs(
			int inAlignmentMode,
			int inAlignmentPrimitive,
			Eigen::Vector3d &inAlignmentOrientation,
			vector<TSlot> &inSourceSlots,
			vector<TSlot> &inTargetSlots,
			Eigen::Affine3d &outTransformation,
			double &outError);

		static bool updateVirtualSlots(
			vector<Eigen::Matrix3Xd> &inSourcePoints, // transformed points on source part : # of matching pairs
			vector<vec2i> &inMatchNodePairs,
			vector<vector<TSlot>> &inWorkingSlots);

		static bool checkSlotsAlignment(
			vector<vector<int>> &inAlignmentGroups,
			vector<vector<TSlot>> &inWorkingSlots,
			Eigen::VectorXd &outError);

		static bool regularizeTransformation(
			Eigen::Affine3d &inTransform,
			Eigen::Affine3d &outTransform,
			double eps);

		static bool regularizeNodeTransformation(
			Eigen::Affine3d &inTransform,
			Eigen::Affine3d &outTransform,
			Eigen::Vector3d &inCenter, double eps);

	protected:

		static bool gDebugVisualization;

	private:

		inline static bool error(string s) { cout << "Error: " << s << endl; return false; }
	};
}