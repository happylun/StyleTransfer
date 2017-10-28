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

#include "ContextPartGraphAssembleUtil.h"

#include "Data/StyleSynthesisTypes.h"

using namespace std;

namespace StyleSynthesis {

	class ContextPartGraph;
	class ContextPartGraphNodeGenerator;
	class ContextPartGraphAssembleResult;

	class ContextPartGraphAssemble {

	public:

		ContextPartGraphAssemble();
		~ContextPartGraphAssemble();

		friend class ContextPartGraphAssembleResult;

	private:

		typedef ContextPartGraph TGraph;
		typedef ContextPartGraphNode TNode;
		typedef ContextPartGraphNodeGenerator TNodeGen;
		typedef ContextPartGraphAssembleUtil TUtil;
		typedef ContextPartGraphAssembleResult TSolution;
		typedef TUtil::TSlot TSlot;
		typedef TUtil::TNodeType TNodeType;

	public:

		bool loadGraph(TNodeGen *nodeGen, TGraph *source, TGraph *target);
		bool loadOldSolution(TSolution *initSolution, TSolution *prevSolution);
		bool loadMatching(vector<vec2i> &matching);
		bool loadRemoving(vector<int> &removing);
		bool loadAdding(vector<int> &adding);
		bool process(bool &outSuccessFlag);
		bool postProcess(bool &outSuccessFlag);
		bool visualize(string folderName);

	private:

		bool initData();
		bool findSlots();
		bool transformParts();
		bool matchSlots();
		bool reuseSlots();
		bool alignSlots(double &outError);
		bool alignAddings(double &outError);
		bool postAlignSlots(double &outError);

		bool visualizeMatchings(string fileName);
		bool visualizeNodes(string fileName);
		bool visualizeSlots(
			string fileName,
			vector<int> &workingNodes,
			vector<vector<TSlot>> &nodeSlots);

	private:

		inline static bool error(string s) { cout << "Error: " << s << endl; return false; }

	protected:

		TNodeGen *mpNodeGenerator;
		TGraph *mpSourceGraph;
		TGraph *mpTargetGraph;
		TSolution *mpPreviousSolution;
		TSolution *mpInitialSolution;
		vector<vec2i> mMatchPairs; // (source node ID, target node ID)
		vector<int> mRemovingNodes; // target node ID
		vector<int> mAddingNodes; // source node ID

		vector<TNodeType> mSourceNodeTypes; // node type : # of nodes on source shape
		vector<TNodeType> mTargetNodeTypes; // node type : # of nodes on target shape
		vector<int> mSourceWorkingNodes; // node ID : # of active/replaced nodes on source shape
		vector<int> mTargetWorkingNodes; // node ID : # of active/replaced nodes on target shape
		vector<int> mSourceAddingWorkingNodes; // node ID : # of nodes to be added to result and their neighbors which are previously replacement nodes

		// only valid for working nodes
		vector<TSampleSet> mSourceNodeSamples; // sample points on node : # of source nodes
		vector<TSampleSet> mTargetNodeSamples; // sample points on node : # of target nodes

		// only valid for working nodes
		vector<vector<TSlot>> mSourceNodeSlots; // slot data : # of slots on source node : # of source nodes
		vector<vector<TSlot>> mTargetNodeSlots; // slot data : # of slots on target node : # of target nodes
		vector<vector<TSlot>> mAlignedNodeSlots; // slot data : # of slots on "assemble" node : # of target nodes
		vector<vector<TSlot>> mAddingNodeSlots; // slot data : # of slots on source node : # of source nodes

		int mAlignmentMode;
		int mNumAlignmentMode; // depend on different type of primitives
		bool mUseBestGuess;

		vector<Eigen::Affine3d> mFittingTransformations; // transformation applied on source part to fit target part (only used for matching slots) : # of matches
		vector<Eigen::Affine3d> mMatchTransformations; // transformation should be applied to source part : # of matches

		// only valid for working nodes
		vector<Eigen::Affine3d> mNodeTransformations; // transformation should be applied to target part : # of target nodes

		vector<int> mAddingNodeIndices; // node ID : # of added source nodes (may be duplicated)
		vector<Eigen::Affine3d> mAddingTransformations; // transformation should be applied to source node : # of added source nodes
	};
}