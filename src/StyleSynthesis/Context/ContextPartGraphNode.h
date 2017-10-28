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

#include "Context/ContextPartGraphNodeDescriptors.h"

#include "Data/StyleSynthesisTypes.h"

using namespace std;

namespace StyleSynthesis {

	class ContextPartGraphNode {

	public:

		ContextPartGraphNode();
		~ContextPartGraphNode();

		friend class ContextPartGraph;
		friend class ContextPartGraphMatch;
		friend class ContextPartGraphAssemble;
		friend class ContextPartGraphAssembleUtil;
		friend class ContextPartGraphAssembleResult;
		friend class ContextPartGraphTabuSearchPart;
		friend class ContextPartGraphMatchCurve;
		friend class ContextPartGraphTabuSearchCurve;
		friend class ContextPartGraphTrain;
		friend class ContextPartGraphCrossValidation;
		friend class FeatureSaliency;
		friend class PipelineTrainLearningIO;
		friend class PipelineMatchPartIO;
		friend class PipelineTrainTreeIO;

	public:

		bool clearNode();
		bool computeDescriptor(int mode = 0);
		bool copyMetaData(ContextPartGraphNode *otherNode);
		bool copyNode(ContextPartGraphNode *otherNode);
		bool transformNode(ContextPartGraphNode *otherNode, Eigen::Affine3d &transformation);

		bool saveHierarchy(ostream &fileStream);
		bool loadHierarchy(istream &fileStream, vector<ContextPartGraphNode*> *allNodes, ContextPartGraphNode *rootNode);
		bool saveDescriptor(ostream &fileStream);
		bool loadDescriptor(istream &fileStream);
		bool saveContext(ostream &fileStream);
		bool loadContext(istream &fileStream, vector<ContextPartGraphNode*> *allNodes);

	public:

		static bool detectSymmetry(ContextPartGraphNode *node1, ContextPartGraphNode *node2);
		static bool detectSymmetryRigid(ContextPartGraphNode *node1, ContextPartGraphNode *node2);
		static bool detectSymmetryNonUniform(ContextPartGraphNode *node1, ContextPartGraphNode *node2);
		static bool detectCoCentricity(ContextPartGraphNode *node1, ContextPartGraphNode *node2);
		static bool detectAdjacency(ContextPartGraphNode *node1, ContextPartGraphNode *node2);
		static bool detectContact(ContextPartGraphNode *node1, ContextPartGraphNode *node2);
		static bool detectSupport(ContextPartGraphNode *node1, ContextPartGraphNode *node2);

		static bool computeNodeDistance(
			ContextPartGraphNode *node1, ContextPartGraphNode *node2,
			vector<double> &distance);
		static bool computeEdgeDistance(
			ContextPartGraphNode *node11, ContextPartGraphNode *node12,
			ContextPartGraphNode *node21, ContextPartGraphNode *node22,
			vector<double> &distance);
		static bool computeSimilarity(
			vector<double> &distance,
			vector<double> &sigma,
			vector<double> &sigmaMultipliers,
			vector<double> &weights,
			double &similarity);

	protected:

		// meta data

		int mUID; // unique ID (order within all nodes)

		// context data

		ContextPartGraphNode* mParent;
		vector<ContextPartGraphNode*> mChildren;
		vector<ContextPartGraphNode*> mSymmetry;
		vector<ContextPartGraphNode*> mCocentric;
		vector<ContextPartGraphNode*> mAdjacent;
		vector<ContextPartGraphNode*> mContact; // side contact
		vector<ContextPartGraphNode*> mSupport; // horizontal support

	private:

		// meta data

		TTriangleMesh *mpGraphMesh;
		vector<vector<vector<int>>> *mpGraphSegments;
		int mPartLevelID;
		int mPartSegmentID; // partSegment = graphSegment[levelID][segmentID]

		// descriptor data

		ContextPartGraphNodeDescriptors mNodeDescriptors;
	};
}