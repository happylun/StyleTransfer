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

#include "ContextPartGraphNode.h"
#include "ContextPartGraphNodeGenerator.h"

#include "Data/StyleSynthesisTypes.h"

using namespace std;

namespace StyleSynthesis {

	class ContextPartGraph {

	public:

		ContextPartGraph();
		~ContextPartGraph();

		typedef ContextPartGraphNode TNode;
		typedef ContextPartGraphNodeGenerator TNodeGen;

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

		bool buildGraphHierarchy(
			TNodeGen &nodeGenerator,
			TTriangleMesh *mesh,
			vector<vector<vector<int>>> *segments);
		bool saveGraphHierarchy(string fileName);
		bool loadGraphHierarchy(string fileName,
			TNodeGen &nodeGenerator,
			TTriangleMesh *mesh,
			vector<vector<vector<int>>> *segments);

		bool buildGraphDescriptor(int mode = 0); // call it after *Hierarchy
		bool saveGraphDescriptor(string fileName);
		bool loadGraphDescriptor(string fileName);

		bool buildGraphContext(int mode = 0); // call it after *Hierarchy and *Descriptor
		bool saveGraphContext(string fileName);
		bool loadGraphContext(string fileName);

	private:

		bool extractGraphLevels(vector<vector<TNode*>> &graphLevels);

		inline bool error(string s) { cout << "Error: " << s << endl; return false; }

	protected:

		// actual node data is stored in ContextPartGraphNodeGenerator instance
		// nodes are organized by layer (parent nodes always precede children nodes)
		// node's UID is the order in mAllNodes

		TNode* mRootNode; // empty root node (UID=-1)
		vector<TNode*> mAllNodes;
	};
}