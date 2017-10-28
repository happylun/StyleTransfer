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

	class ContextPartGraphNodeDescriptors {

	public:

		ContextPartGraphNodeDescriptors();
		~ContextPartGraphNodeDescriptors();

		friend class ContextPartGraph;
		friend class ContextPartGraphNode;
		friend class ContextPartGraphMatch;
		friend class ContextPartGraphAssemble;
		friend class ContextPartGraphAssembleUtil;
		friend class ContextPartGraphTabuSearchPart;
		friend class ContextPartGraphMatchCurve;
		friend class ContextPartGraphTrain;

	public:

		bool preprocess(TTriangleMesh &mesh, vector<int> &segment);
		bool preprocess(); // only call it when mNodeMesh is established
		bool compute(TTriangleMesh &mesh, int mode = 0);

		bool saveData(ostream &fileStream);
		bool loadData(istream &fileStream);

	private:

		static bool computeMassDistributions(TTriangleMesh &wholeMesh, TTriangleMesh &nodeMesh, vector<double> &histogram);
		static bool computeShapeDistributions(TSampleSet &samples, vector<double> &histogram);
		static const int NUM_MASS_DIST_GRIDS = 4;

		static bool error(string s) { cout << "Error: " << s << endl; return false; }

	protected:

		//////////////// pre-processed data ////////////////

		TTriangleMesh mNodeMesh;
		TSampleSet mNodeSamples;
		Eigen::Matrix3Xd mSampleMatP;
		Eigen::Matrix3Xd mSampleMatN;

		//////////////// descriptor data ////////////////

		// position
		Eigen::Vector3d mMassCenter; // mass center
		double mCenterHeight; // height of mass center
		double mLowHeight; // height of lowest point
		double mHighHeight; // height of highest point
		double mRadialDistance; // horizontal radius of mass center

		// size
		double mSizeEpsilon; // epsilon for size metrics (sampling radius)
		double mMeshArea; // total area of all mesh faces
		Eigen::AlignedBox3d mBoundingBox; // axis-aligned bounding box
		Eigen::Vector3d mAxisAlignedVariance; // axis-aligned variance
		Eigen::Vector3d mPrincipalVariance; // principal variance by PCA

		// orientation
		int mPrimitive; // 0: stick; 1: plane; 2: sphere
		Eigen::Vector3d mMajorOrientation; // major orientation (PCA axis)

		// mass distribution
		vector<double> mMassDistribution; // histogram of 4x4x4 mass distribution

		// shape distributions (only used for comparison with Laga's descriptors)

		vector<double> mShapeDistributions; // histogram of D2 distances
	};
}