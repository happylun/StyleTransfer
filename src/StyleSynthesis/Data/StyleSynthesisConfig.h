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

#include <string>
#include <vector>
#include <ostream>

using namespace std;

namespace StyleSynthesis {

	template <class T>
	class TAbstractList {
	public:
		vector<T> values;
		friend class StyleSynthesisConfig;
		friend ostream& operator<< (ostream& stream, const TAbstractList<T>& list) {
			for(auto &value : list.values) stream << value << " ";
			return stream;
		}

		// NOTE: 3 ways for initialization
		//     1) TDList a = TDList(1.0, 2.0, 3.0);
		//     2) TDList a (1.0, 2.0, 3.0);
		//     3) TDList a {1.0, 2.0, 3.0};
		// ref: http://stackoverflow.com/questions/28866559/writing-variadic-template-constructor
		TAbstractList() = default; // terminate recursion
		template <class T1, class... T2> TAbstractList(T1 n, T2... rest)
			: TAbstractList(rest...) { values.insert(values.begin(), n); } // recursive variadic template initialization
		TAbstractList(initializer_list<T> il) : values(il) {} // std list initializer
	};

	typedef TAbstractList<int>    TIList;
	typedef TAbstractList<double> TDList;

	class StyleSynthesisConfig {

	private:

		// make it non-instantiable
		StyleSynthesisConfig() {}
		~StyleSynthesisConfig() {}

	public:

		static bool loadConfig(string fileName);
		static bool saveConfig(string fileName);

	public:

		static bool   mPipeline_PipelineMesh;
		static bool   mPipeline_PipelineCurve;
		static bool   mPipeline_PipelineSegment;
		static bool   mPipeline_PipelineGraph;
		static bool   mPipeline_PipelineFeature;
		static bool   mPipeline_PipelineSimilarity;
		static bool   mPipeline_PipelineTrainPart;
		static bool   mPipeline_PipelineTrainCurve;
		static bool   mPipeline_PipelineTrainValidation;
		static bool   mPipeline_PipelineTrainLearning;
		static bool   mPipeline_PipelineMatchPart;
		static bool   mPipeline_PipelineMatchCurve;
		static bool   mPipeline_PauseAtFinish;
		static int    mPipeline_MaximumThreads;

		static string mData_DataSetRootFolder;
		static string mData_CustomString1;
		static string mData_CustomString2;
		static string mData_CustomString3;
		static double mData_CustomNumber1;
		static double mData_CustomNumber2;
		static double mData_CustomNumber3;
		static TDList mData_CustomNumberList1;
		static TDList mData_CustomNumberList2;
		static TDList mData_CustomNumberList3;

		static int    mSample_WholeMeshSampleNumber;
		static double mSample_WholeMeshSampleRadius;
		static double mSample_MinimumSampleRate;
		static int    mSample_MaximumFailedCount;
		static int    mSample_MaximumCheckedFaceCount;
		static bool   mSample_ApplyExtraFiltering;
		static bool   mSample_VisibilityChecking;
		static bool   mSample_AddVirtualGround;

		static TDList mSegmentation_VisibilityThresholds;

		static double mCurve_RidgeValleyStrength;
		static double mCurve_RidgeValleyLength;
		static double mCurve_SamplingRadius;
		static double mCurve_MaximumChainingAngle;
		static double mCurve_MinimumSegmentLength;
		static bool   mCurve_SmoothNormal;
		static bool   mCurve_SmoothCurvature;
		static bool   mCurve_SmoothCurvatureDerivative;
		static double mCurve_VisualizationTubeRadius;

		static double mStyle_ElementDistanceFilteringThreshold;
		static double mStyle_ReplaceableContributionThreshold;
		static double mStyle_TransferrableContributionThreshold;
		static bool   mStyle_UseContextualSaliencyFeature;

		static int    mContext_GraphPropagationIteration;
		static TDList mContext_GraphNodePrimitiveAspectWeights;
		static int    mContext_FunctionalityLearningIteration;
		static double mContext_FunctionalityLearningStepSize;
		static double mContext_FunctionalityLearningRegularization;
		static int    mContext_FunctionalityLearningCrossValidationSplits;
		static double mContext_MatchNodeSimilarityThreshold;
		static int    mContext_MatchNodeLevels;
		static bool   mContext_SymmetryCheckingByVertex;
		static double mContext_SymmetryCheckingThreshold;
		static bool   mContext_SymmetryCheckingPruneByHeight;
		static bool   mContext_HandleRotationalSymmetry;
		static bool   mContext_UseLagaDescriptors;

		static int    mAssemble_SlotsAlignmentIteration;
		static bool   mAssemble_SlotsAlignmentAllowEarlyQuit;
		static bool   mAssemble_SlotsMatchingAggressiveFiltering;
		static double mAssemble_SlotPointsPruningMultiplier;
		static double mAssemble_SlotPointsPruningPercentile;
		static double mAssemble_SlotsMatchingNearestDistance;
		static double mAssemble_SlotsMatchingFarthestDistance;
		static double mAssemble_SlotsReusingDistanceThreshold;
		static double mAssemble_SlotsReusingWeightFactor;
		static double mAssemble_SlotsAlignmentErrorThreshold;
		static double mAssemble_SlotsAddingErrorThreshold;
		static double mAssemble_SlotsSamplingRadius;
		static bool   mAssemble_AddFloorSlot;
		static bool   mAssemble_AddCeilingSlot;
		static double mAssemble_PrimitiveAngleThreshold;
		static double mAssemble_OrientationAngleThreshold;
		static bool   mAssemble_AllowSphereAlignment;
		static bool   mAssemble_AllowMinimalAlignment;
		static bool   mAssemble_AllowGlobalAlignment;
		static bool   mAssemble_AllowPostAlignment;
		static bool   mAssemble_AllowRemovingParts;
		static bool   mAssemble_AllowAddingParts;
		static double mAssemble_PostRegularizationEpsilon;
		static bool   mAssemble_PostRegularizationRounding;
		static bool   mAssemble_PostAlignCandidateSlots;
		static int    mAssemble_InitialAlignmentMode;
		static int    mAssemble_GlobalAlignmentMode;
		static TIList mAssemble_BestGuessAlignmentMode;
		static bool   mAssemble_DebugVisualization;
		static bool   mAssemble_DebugNoBranching;

		static string mDeform_Method;
		static double mDeform_TetGenRadiusEdgeRatio;
		static double mDeform_CurveImpactRadius;
		static double mDeform_HandleImpactRadius;
		static double mDeform_MeshSubdivisionRadius;
		static double mDeform_InterpolationInterval;
		static double mDeform_MinimumMatchedCurveLength;
		static double mDeform_MaximumStraightCurveCurvature;
		static bool   mDeform_AlignCurveEndPoints;
		static int    mDeform_PostSmoothingIterations;
		static int    mDeform_Steps;

	private:

		static string trim(string s);

		static int parseInt(string s);
		static double parseDouble(string s);
		static bool parseBool(string s);
		static TIList parseIntList(string s);
		static TDList parseDoubleList(string s);

	};

}