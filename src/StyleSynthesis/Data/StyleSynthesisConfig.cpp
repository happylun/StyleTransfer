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

#include "StyleSynthesisConfig.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>

using namespace StyleSynthesis;

// default configuration parameters

bool   StyleSynthesisConfig::mPipeline_PipelineMesh            = false;
bool   StyleSynthesisConfig::mPipeline_PipelineCurve           = false;
bool   StyleSynthesisConfig::mPipeline_PipelineSegment         = false;
bool   StyleSynthesisConfig::mPipeline_PipelineGraph           = false;
bool   StyleSynthesisConfig::mPipeline_PipelineFeature         = false;
bool   StyleSynthesisConfig::mPipeline_PipelineSimilarity      = false;
bool   StyleSynthesisConfig::mPipeline_PipelineTrainPart       = false;
bool   StyleSynthesisConfig::mPipeline_PipelineTrainCurve      = false;
bool   StyleSynthesisConfig::mPipeline_PipelineTrainValidation = false;
bool   StyleSynthesisConfig::mPipeline_PipelineTrainLearning   = false;
bool   StyleSynthesisConfig::mPipeline_PipelineMatchPart       = false;
bool   StyleSynthesisConfig::mPipeline_PipelineMatchCurve      = false;
bool   StyleSynthesisConfig::mPipeline_PauseAtFinish           = false;
int    StyleSynthesisConfig::mPipeline_MaximumThreads          = 0;

string StyleSynthesisConfig::mData_DataSetRootFolder = "";
string StyleSynthesisConfig::mData_CustomString1 = "";
string StyleSynthesisConfig::mData_CustomString2 = "";
string StyleSynthesisConfig::mData_CustomString3 = "";
double StyleSynthesisConfig::mData_CustomNumber1 = 0;
double StyleSynthesisConfig::mData_CustomNumber2 = 0;
double StyleSynthesisConfig::mData_CustomNumber3 = 0;
TDList StyleSynthesisConfig::mData_CustomNumberList1 = TDList();
TDList StyleSynthesisConfig::mData_CustomNumberList2 = TDList();
TDList StyleSynthesisConfig::mData_CustomNumberList3 = TDList();

int    StyleSynthesisConfig::mSample_WholeMeshSampleNumber   = 20000;
double StyleSynthesisConfig::mSample_WholeMeshSampleRadius   = 0.02;
double StyleSynthesisConfig::mSample_MinimumSampleRate       = 0.8;
int    StyleSynthesisConfig::mSample_MaximumFailedCount      = 10000;
int    StyleSynthesisConfig::mSample_MaximumCheckedFaceCount = 10000;
bool   StyleSynthesisConfig::mSample_ApplyExtraFiltering     = true;
bool   StyleSynthesisConfig::mSample_VisibilityChecking      = true;
bool   StyleSynthesisConfig::mSample_AddVirtualGround        = false;

TDList StyleSynthesisConfig::mSegmentation_VisibilityThresholds = TDList(1, 0.7, 0.5, 0.3);

double StyleSynthesisConfig::mCurve_RidgeValleyStrength       = 10.0;
double StyleSynthesisConfig::mCurve_RidgeValleyLength         = 0.0;
double StyleSynthesisConfig::mCurve_SamplingRadius            = 0.01;
double StyleSynthesisConfig::mCurve_MaximumChainingAngle      = 45.0;
double StyleSynthesisConfig::mCurve_MinimumSegmentLength      = 0.1;
bool   StyleSynthesisConfig::mCurve_SmoothNormal              = true;
bool   StyleSynthesisConfig::mCurve_SmoothCurvature           = false;
bool   StyleSynthesisConfig::mCurve_SmoothCurvatureDerivative = true;
double StyleSynthesisConfig::mCurve_VisualizationTubeRadius   = 0.003;

double StyleSynthesisConfig::mStyle_ElementDistanceFilteringThreshold  = 0.2;
double StyleSynthesisConfig::mStyle_ReplaceableContributionThreshold   = 0.01;
double StyleSynthesisConfig::mStyle_TransferrableContributionThreshold = 0.02;
bool   StyleSynthesisConfig::mStyle_UseContextualSaliencyFeature       = true;

int    StyleSynthesisConfig::mContext_GraphPropagationIteration           = 5;
TDList StyleSynthesisConfig::mContext_GraphNodePrimitiveAspectWeights     = TDList(1.0, 2.0, 3.0);
int    StyleSynthesisConfig::mContext_FunctionalityLearningIteration      = 1000;
double StyleSynthesisConfig::mContext_FunctionalityLearningStepSize       = 0.1;
double StyleSynthesisConfig::mContext_FunctionalityLearningRegularization = 0.001;
int    StyleSynthesisConfig::mContext_FunctionalityLearningCrossValidationSplits = 10;
double StyleSynthesisConfig::mContext_MatchNodeSimilarityThreshold        = 0.2;
int    StyleSynthesisConfig::mContext_MatchNodeLevels                     = 1;
bool   StyleSynthesisConfig::mContext_SymmetryCheckingByVertex            = true;
double StyleSynthesisConfig::mContext_SymmetryCheckingThreshold           = 2.0;
bool   StyleSynthesisConfig::mContext_SymmetryCheckingPruneByHeight       = true;
bool   StyleSynthesisConfig::mContext_HandleRotationalSymmetry            = true;
bool   StyleSynthesisConfig::mContext_UseLagaDescriptors                  = false;

int    StyleSynthesisConfig::mAssemble_SlotsAlignmentIteration          = 10;
bool   StyleSynthesisConfig::mAssemble_SlotsAlignmentAllowEarlyQuit     = false;
bool   StyleSynthesisConfig::mAssemble_SlotsMatchingAggressiveFiltering = false;
double StyleSynthesisConfig::mAssemble_SlotPointsPruningMultiplier      = 2.0;
double StyleSynthesisConfig::mAssemble_SlotPointsPruningPercentile      = 0.7;
double StyleSynthesisConfig::mAssemble_SlotsMatchingNearestDistance     = 0.2;
double StyleSynthesisConfig::mAssemble_SlotsMatchingFarthestDistance    = 0.6;
double StyleSynthesisConfig::mAssemble_SlotsReusingDistanceThreshold    = 0.5;
double StyleSynthesisConfig::mAssemble_SlotsReusingWeightFactor         = 5.0;
double StyleSynthesisConfig::mAssemble_SlotsAlignmentErrorThreshold     = 0.04;
double StyleSynthesisConfig::mAssemble_SlotsAddingErrorThreshold        = 0.001;
double StyleSynthesisConfig::mAssemble_SlotsSamplingRadius              = 0.005;
bool   StyleSynthesisConfig::mAssemble_AddFloorSlot                     = false;
bool   StyleSynthesisConfig::mAssemble_AddCeilingSlot                   = false;
double StyleSynthesisConfig::mAssemble_PrimitiveAngleThreshold          = 30.0;
double StyleSynthesisConfig::mAssemble_OrientationAngleThreshold        = 10.0;
bool   StyleSynthesisConfig::mAssemble_AllowSphereAlignment             = false;
bool   StyleSynthesisConfig::mAssemble_AllowMinimalAlignment            = false;
bool   StyleSynthesisConfig::mAssemble_AllowGlobalAlignment             = true;
bool   StyleSynthesisConfig::mAssemble_AllowPostAlignment               = false;
bool   StyleSynthesisConfig::mAssemble_AllowRemovingParts               = true;
bool   StyleSynthesisConfig::mAssemble_AllowAddingParts                 = true;
double StyleSynthesisConfig::mAssemble_PostRegularizationEpsilon        = 0.4;
bool   StyleSynthesisConfig::mAssemble_PostRegularizationRounding       = true;
bool   StyleSynthesisConfig::mAssemble_PostAlignCandidateSlots          = false;
int    StyleSynthesisConfig::mAssemble_InitialAlignmentMode             = -1;
int    StyleSynthesisConfig::mAssemble_GlobalAlignmentMode              = -1;
TIList StyleSynthesisConfig::mAssemble_BestGuessAlignmentMode           = TIList(0, 3, 1);
bool   StyleSynthesisConfig::mAssemble_DebugVisualization               = false;
bool   StyleSynthesisConfig::mAssemble_DebugNoBranching                 = false;

string StyleSynthesisConfig::mDeform_Method                        = "None";
double StyleSynthesisConfig::mDeform_TetGenRadiusEdgeRatio         = 1.414;
double StyleSynthesisConfig::mDeform_CurveImpactRadius             = 0.1;
double StyleSynthesisConfig::mDeform_HandleImpactRadius            = 0.01;
double StyleSynthesisConfig::mDeform_MeshSubdivisionRadius         = 0.005;
double StyleSynthesisConfig::mDeform_InterpolationInterval         = 1.0;
double StyleSynthesisConfig::mDeform_MinimumMatchedCurveLength     = 0.2;
double StyleSynthesisConfig::mDeform_MaximumStraightCurveCurvature = 0.001;
bool   StyleSynthesisConfig::mDeform_AlignCurveEndPoints           = true;
int    StyleSynthesisConfig::mDeform_PostSmoothingIterations       = 0;
int    StyleSynthesisConfig::mDeform_Steps                         = 1;


bool StyleSynthesisConfig::loadConfig(string fileName) {

	ifstream cfgFile(fileName);
	if(!cfgFile.is_open()) {
		cout << "Error: cannot load config file " << fileName << endl;
		return false;
	}

	string category;
	while(!cfgFile.eof()) {

		string line;
		getline(cfgFile, line);
		line = trim(line);

		if(line.length() == 0) continue;
		if(line.substr(0,2) == "//") continue;

		if(line.front() == '[' && line.back() == ']') {
			category = trim(line.substr(1, line.length()-2));
			continue;
		}
		
		int pos = (int)line.find_first_of('=');
		string key = trim(line.substr(0, pos));
		string value = trim(line.substr(pos+1));

		if(category == "PIPELINE") {
			if (key == "Pipeline Mesh") {
				mPipeline_PipelineMesh = parseBool(value);
			} else if (key == "Pipeline Curve") {
				mPipeline_PipelineCurve = parseBool(value);
			} else if (key == "Pipeline Segment") {
				mPipeline_PipelineSegment = parseBool(value);
			} else if (key == "Pipeline Graph") {
				mPipeline_PipelineGraph = parseBool(value);
			} else if (key == "Pipeline Feature") {
				mPipeline_PipelineFeature = parseBool(value);
			} else if (key == "Pipeline Similarity") {
				mPipeline_PipelineSimilarity = parseBool(value);
			} else if (key == "Pipeline Train Part") {
				mPipeline_PipelineTrainPart = parseBool(value);
			} else if (key == "Pipeline Train Curve") {
				mPipeline_PipelineTrainCurve = parseBool(value);
			} else if (key == "Pipeline Train Validation") {
				mPipeline_PipelineTrainValidation = parseBool(value);
			} else if (key == "Pipeline Train Learning") {
				mPipeline_PipelineTrainLearning = parseBool(value);
			} else if (key == "Pipeline Match Part") {
				mPipeline_PipelineMatchPart = parseBool(value);
			} else if (key == "Pipeline Match Curve") {
				mPipeline_PipelineMatchCurve = parseBool(value);
			} else if (key == "Pause At Finish") {
				mPipeline_PauseAtFinish = parseBool(value);
			} else if (key == "Maximum Threads") {
				mPipeline_MaximumThreads = parseInt(value);
			} else {
				cout << "Error: unrecognized config => " << line << endl;
				system("pause");
			}

		} else if (category == "DATA") {

			if (key == "Data Set Root Folder") {
				mData_DataSetRootFolder = value + "/";
			} else if (key == "Custom String 1") {
				mData_CustomString1 = value;
			} else if (key == "Custom String 2") {
				mData_CustomString2 = value;
			} else if (key == "Custom String 3") {
				mData_CustomString3 = value;
			} else if (key == "Custom Number 1") {
				mData_CustomNumber1 = parseDouble(value);
			} else if (key == "Custom Number 2") {
				mData_CustomNumber2 = parseDouble(value);
			} else if (key == "Custom Number 3") {
				mData_CustomNumber3 = parseDouble(value);
			} else if (key == "Custom Number List 1") {
				mData_CustomNumberList1 = parseDoubleList(value);
			} else if (key == "Custom Number List 2") {
				mData_CustomNumberList2 = parseDoubleList(value);
			} else if (key == "Custom Number List 3") {
				mData_CustomNumberList3 = parseDoubleList(value);
			} else {
				cout << "Error: unrecognized config => " << line << endl;
				system("pause");
			}

		} else if(category == "SAMPLE") {

			if (key == "Whole Mesh Sample Number") {
				mSample_WholeMeshSampleNumber = parseInt(value);
			} else if (key == "Whole Mesh Sample Radius") {
				mSample_WholeMeshSampleRadius = parseDouble(value);
			} else if( key == "Minimum Success Rate" ) {
				mSample_MinimumSampleRate = parseDouble(value);
			} else if( key == "Maximum Failed Count" ) {
				mSample_MaximumFailedCount = parseInt(value);
			} else if (key == "Maximum Checked Face Count") {
				mSample_MaximumCheckedFaceCount = parseInt(value);
			} else if (key == "Apply Extra Filtering") {
				mSample_ApplyExtraFiltering = parseBool(value);
			} else if (key == "Visibility Checking") {
				mSample_VisibilityChecking = parseBool(value);
			} else if (key == "Add Virtual Ground") {
				mSample_AddVirtualGround = parseBool(value);
			} else {
				cout << "Error: unrecognized config => " << line << endl;
				system("pause");
			}

		} else if (category == "SEGMENTATION") {

			if (key == "Visibility Thresholds") {
					mSegmentation_VisibilityThresholds = parseDoubleList(value);
			} else {
				cout << "Error: unrecognized config => " << line << endl;
				system("pause");
			}

		} else if (category == "CURVE") {

			if (key == "Ridge Valley Strength") {
				mCurve_RidgeValleyStrength = parseDouble(value);
			} else if (key == "Ridge Valley Length") {
				mCurve_RidgeValleyLength = parseDouble(value);
			} else if (key == "Sampling Radius") {
				mCurve_SamplingRadius = parseDouble(value);
			} else if (key == "Maximum Chaining Angle") {
				mCurve_MaximumChainingAngle = parseDouble(value);
			} else if (key == "Minimum Segment Length") {
				mCurve_MinimumSegmentLength = parseDouble(value);
			} else if (key == "Smooth Normal") {
				mCurve_SmoothNormal = parseBool(value);
			} else if (key == "Smooth Curvature") {
				mCurve_SmoothCurvature = parseBool(value);
			} else if (key == "Smooth Curvature Derivative") {
				mCurve_SmoothCurvatureDerivative = parseBool(value);
			} else if (key == "Visualization Tube Radius") {
				mCurve_VisualizationTubeRadius = parseDouble(value);
			} else {
				cout << "Error: unrecognized config => " << line << endl;
				system("pause");
			}

		} else if (category == "STYLE") {

			if (key == "Element Distance Filtering Threshold") {
				mStyle_ElementDistanceFilteringThreshold = parseDouble(value);
			} else if (key == "Replaceable Contribution Threshold") {
				mStyle_ReplaceableContributionThreshold = parseDouble(value);
			} else if (key == "Transferrable Contribution Threshold") {
				mStyle_TransferrableContributionThreshold = parseDouble(value);
			} else if (key == "Use Contextual Saliency Feature") {
				mStyle_UseContextualSaliencyFeature = parseBool(value);
			} else {
				cout << "Error: unrecognized config => " << line << endl;
				system("pause");
			}

		} else if (category == "CONTEXT") {

			if (key == "Graph Propagation Iteration") {
				mContext_GraphPropagationIteration = parseInt(value);
			} else if (key == "Graph Node Primitive Aspect Weights") {
				mContext_GraphNodePrimitiveAspectWeights = parseDoubleList(value);
			} else if (key == "Functionality Learning Iteration") {
				mContext_FunctionalityLearningIteration = parseInt(value);
			} else if (key == "Functionality Learning Step Size") {
				mContext_FunctionalityLearningStepSize = parseDouble(value);
			} else if (key == "Functionality Learning Regularization") {
				mContext_FunctionalityLearningRegularization = parseDouble(value);
			} else if (key == "Functionality Learning Cross Validation Splits") {
				mContext_FunctionalityLearningCrossValidationSplits = parseInt(value);
			} else if (key == "Match Node Similarity Threshold") {
				mContext_MatchNodeSimilarityThreshold = parseDouble(value);
			} else if (key == "Match Node Levels") {
				mContext_MatchNodeLevels = parseInt(value);
			} else if (key == "Symmetry Checking By Vertex") {
				mContext_SymmetryCheckingByVertex = parseBool(value);
			} else if (key == "Symmetry Checking Threshold") {
				mContext_SymmetryCheckingThreshold = parseDouble(value);
			} else if (key == "Symmetry Checking Prune By Height") {
				mContext_SymmetryCheckingPruneByHeight = parseBool(value);
			} else if (key == "Handle Rotational Symmetry") {
				mContext_HandleRotationalSymmetry = parseBool(value);
			} else if (key == "Use Laga Descriptors") {
				mContext_UseLagaDescriptors = parseBool(value);
			} else {
				cout << "Error: unrecognized config => " << line << endl;
				system("pause");
			}

		} else if (category == "ASSEMBLE") {

			if (key == "Slots Alignment Iteration") {
				mAssemble_SlotsAlignmentIteration = parseInt(value);
			} else if (key == "Slots Alignment Allow Early Quit") {
				mAssemble_SlotsAlignmentAllowEarlyQuit = parseBool(value);
			} else if (key == "Slots Matching Aggressive Filtering") {
				mAssemble_SlotsMatchingAggressiveFiltering = parseBool(value);
			} else if (key == "Slot Points Pruning Multiplier") {
				mAssemble_SlotPointsPruningMultiplier = parseDouble(value);
			} else if (key == "Slot Points Pruning Percentile") {
				mAssemble_SlotPointsPruningPercentile = parseDouble(value);
			} else if (key == "Slots Matching Nearest Distance") {
				mAssemble_SlotsMatchingNearestDistance = parseDouble(value);
			} else if (key == "Slots Matching Farthest Distance") {
				mAssemble_SlotsMatchingFarthestDistance = parseDouble(value);
			} else if (key == "Slots Reusing Distance Threshold") {
				mAssemble_SlotsReusingDistanceThreshold = parseDouble(value);
			} else if (key == "Slots Reusing Weight Factor") {
				mAssemble_SlotsReusingWeightFactor = parseDouble(value);
			} else if (key == "Slots Alignment Error Threshold") {
				mAssemble_SlotsAlignmentErrorThreshold = parseDouble(value);
			} else if (key == "Slots Adding Error Threshold") {
				mAssemble_SlotsAddingErrorThreshold = parseDouble(value);
			} else if (key == "Slots Sampling Radius") {
				mAssemble_SlotsSamplingRadius = parseDouble(value);
			} else if (key == "Add Floor Slot") {
				mAssemble_AddFloorSlot = parseBool(value);
			} else if (key == "Add Ceiling Slot") {
				mAssemble_AddCeilingSlot = parseBool(value);
			} else if (key == "Primitive Angle Threshold") {
				mAssemble_PrimitiveAngleThreshold = parseDouble(value);
			} else if (key == "Orientation Angle Threshold") {
				mAssemble_OrientationAngleThreshold = parseDouble(value);
			} else if (key == "Allow Sphere Alignment") {
				mAssemble_AllowSphereAlignment = parseBool(value);
			} else if (key == "Allow Minimal Alignment") {
				mAssemble_AllowMinimalAlignment = parseBool(value);
			} else if (key == "Allow Global Alignment") {
				mAssemble_AllowGlobalAlignment = parseBool(value);
			} else if (key == "Allow Post Alignment") {
				mAssemble_AllowPostAlignment = parseBool(value);
			} else if (key == "Allow Removing Parts") {
				mAssemble_AllowRemovingParts = parseBool(value);
			} else if (key == "Allow Adding Parts") {
				mAssemble_AllowAddingParts = parseBool(value);
			} else if (key == "Post Regularization Epsilon") {
				mAssemble_PostRegularizationEpsilon = parseDouble(value);
			} else if (key == "Post Regularization Rounding") {
				mAssemble_PostRegularizationRounding = parseBool(value);
			} else if (key == "Post Align Candidate Slots") {
				mAssemble_PostAlignCandidateSlots = parseBool(value);
			} else if (key == "Initial Alignment Mode") {
				mAssemble_InitialAlignmentMode = parseInt(value);
			} else if (key == "Global Alignment Mode") {
				mAssemble_GlobalAlignmentMode = parseInt(value);
			} else if (key == "Best Guess Alignment Mode") {
				mAssemble_BestGuessAlignmentMode = parseIntList(value);
			} else if (key == "Debug Visualization") {
				mAssemble_DebugVisualization = parseBool(value);
			} else if (key == "Debug No Branching") {
				mAssemble_DebugNoBranching = parseBool(value);
			} else {
				cout << "Error: unrecognized config => " << line << endl;
				system("pause");
			}

		} else if (category == "DEFORM") {

			if (key == "Method") {
				mDeform_Method = value;
			} else if (key == "TetGen Radius Edge Ratio") {
				mDeform_TetGenRadiusEdgeRatio = parseDouble(value);
			} else if (key == "Curve Impact Radius") {
				mDeform_CurveImpactRadius = parseDouble(value);
			} else if (key == "Handle Impact Radius") {
				mDeform_HandleImpactRadius = parseDouble(value);
			} else if (key == "Mesh Subdivision Radius") {
				mDeform_MeshSubdivisionRadius = parseDouble(value);
			} else if (key == "Interpolation Interval") {
				mDeform_InterpolationInterval = parseDouble(value);
			} else if (key == "Minimum Matched Curve Length") {
				mDeform_MinimumMatchedCurveLength = parseDouble(value);
			} else if (key == "Maximum Straight Curve Curvature") {
				mDeform_MaximumStraightCurveCurvature = parseDouble(value);
			} else if (key == "Align Curve End Points") {
				mDeform_AlignCurveEndPoints = parseBool(value);
			} else if (key == "Post Smoothing Iterations") {
				mDeform_PostSmoothingIterations = parseInt(value);
			} else if (key == "Steps") {
				mDeform_Steps = parseInt(value);
			} else {
				cout << "Error: unrecognized config => " << line << endl;
				system("pause");
			}

		} else {
			cout << "Error: unrecognized config => " << line << endl;
			system("pause");
		}

	}

	cfgFile.close();

	return true;
}

bool StyleSynthesisConfig::saveConfig(string fileName) {

	ofstream cfgFile(fileName);

	cfgFile << endl << "[PIPELINE]" << endl << endl;

	cfgFile << "Pipeline Mesh"             << " = " << (mPipeline_PipelineMesh ? "true" : "false") << endl;
	cfgFile << "Pipeline Curve"            << " = " << (mPipeline_PipelineCurve ? "true" : "false") << endl;
	cfgFile << "Pipeline Segment"          << " = " << (mPipeline_PipelineSegment ? "true" : "false") << endl;
	cfgFile << "Pipeline Graph"            << " = " << (mPipeline_PipelineGraph ? "true" : "false") << endl;
	cfgFile << "Pipeline Feature"          << " = " << (mPipeline_PipelineFeature ? "true" : "false") << endl;
	cfgFile << "Pipeline Similarity"       << " = " << (mPipeline_PipelineSimilarity ? "true" : "false") << endl;
	cfgFile << "Pipeline Train Part"       << " = " << (mPipeline_PipelineTrainPart ? "true" : "false") << endl;
	cfgFile << "Pipeline Train Curve"      << " = " << (mPipeline_PipelineTrainCurve ? "true" : "false") << endl;
	cfgFile << "Pipeline Train Validation" << " = " << (mPipeline_PipelineTrainValidation ? "true" : "false") << endl;
	cfgFile << "Pipeline Train Learning"   << " = " << (mPipeline_PipelineTrainLearning ? "true" : "false") << endl;
	cfgFile << "Pipeline Match Part"       << " = " << (mPipeline_PipelineMatchPart ? "true" : "false") << endl;
	cfgFile << "Pipeline Match Curve"      << " = " << (mPipeline_PipelineMatchCurve ? "true" : "false") << endl;
	cfgFile << "Pause At Finish"           << " = " << (mPipeline_PauseAtFinish ? "true" : "false") << endl;
	cfgFile << "Maximum Threads"           << " = " << mPipeline_MaximumThreads << endl;

	cfgFile << endl << "[DATA]" << endl << endl;

	cfgFile << "Data Set Root Folder" << " = " << mData_DataSetRootFolder << endl;
	cfgFile << "Custom String 1"      << " = " << mData_CustomString1 << endl;
	cfgFile << "Custom String 2"      << " = " << mData_CustomString2 << endl;
	cfgFile << "Custom String 3"      << " = " << mData_CustomString3 << endl;
	cfgFile << "Custom Number 1"      << " = " << mData_CustomNumber1 << endl;
	cfgFile << "Custom Number 2"      << " = " << mData_CustomNumber2 << endl;
	cfgFile << "Custom Number 3"      << " = " << mData_CustomNumber3 << endl;
	cfgFile << "Custom Number List 1" << " = " << mData_CustomNumberList1 << endl;
	cfgFile << "Custom Number List 2" << " = " << mData_CustomNumberList2 << endl;
	cfgFile << "Custom Number List 3" << " = " << mData_CustomNumberList3 << endl;

	cfgFile << endl << "[SAMPLE]" << endl << endl;

	cfgFile << "Whole Mesh Sample Number"   << " = " << mSample_WholeMeshSampleNumber << endl;
	cfgFile << "Whole Mesh Sample Radius"   << " = " << mSample_WholeMeshSampleRadius << endl;
	cfgFile << "Minimum Success Rate"       << " = " << mSample_MinimumSampleRate << endl;
	cfgFile << "Maximum Failed Count"       << " = " << mSample_MaximumFailedCount << endl;
	cfgFile << "Maximum Checked Face Count" << " = " << mSample_MaximumCheckedFaceCount << endl;
	cfgFile << "Apply Extra Filtering"      << " = " << (mSample_ApplyExtraFiltering ? "true" : "false") << endl;
	cfgFile << "Visibility Checking"        << " = " << (mSample_VisibilityChecking ? "true" : "false") << endl;
	cfgFile << "Add Virtual Ground"         << " = " << (mSample_AddVirtualGround ? "true" : "false") << endl;

	cfgFile << endl << "[SEGMENTATION]" << endl << endl;

	cfgFile << "Visibility Thresholds" << " = " << mSegmentation_VisibilityThresholds << endl;

	cfgFile << endl << "[CURVE]" << endl << endl;

	cfgFile << "Ridge Valley Strength"       << " = " << mCurve_RidgeValleyStrength << endl;
	cfgFile << "Ridge Valley Length"         << " = " << mCurve_RidgeValleyLength << endl;
	cfgFile << "Sampling Radius"             << " = " << mCurve_SamplingRadius << endl;
	cfgFile << "Maximum Chaining Angle"      << " = " << mCurve_MaximumChainingAngle << endl;
	cfgFile << "Minimum Segment Length"      << " = " << mCurve_MinimumSegmentLength << endl;
	cfgFile << "Smooth Normal"               << " = " << (mCurve_SmoothNormal ? "true" : "false") << endl;
	cfgFile << "Smooth Curvature"            << " = " << (mCurve_SmoothCurvature ? "true" : "false") << endl;
	cfgFile << "Smooth Curvature Derivative" << " = " << (mCurve_SmoothCurvatureDerivative ? "true" : "false") << endl;
	cfgFile << "Visualization Tube Radius"   << " = " << mCurve_VisualizationTubeRadius << endl;

	cfgFile << endl << "[STYLE]" << endl << endl;

	cfgFile << "Element Distance Filtering Threshold" << " = " << mStyle_ElementDistanceFilteringThreshold << endl;
	cfgFile << "Replaceable Contribution Threshold"   << " = " << mStyle_ReplaceableContributionThreshold << endl;
	cfgFile << "Transferrable Contribution Threshold" << " = " << mStyle_TransferrableContributionThreshold << endl;
	cfgFile << "Use Contextual Saliency Feature"      << " = " << (mStyle_UseContextualSaliencyFeature ? "true" : "false") << endl;

	cfgFile << endl << "[CONTEXT]" << endl << endl;

	cfgFile << "Graph Propagation Iteration"           << " = " << mContext_GraphPropagationIteration << endl;
	cfgFile << "Graph Node Primitive Aspect Weights"   << " = " << mContext_GraphNodePrimitiveAspectWeights << endl;
	cfgFile << "Functionality Learning Iteration"      << " = " << mContext_FunctionalityLearningIteration << endl;
	cfgFile << "Functionality Learning Step Size"      << " = " << mContext_FunctionalityLearningStepSize << endl;
	cfgFile << "Functionality Learning Regularization" << " = " << mContext_FunctionalityLearningRegularization << endl;
	cfgFile << "Functionality Learning Cross Validation Splits" << " = " << mContext_FunctionalityLearningCrossValidationSplits << endl;
	cfgFile << "Match Node Similarity Threshold"       << " = " << mContext_MatchNodeSimilarityThreshold << endl;
	cfgFile << "Match Node Levels"                     << " = " << mContext_MatchNodeLevels << endl;
	cfgFile << "Symmetry Checking By Vertex"           << " = " << (mContext_SymmetryCheckingByVertex ? "true" : "false") << endl;
	cfgFile << "Symmetry Checking Threshold"           << " = " << mContext_SymmetryCheckingThreshold << endl;
	cfgFile << "Symmetry Checking Prune By Height"     << " = " << (mContext_SymmetryCheckingPruneByHeight ? "true" : "false") << endl;
	cfgFile << "Handle Rotational Symmetry"            << " = " << (mContext_HandleRotationalSymmetry ? "true" : "false") << endl;
	cfgFile << "Use Laga Descriptors"                  << " = " << (mContext_UseLagaDescriptors ? "true" : "false") << endl;

	cfgFile << endl << "[ASSEMBLE]" << endl << endl;

	cfgFile << "Slots Alignment Iteration"           << " = " << mAssemble_SlotsAlignmentIteration << endl;
	cfgFile << "Slots Alignment Allow Early Quit"    << " = " << (mAssemble_SlotsAlignmentAllowEarlyQuit ? "true" : "false") << endl;
	cfgFile << "Slots Matching Aggressive Filtering" << " = " << (mAssemble_SlotsMatchingAggressiveFiltering ? "true" : "false") << endl;
	cfgFile << "Slot Points Pruning Multiplier"      << " = " << mAssemble_SlotPointsPruningMultiplier << endl;
	cfgFile << "Slot Points Pruning Percentile"      << " = " << mAssemble_SlotPointsPruningPercentile << endl;
	cfgFile << "Slots Matching Nearest Distance"     << " = " << mAssemble_SlotsMatchingNearestDistance << endl;
	cfgFile << "Slots Matching Farthest Distance"    << " = " << mAssemble_SlotsMatchingFarthestDistance << endl;
	cfgFile << "Slots Reusing Distance Threshold"    << " = " << mAssemble_SlotsReusingDistanceThreshold << endl;
	cfgFile << "Slots Reusing Weight Factor"         << " = " << mAssemble_SlotsReusingWeightFactor << endl;
	cfgFile << "Slots Alignment Error Threshold"     << " = " << mAssemble_SlotsAlignmentErrorThreshold << endl;
	cfgFile << "Slots Adding Error Threshold"        << " = " << mAssemble_SlotsAddingErrorThreshold << endl;
	cfgFile << "Slots Sampling Radius"               << " = " << mAssemble_SlotsSamplingRadius << endl;
	cfgFile << "Add Floor Slot"                      << " = " << (mAssemble_AddFloorSlot ? "true" : "false") << endl;
	cfgFile << "Add Ceiling Slot"                    << " = " << (mAssemble_AddCeilingSlot ? "true" : "false") << endl;
	cfgFile << "Primitive Angle Threshold"           << " = " << mAssemble_PrimitiveAngleThreshold << endl;
	cfgFile << "Orientation Angle Threshold"         << " = " << mAssemble_OrientationAngleThreshold << endl;
	cfgFile << "Allow Sphere Alignment"              << " = " << (mAssemble_AllowSphereAlignment ? "true" : "false") << endl;
	cfgFile << "Allow Minimal Alignment"             << " = " << (mAssemble_AllowMinimalAlignment ? "true" : "false") << endl;
	cfgFile << "Allow Global Alignment"              << " = " << (mAssemble_AllowGlobalAlignment ? "true" : "false") << endl;
	cfgFile << "Allow Post Alignment"                << " = " << (mAssemble_AllowPostAlignment ? "true" : "false") << endl;
	cfgFile << "Allow Removing Parts"                << " = " << (mAssemble_AllowRemovingParts ? "true" : "false") << endl;
	cfgFile << "Allow Adding Parts"                  << " = " << (mAssemble_AllowAddingParts ? "true" : "false") << endl;
	cfgFile << "Post Regularization Epsilon"         << " = " << mAssemble_PostRegularizationEpsilon << endl;
	cfgFile << "Post Regularization Rounding"        << " = " << (mAssemble_PostRegularizationRounding ? "true" : "false") << endl;
	cfgFile << "Post Align Candidate Slots"          << " = " << (mAssemble_PostAlignCandidateSlots ? "true" : "false") << endl;
	cfgFile << "Initial Alignment Mode"              << " = " << mAssemble_InitialAlignmentMode << endl;
	cfgFile << "Global Alignment Mode"               << " = " << mAssemble_GlobalAlignmentMode << endl;
	cfgFile << "Best Guess Alignment Mode"           << " = " << mAssemble_BestGuessAlignmentMode << endl;
	cfgFile << "Debug Visualization"                 << " = " << (mAssemble_DebugVisualization ? "true" : "false") << endl;
	cfgFile << "Debug No Branching"                  << " = " << (mAssemble_DebugNoBranching ? "true" : "false") << endl;

	cfgFile << endl << "[DEFORM]" << endl << endl;

	cfgFile << "Method"                           << " = " << mDeform_Method << endl;
	cfgFile << "TetGen Radius Edge Ratio"         << " = " << mDeform_TetGenRadiusEdgeRatio << endl;
	cfgFile << "Curve Impact Radius"              << " = " << mDeform_CurveImpactRadius << endl;
	cfgFile << "Handle Impact Radius"             << " = " << mDeform_HandleImpactRadius << endl;
	cfgFile << "Mesh Subdivision Radius"          << " = " << mDeform_MeshSubdivisionRadius << endl;
	cfgFile << "Interpolation Interval"           << " = " << mDeform_InterpolationInterval << endl;
	cfgFile << "Minimum Matched Curve Length"     << " = " << mDeform_MinimumMatchedCurveLength << endl;
	cfgFile << "Maximum Straight Curve Curvature" << " = " << mDeform_MaximumStraightCurveCurvature << endl;
	cfgFile << "Align Curve End Points"           << " = " << (mDeform_AlignCurveEndPoints ? "true" : "false") << endl;
	cfgFile << "Post Smoothing Iterations"        << " = " << mDeform_PostSmoothingIterations << endl;
	cfgFile << "Steps"                            << " = " << mDeform_Steps << endl;

	cfgFile.close();

	return true;
}

string StyleSynthesisConfig::trim(string s) {

	// ref: http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring

	s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
	return s;
}

int StyleSynthesisConfig::parseInt(string s) {
	stringstream ss(s);
	int v;
	ss >> v;
	if (ss.fail()) {
		cout << "Error: incorrect integer value: " << s << endl;
		system("pause");
	}
	return v;
}

double StyleSynthesisConfig::parseDouble(string s) {
	stringstream ss(s);
	double v;
	ss >> v;
	if (ss.fail()) {
		cout << "Error: incorrect double value: " << s << endl;
		system("pause");
	}
	return v;
}

bool StyleSynthesisConfig::parseBool(string s) {
	
	string name = s;
	transform(name.begin(), name.end(), name.begin(), ::tolower);
	if(s == "true") {
		return true;
	} else if(s == "false") {
		return false;
	}
	
	cout << "Error: incorrect boolean value: " << s << endl;
	system("pause");
	return false;
}

TIList StyleSynthesisConfig::parseIntList(string s) {

	stringstream ss(s);
	TIList list;
	list.values.clear();
	while(!ss.eof()) {
		int v;
		ss >> v;
		if(!ss.eof() && ss.fail()) {
			cout << "Error: incorrect list \'" << s << "\'" << endl;
			system("pause");
		}
		list.values.push_back(v);
	}
	return list;
}

TDList StyleSynthesisConfig::parseDoubleList(string s) {

	stringstream ss(s);
	TDList list;
	list.values.clear();
	while(!ss.eof()) {
		double v;
		ss >> v;
		if(!ss.eof() && ss.fail()) {
			cout << "Error: incorrect list \'" << s << "\'" << endl;
			system("pause");
		}
		list.values.push_back(v);
	}
	return list;
}
