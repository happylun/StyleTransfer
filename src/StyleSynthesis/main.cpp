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

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <string>
#include <omp.h>

#include "Eigen/Core"

#include "Data/StyleSynthesisConfig.h"

#include "IO/Pipeline/PipelineMeshIO.h"
#include "IO/Pipeline/PipelineCurveIO.h"
#include "IO/Pipeline/PipelineSegmentIO.h"
#include "IO/Pipeline/PipelineFeatureIO.h"
#include "IO/Pipeline/PipelineGraphIO.h"
#include "IO/Pipeline/PipelineSimilarityIO.h"
#include "IO/Pipeline/PipelineTrainPartIO.h"
#include "IO/Pipeline/PipelineTrainCurveIO.h"
#include "IO/Pipeline/PipelineTrainValidationIO.h"
#include "IO/Pipeline/PipelineTrainLearningIO.h"
#include "IO/Pipeline/PipelineMatchPartIO.h"
#include "IO/Pipeline/PipelineMatchCurveIO.h"

#include "Utility/Timer.h"

using namespace std;
using namespace StyleSynthesis;

void error(int e) {

	cout << "Error in main routine: " << e << endl;
	cin.get();
	exit(e);
}

int main(int argc, char** argv) {

	//return StyleSynthesisConfig::saveConfig("Debug/a.data/new-params.cfg");

	TTimer startTime = Timer::tic();

	// init
	Eigen::initParallel();
	//Eigen::setNbThreads(1);
	srand((unsigned int)time(0));
	
	// load all configs
	for(int cfgID=1; cfgID<argc; cfgID++) {
		if (!StyleSynthesisConfig::loadConfig(argv[cfgID])) error(-11);
	}

#ifdef _OPENMP
	// limit number of threads
	if (StyleSynthesisConfig::mPipeline_MaximumThreads) {
		omp_set_num_threads(StyleSynthesisConfig::mPipeline_MaximumThreads);
	} else {
		StyleSynthesisConfig::mPipeline_MaximumThreads = omp_get_max_threads();
	}
#endif

	// pipeline
	if (StyleSynthesisConfig::mPipeline_PipelineMesh            && !PipelineMeshIO::process()) error(-41);
	if (StyleSynthesisConfig::mPipeline_PipelineCurve           && !PipelineCurveIO::process()) error(-42);
	if (StyleSynthesisConfig::mPipeline_PipelineSegment         && !PipelineSegmentIO::process()) error(-43);
	if (StyleSynthesisConfig::mPipeline_PipelineGraph           && !PipelineGraphIO::process()) error(-44);
	if (StyleSynthesisConfig::mPipeline_PipelineFeature         && !PipelineFeatureIO::process()) error(-45);
	if (StyleSynthesisConfig::mPipeline_PipelineSimilarity      && !PipelineSimilarityIO::process()) error(-46);
	if (StyleSynthesisConfig::mPipeline_PipelineTrainPart       && !PipelineTrainPartIO::process()) error(-47);
	if (StyleSynthesisConfig::mPipeline_PipelineTrainCurve      && !PipelineTrainCurveIO::process()) error(-48);
	if (StyleSynthesisConfig::mPipeline_PipelineTrainValidation && !PipelineTrainValidationIO::process()) error(-49);
	if (StyleSynthesisConfig::mPipeline_PipelineTrainLearning   && !PipelineTrainLearningIO::process()) error(-410);
	if (StyleSynthesisConfig::mPipeline_PipelineMatchPart       && !PipelineMatchPartIO::process()) error(-411);
	if (StyleSynthesisConfig::mPipeline_PipelineMatchCurve      && !PipelineMatchCurveIO::process()) error(-412);
	
	cout << "Total time: " << Timer::toString(Timer::toc(startTime)) << endl;
	if (StyleSynthesisConfig::mPipeline_PauseAtFinish) system("pause");

	return 0;
}
