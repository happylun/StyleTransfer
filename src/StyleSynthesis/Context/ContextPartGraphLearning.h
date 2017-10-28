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
#include <fstream>

#include "ContextPartGraph.h"
#include "Context/ContextPartGraphMatch.h"

#include "Data/StyleSynthesisTypes.h"

#include "Library/cppoptlibHelper.h"

using namespace std;

namespace StyleSynthesis {


class ContextPartGraphLearning : public cppoptlib::Problem < double > {

	public:
		ContextPartGraphLearning();
		~ContextPartGraphLearning();

	public:
		bool setTrainingData(const vector <ContextPartGraphMatch*>&, const vector < vector < tuple< int, int, int> > >&);
		bool process();
		bool output(const string&);
		double computeThresholdSimilarity();

	private:
		void projectedGradientDescent();
		void bfgsSolver();
		double value(const Eigen::VectorXd &);
		void gradient(Eigen::MatrixXd&);  // same as Eigen::MatrixXd
		bool updateWeightsInGraph();
		int numErrors();

		inline bool error(string s) { cout << "Error: " << s << endl; return false; }
		inline double logsigmoid(double in);
		inline double sigmoid(double in);
		inline double sgn(double val) { return (0.0 < val) - (val < 0.0); };

		vector <ContextPartGraphMatch*> matched_graphs;
		vector < vector < tuple< int, int, int> > > training_meshes_triplets;
		bool verbose_mode;

	protected:
		Eigen::VectorXd weights;

	};
}