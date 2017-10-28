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

#include "ContextPartGraphLearning.h"

#include "Data/StyleSynthesisConfig.h"

#include "Utility/Timer.h"

namespace StyleSynthesis {

	ContextPartGraphLearning::ContextPartGraphLearning() {
		verbose_mode = true;
	}

	ContextPartGraphLearning::~ContextPartGraphLearning() {
	}

	bool ContextPartGraphLearning::setTrainingData(const vector <ContextPartGraphMatch*>& _matched_graphs, const vector < vector < tuple< int, int, int> > >& _training_meshes_triplets)
	{
		if (_matched_graphs.empty())
			return error("No training data were found!");
		matched_graphs = _matched_graphs;
		training_meshes_triplets = _training_meshes_triplets;
		return true;
	}

	bool ContextPartGraphLearning::output(const string& learning_folder)
	{
		string functionality_learning_filename = learning_folder + "/weights.txt";
		ofstream functionality_learning_file(functionality_learning_filename);
		for (int w = 0; w < weights.size(); w++)
			functionality_learning_file << weights[w] << std::endl;
		if (!functionality_learning_file.good())
			return error("Cannot write output weights to " + functionality_learning_filename);
		functionality_learning_file.close();

		return true;
	}


	bool ContextPartGraphLearning::process()
	{
		int dimNodeWeights = (int)matched_graphs[0]->mNodeWeights.size();
		int dimEdgeWeights = (int)matched_graphs[0]->mEdgeWeights.size();
		int dimPathWeights = (int)matched_graphs[0]->mPathWeights.size();
		double initialNodeWeights = 1.0 / (double)dimNodeWeights;
		double initialEdgeWeights = 1.0 / (double)dimEdgeWeights;
		double initialPathWeights = 1.0 / (double)dimPathWeights;
		weights.resize(2*dimNodeWeights + 2*dimEdgeWeights + dimPathWeights);
		weights << Eigen::VectorXd::Ones(dimNodeWeights) * initialNodeWeights,
				   Eigen::VectorXd::Ones(dimNodeWeights),  // sigma multipliers
				   Eigen::VectorXd::Ones(dimEdgeWeights) * initialEdgeWeights,
				   Eigen::VectorXd::Ones(dimEdgeWeights),  // sigma multipliers
			       Eigen::VectorXd::Ones(dimPathWeights) * initialPathWeights;

		//projectedGradientDescent();
		bfgsSolver();

		return true;
	}

	void ContextPartGraphLearning::bfgsSolver()
	{
		verbose_mode = false;
		Eigen::VectorXd zero_weights = 0 * weights.array() + 1e-7; 
		setLowerBound(zero_weights);
		cppoptlib::LbfgsbSolver<double> solver;
		solver.settings_.maxIter = StyleSynthesisConfig::mContext_FunctionalityLearningIteration; // many fewer iterations needed here
		solver.settings_.rate = StyleSynthesisConfig::mContext_FunctionalityLearningStepSize;
		solver.minimize(*this, weights);

		std::cout << "Solved. loss:" << value(weights) << std::endl;
		std::cout << "#errors:" << numErrors() << std::endl;
		std::cout << "weights =\n " << weights << std::endl;
	}


	void ContextPartGraphLearning::projectedGradientDescent()
	{
		verbose_mode = false;
		Eigen::VectorXd zero_weights = weights * 0; 
		double objective = value(weights);
		std::cout << "At the beginning - objective: " << objective << ", #errors: " << numErrors() << std::endl;
		double old_objective = objective;

		for (int iter = 1; iter <= StyleSynthesisConfig::mContext_FunctionalityLearningIteration; iter++)
		{
			//////// for debugging derivatives			
			//Eigen::MatrixXd _grad;
			//gradient(_grad);
			//std::cout << _grad << std::endl;
			//double objective1 = value(weights);
			//double step = 1e-3;

			//for (int wi = 0; wi < weights.size(); wi++)
			//{
			//	weights(wi) = weights(wi) + step;
			//	updateWeightsInGraph();
			//	double objective2 = value(weights);
			//	std::cout << (objective2 - objective1) / step << std::endl;
			//	weights(wi) = weights(wi) - step;
			//}
			//system("pause"); 

			double step_size = StyleSynthesisConfig::mContext_FunctionalityLearningStepSize;
			Eigen::MatrixXd grad;
			gradient(grad);
			weights = weights - step_size * grad;
			weights = (weights.array() < 0.0).select(zero_weights, weights);
			std::cout << weights << std::endl;

			objective = value(weights);  // the argument is not used, since it is a class member  updated by this function directly
			std::cout << "Iteration " << iter << " - objective: " << objective << ", #errors: " << numErrors() << std::endl;
			if (fabs(objective - old_objective) < 1e-6)
			{
				std::cout << "Converged / Stuck." << std::endl;
				std::cout << "Change of objective = " << fabs(objective - old_objective) << std::endl;
				break;
			}
			old_objective = objective;
		}
	}

	double ContextPartGraphLearning::value(const Eigen::VectorXd &weights_)
	{
		weights = weights_;
		updateWeightsInGraph();

		double objective = 0.0;
		double num_triplets = 0.0;
		int num_errors = 0;
		double l1regularization = StyleSynthesisConfig::mContext_FunctionalityLearningRegularization;

		for (int m = 0; m < (int)training_meshes_triplets.size(); m++)
		{
			Eigen::MatrixXd K;
			matched_graphs[m]->exportSimilarityMatrix(K);

			for (int t = 0; t < (int)training_meshes_triplets[m].size(); t++)
			{
				double Kab = K(std::get<0>(training_meshes_triplets[m][t]), std::get<1>(training_meshes_triplets[m][t]));
				double Kac = K(std::get<0>(training_meshes_triplets[m][t]), std::get<2>(training_meshes_triplets[m][t]));
				num_errors += int(Kac > Kab);
				objective -= logsigmoid(Kab - Kac);
				num_triplets++;

				int dd = 0;
				for (int d = 0; d < (int)matched_graphs[m]->mNodeWeights.size(); d++)
				{
					objective += l1regularization * abs(weights(dd));
					dd++;
				}
				for (int d = 0; d < (int)matched_graphs[m]->mNodeWeights.size(); d++)
				{
					dd++;
				}
				for (int d = 0; d < (int)matched_graphs[m]->mEdgeWeights.size(); d++)
				{
					objective += l1regularization * abs(weights(dd));
					dd++;
				}
			}
		}
		objective /= num_triplets;

		if (verbose_mode)
			std::cout << "Objective function called (#loss=" << objective << ", #err=" << num_errors << ")" << std::endl;

		return objective;
	}



	void ContextPartGraphLearning::gradient(Eigen::MatrixXd& dweights)
	{
		double l1regularization = StyleSynthesisConfig::mContext_FunctionalityLearningRegularization;
		dweights = Eigen::MatrixXd::Zero(weights.size(), 1);
		double num_triplets = 0.0;

		for (int m = 0; m < (int)training_meshes_triplets.size(); m++)
		{
			Eigen::MatrixXd K;
			matched_graphs[m]->exportSimilarityMatrix(K);
			vector< Eigen::MatrixXd > node_derivative_weights;
			vector< Eigen::MatrixXd > node_sigma_derivatives;
			vector< Eigen::MatrixXd > edge_derivative_weights;
			vector< Eigen::MatrixXd > edge_sigma_derivatives;
			vector< Eigen::MatrixXd > path_derivative_weights;
			matched_graphs[m]->exportDerivatives(node_derivative_weights, node_sigma_derivatives, edge_derivative_weights, edge_sigma_derivatives, path_derivative_weights);
			if (node_derivative_weights.empty() || node_sigma_derivatives.empty() || edge_derivative_weights.empty() || edge_sigma_derivatives.empty() || path_derivative_weights.empty() )
				return;

			for (int t = 0; t < (int)training_meshes_triplets[m].size(); t++)
			{
				double Kab = K(std::get<0>(training_meshes_triplets[m][t]), std::get<1>(training_meshes_triplets[m][t]));
				double Kac = K(std::get<0>(training_meshes_triplets[m][t]), std::get<2>(training_meshes_triplets[m][t]));
				double dk = sigmoid(Kab - Kac) - 1.0;
				int dd = 0;
				for (int d = 0; d < (int)matched_graphs[m]->mNodeWeights.size(); d++)
				{
					double dKab = node_derivative_weights[d](std::get<0>(training_meshes_triplets[m][t]), std::get<1>(training_meshes_triplets[m][t]));
					double dKac = node_derivative_weights[d](std::get<0>(training_meshes_triplets[m][t]), std::get<2>(training_meshes_triplets[m][t]));
					dweights(dd) += l1regularization * sgn( weights(dd) );
					dweights(dd++) += dk * (dKab - dKac);
				}
				for (int d = 0; d < (int)matched_graphs[m]->mNodeWeights.size(); d++)
				{
					double dKab = node_sigma_derivatives[d](std::get<0>(training_meshes_triplets[m][t]), std::get<1>(training_meshes_triplets[m][t]));
					double dKac = node_sigma_derivatives[d](std::get<0>(training_meshes_triplets[m][t]), std::get<2>(training_meshes_triplets[m][t]));
					dweights(dd++) += dk * (dKab - dKac);
				}
				for (int d = 0; d < (int)matched_graphs[m]->mEdgeWeights.size(); d++)
				{
					double dKab = edge_derivative_weights[d](std::get<0>(training_meshes_triplets[m][t]), std::get<1>(training_meshes_triplets[m][t]));
					double dKac = edge_derivative_weights[d](std::get<0>(training_meshes_triplets[m][t]), std::get<2>(training_meshes_triplets[m][t]));
					dweights(dd) += l1regularization * sgn( weights(dd) ); 
					dweights(dd++) += dk * (dKab - dKac);
				}
				for (int d = 0; d < (int)matched_graphs[m]->mEdgeWeights.size(); d++)
				{
					double dKab = edge_sigma_derivatives[d](std::get<0>(training_meshes_triplets[m][t]), std::get<1>(training_meshes_triplets[m][t]));
					double dKac = edge_sigma_derivatives[d](std::get<0>(training_meshes_triplets[m][t]), std::get<2>(training_meshes_triplets[m][t]));
					dweights(dd++) += dk * (dKab - dKac);
				}
				for (int d = 0; d < (int)matched_graphs[m]->mPathWeights.size(); d++)
				{
					double dKab = path_derivative_weights[d](std::get<0>(training_meshes_triplets[m][t]), std::get<1>(training_meshes_triplets[m][t]));
					double dKac = path_derivative_weights[d](std::get<0>(training_meshes_triplets[m][t]), std::get<2>(training_meshes_triplets[m][t]));
					dweights(dd++) += dk * (dKab - dKac);
				}
				num_triplets++;
			}
		}

		dweights /= num_triplets;
	}

	bool ContextPartGraphLearning::updateWeightsInGraph()
	{
		if (weights.size() == 0)
			return false;

		for (int m = 0; m < (int)training_meshes_triplets.size(); m++)
		{
			int dd = 0;
			for (int d = 0; d < (int)matched_graphs[m]->mNodeWeights.size(); d++)
				matched_graphs[m]->mNodeWeights[d] = weights(dd++);
			for (int d = 0; d < (int)matched_graphs[m]->mNodeWeights.size(); d++)
				matched_graphs[m]->mNodeSigmaMultipliers[d] = weights(dd++); 
			for (int d = 0; d < (int)matched_graphs[m]->mEdgeWeights.size(); d++)
				matched_graphs[m]->mEdgeWeights[d] = weights(dd++);
			for (int d = 0; d < (int)matched_graphs[m]->mEdgeWeights.size(); d++)
				matched_graphs[m]->mEdgeSigmaMultipliers[d] = weights(dd++);
			for (int d = 0; d < (int)matched_graphs[m]->mPathWeights.size(); d++)
				matched_graphs[m]->mPathWeights[d] = weights(dd++);

			if (!matched_graphs[m]->evaluateNodeSimilarity()) return false;
			if (!matched_graphs[m]->evaluateEdgeSimilarity()) return false;
			if (!matched_graphs[m]->computeGraphSimilarity()) return false;
			if (verbose_mode) cout << ".";
		}
		//if (verbose_mode) cout << endl;

		return true;
	}


	int ContextPartGraphLearning::numErrors()
	{
		int num_errors = 0;
		for (int m = 0; m < (int)training_meshes_triplets.size(); m++)
		{
			Eigen::MatrixXd K;
			matched_graphs[m]->exportSimilarityMatrix(K);

			for (int t = 0; t < (int)training_meshes_triplets[m].size(); t++)
			{
				double Kab = K(std::get<0>(training_meshes_triplets[m][t]), std::get<1>(training_meshes_triplets[m][t]));
				double Kac = K(std::get<0>(training_meshes_triplets[m][t]), std::get<2>(training_meshes_triplets[m][t]));
				num_errors += int(Kac > Kab);
			}
		}
		return num_errors;
	}

	double ContextPartGraphLearning::computeThresholdSimilarity()
	{
		if (!updateWeightsInGraph())
		{
			error("computeThresholdSimilarity() was called before weights are learned!");
			return -1.0;
		}

		vector<double> similarities;
		for (int m = 0; m < (int)training_meshes_triplets.size(); m++)
		{
			ContextPartGraphMatch cpgm1, cpgm2;
			cpgm1.mNodeWeights = matched_graphs[m]->mNodeWeights;
			cpgm1.mNodeSigmaMultipliers = matched_graphs[m]->mNodeSigmaMultipliers;
			cpgm1.mEdgeWeights = matched_graphs[m]->mEdgeWeights;
			cpgm1.mEdgeSigmaMultipliers = matched_graphs[m]->mEdgeSigmaMultipliers;
			cpgm1.mNodeSigma = matched_graphs[m]->mNodeSigma;
			cpgm1.mEdgeSigma = matched_graphs[m]->mEdgeSigma;
			cpgm1.mPathWeights = matched_graphs[m]->mPathWeights;

			cpgm2.mNodeWeights = matched_graphs[m]->mNodeWeights;
			cpgm2.mNodeSigmaMultipliers = matched_graphs[m]->mNodeSigmaMultipliers;
			cpgm2.mEdgeWeights = matched_graphs[m]->mEdgeWeights;
			cpgm2.mEdgeSigmaMultipliers = matched_graphs[m]->mEdgeSigmaMultipliers;
			cpgm2.mNodeSigma = matched_graphs[m]->mNodeSigma;
			cpgm2.mEdgeSigma = matched_graphs[m]->mEdgeSigma;
			cpgm2.mPathWeights = matched_graphs[m]->mPathWeights;

			cpgm1.loadGraph(*(matched_graphs[m]->mpSourceGraph), *(matched_graphs[m]->mpSourceGraph));
			cpgm2.loadGraph(*(matched_graphs[m]->mpTargetGraph), *(matched_graphs[m]->mpTargetGraph));
			cpgm1.process();
			cpgm2.process();

			Eigen::MatrixXd K1;
			cpgm1.exportSimilarityMatrix(K1);
			Eigen::MatrixXd K2;
			cpgm2.exportSimilarityMatrix(K2);

			Eigen::MatrixXd K;
			matched_graphs[m]->normalizeGraphSimilarity(K1, K2);
			matched_graphs[m]->exportSimilarityMatrix(K);

			for (int t = 0; t < (int)training_meshes_triplets[m].size(); t++)
			{
				double Kab = K(std::get<0>(training_meshes_triplets[m][t]), std::get<1>(training_meshes_triplets[m][t]));
				similarities.push_back(Kab);
			}
		}
		std::sort(similarities.begin(), similarities.end());

		//return similarities[similarities.size() - 1];
		// another possibility:
		 return similarities[ similarities.size() / 2 ];
	}

	double ContextPartGraphLearning::logsigmoid(double in)
	{
		in = min(in, 10.0);
		in = max(in, -10.0);

		return -log(1.0 + exp(-in));
	}

	double ContextPartGraphLearning::sigmoid(double in)
	{
		in = min(in, 10.0);
		in = max(in, -10.0);

		return 1.0 / (1.0 + exp(-in));
	}
}

