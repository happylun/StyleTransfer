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

#include "PipelineTrainLearningIO.h"

#include <iostream>
#include <fstream>
#include <utility> 
#include <random>
#include <set>

#include "Data/StyleSynthesisConfig.h"
#include "Data/StyleSynthesisTypes.h"

#include "Segment/SegmentGroupApxCvx.h"

#include "Context/ContextPartGraph.h"
#include "Context/ContextPartGraphNodeGenerator.h"
#include "Context/ContextPartGraphMatch.h"
#include "Context/ContextPartGraphLearning.h"

#include "Mesh/MeshUtil.h"
#include "Sample/SampleUtil.h"
#include "Segment/SegmentUtil.h"
#include "Data/DataUtil.h"

using namespace std;
using namespace StyleSynthesis;

bool PipelineTrainLearningIO::process(string learningFolder) {
	string dataset_prefix = StyleSynthesisConfig::mData_DataSetRootFolder;
	string functionality_learning_path = dataset_prefix + learningFolder;
	string functionality_learning_filename = functionality_learning_path + "matching_pairs.txt";
	vector< pair<string, string> > training_meshes_name_pairs;
	vector < vector< pair<int, int> > > training_meshes_matching_node_ids_pairs;

	// build label map
	map<string, int> labelMap;
	int numLabels = 0;
	if (true) {
		vector<int> labelList;
		if (!DataUtil::loadIndexListASCII(functionality_learning_path + "mesh-labels-all.txt", labelList)) return false;
		ifstream meshListFile(functionality_learning_path + "mesh-list-all.txt");
		int meshID = 0;
		while (!meshListFile.eof()) {
			string line;
			getline(meshListFile, line);
			line = StringUtil::trim(line);
			if (meshListFile.fail() || line.empty()) break;
			string meshName = line.substr(line.find_first_of('/') + 1);
			int label = labelList[meshID] - 1;
			numLabels = max(numLabels, label + 1);
			labelMap[meshName] = label;
			meshID++;
		}
		meshListFile.close();
	}

	ifstream functionality_learning_file(functionality_learning_filename);
	while (!functionality_learning_file.eof()) {
		string line;
		getline(functionality_learning_file, line);
		if (line.empty()) break;

		vector<string> tokens;
		StringUtil::split(line, ' ', tokens);

		if (tokens.size() < 4)
		{
			std::cerr << "Each line in matching_pairs.txt should have two mesh filenames and at least one pair of matching part ids..." << std::endl;
			std::cerr << "Instead the last line read was: " << line << std::endl;
			return false;
		}

		// store the pair filenames
		training_meshes_name_pairs.push_back(make_pair(StringUtil::trim(tokens[0]), StringUtil::trim(tokens[1])));

		// read the matching part pairs of these two meshes
		vector< pair<int, int> > current_training_mesh_matching_node_ids_pairs;
		for (size_t pair_id = 2; pair_id < tokens.size(); pair_id += 2)
		{
			current_training_mesh_matching_node_ids_pairs.push_back(make_pair(stoi(tokens[pair_id]), stoi(tokens[pair_id + 1])));
		}
		training_meshes_matching_node_ids_pairs.push_back(current_training_mesh_matching_node_ids_pairs);
	}
	functionality_learning_file.close();

	// don't be scared of these very-high-dimensional data structures...
	int numLabelPairs = numLabels*numLabels;
	vector<vector<vector<tuple<int, int, int>>>> allTriplets(numLabelPairs); // triplet : # of triplet : # of mesh pairs : # of label pairs
	vector<vector<pair<ContextPartGraph*, ContextPartGraph*>>> allGraphs(numLabelPairs); // graph pair : # of mesh pairs : # of label pairs
	vector<vector<ContextPartGraphMatch*>> allMatchings(numLabelPairs); // matching : # of mesh pairs : # of label pairs

	ContextPartGraphNodeGenerator graph_node_generator;

	// generate mesh graphs, training data, translate segment ids to graph node ids
	for (int p = 0; p < (int)training_meshes_name_pairs.size(); p++)
	{
		ContextPartGraphMatch* cpgm = new ContextPartGraphMatch();
		pair<ContextPartGraph*, ContextPartGraph*> cpg;

		int labelPairID = -1;
		if (true) {
			int l1 = labelMap[training_meshes_name_pairs[p].first];
			int l2 = labelMap[training_meshes_name_pairs[p].second];
			if (l1 > l2) swap(l1, l2);
			labelPairID = l1*numLabels + l2;
		}

		std::cout << std::endl << "Loading meshes from pair " << p + 1 << " / " << training_meshes_name_pairs.size() << "...";
		if (loadGraphMatchStructures(training_meshes_name_pairs[p], functionality_learning_path, &graph_node_generator, cpgm, cpg))
		{			
			vector< tuple<int, int, int> > current_training_meshes_triplets;
			std::cout << " Associating training segment ids to graph node ids...";

			if (generateTrainingTriplets(training_meshes_matching_node_ids_pairs[p], cpgm, current_training_meshes_triplets))
			{
				// make sure that all of these will have the same size = #training meshes
				allMatchings[labelPairID].push_back(cpgm);
				allGraphs[labelPairID].push_back(cpg);
				allTriplets[labelPairID].push_back(current_training_meshes_triplets);
			}
			std::cout << " Done!" << std::endl;
		}	
		else
		{
			std::cerr << " Training pair will be ignored!" << std::endl;
			delete cpgm;
		}
	}

	Eigen::MatrixXd simThresholds(numLabels, numLabels);
	simThresholds.setZero();

	// do the learning!
	std::cout << std::endl << "Initiating training..." << std::endl;
	for (int labelID1 = 0; labelID1 < numLabels; labelID1++) {
		for (int labelID2 = labelID1; labelID2 < numLabels; labelID2++) {
			int labelPairID = labelID1*numLabels + labelID2;
			if (allMatchings[labelPairID].empty()) continue;

			string labelPairName = to_string(labelID1) + "--" + to_string(labelID2);
			cout << "||||||||||||||||||||||||||| learning " << labelPairName << " ||||||||||||||||||||||||||||||||" << endl;

			if (true) {
				// print out number of training examples
				int count = 0;
				for (auto &meshPair : allTriplets[labelPairID]) {
					count += (int)meshPair.size();
				}
				cout << count << " / " << allTriplets[labelPairID].size() << endl;
				//continue; // HACK: skip
			}

			string outputFolder = functionality_learning_path + labelPairName + "/";
			if (!FileUtil::makedir(outputFolder)) return false;
			if (FileUtil::existsfile(outputFolder + "weights.txt")) continue; // early quit

			ContextPartGraphLearning learning;
			learning.setTrainingData(allMatchings[labelPairID], allTriplets[labelPairID]);

			if (!learning.process())
			{
				std::cerr << "Training failed." << std::endl;
				return false;
			}
			learning.output(outputFolder);
			double thr = learning.computeThresholdSimilarity();
			std::cout << "Similarity threshold: " << thr << std::endl;
			simThresholds(labelID1, labelID2) = thr;
		}
	}
	//return true; // HACK: skip

	if (!DataUtil::saveMatrixASCII(functionality_learning_path + "thresholds.txt", simThresholds)) return false;

	// clean-up memory! (crash here?)
	for (int p = 0; p < (int)training_meshes_name_pairs.size(); p++)
	{
		//delete matched_graphs[p];
		//delete training_meshes_graphs[p].first;
		//delete training_meshes_graphs[p].second;
	}

	return true;
}



bool PipelineTrainLearningIO::generateTrainingTriplets(const vector< pair<int, int> >& training_mesh_matching_node_ids_pairs, const ContextPartGraphMatch* cpgm, vector< tuple<int, int, int> > & current_training_meshes_triplets)
{
	std::default_random_engine generator;

	for (int p = 0; p < (int)training_mesh_matching_node_ids_pairs.size(); p++)
	{
		int first_mesh_node_id = training_mesh_matching_node_ids_pairs[p].first;
		int second_mesh_node_id = training_mesh_matching_node_ids_pairs[p].second;

		auto second_mesh_graph = cpgm->mpTargetGraph;

		// mark nodes can be used for the third node
		int num_second_mesh_nodes = (int)second_mesh_graph->mAllNodes.size();
		vector<int> node_valid_flag(num_second_mesh_nodes, true);
		node_valid_flag[second_mesh_node_id] = false;
		for (auto node : second_mesh_graph->mAllNodes[second_mesh_node_id]->mSymmetry) {
			node_valid_flag[node->mUID] = false; // skip symmetric nodes
		}
		///////*
		//////for (int node_id = 0; node_id < num_second_mesh_nodes; node_id++) {
		//////	auto node = second_mesh_graph->mAllNodes[node_id];
		//////	if (node->mParent != second_mesh_graph->mRootNode) {
		//////		node_valid_flag[node->mUID] = false; // HACK: skip nodes not on first level
		//////	}
		//////}
		//////*/

		//// get candidate nodes
		vector<int> candidate_nodes_list;
		for (int node_id = 0; node_id < num_second_mesh_nodes; node_id++) {
			if (node_valid_flag[node_id]) 
				candidate_nodes_list.push_back(node_id);
		}
		if (candidate_nodes_list.empty()) {
			// no valid nodes on second mesh
			continue;
		}

		//// pick a VALID target (second) mesh node at random as dissimilar part
		std::uniform_int_distribution<int> distribution( 0, (int)candidate_nodes_list.size()-1 );
		int random_node_id = candidate_nodes_list[distribution(generator)];
		current_training_meshes_triplets.push_back( make_tuple(first_mesh_node_id, second_mesh_node_id, random_node_id) );

		////// pick a VALID target (second) mesh node with largest similarity to the first node as dissimilar part
		//double max_similarity = 0.0;
		//int most_erroneously_similar_node_id = -1;
		//for (int node_id = 0; node_id < num_second_mesh_nodes; node_id++) {
		//	if (node_valid_flag[node_id])
		//	{				
		//		if (cpgm->mGraphSimilarityMatrix(first_mesh_node_id, node_id) > max_similarity)
		//		{
		//			max_similarity = cpgm->mGraphSimilarityMatrix(first_mesh_node_id, node_id);					
		//			most_erroneously_similar_node_id = node_id;
		//		}
		//	}
		//}
		//if (most_erroneously_similar_node_id == -1)
		//	continue;
		//current_training_meshes_triplets.push_back(make_tuple(first_mesh_node_id, second_mesh_node_id, most_erroneously_similar_node_id));		
	}

	return true;
}


bool PipelineTrainLearningIO::loadGraphMatchStructures(const pair<string, string>& training_mesh_pairs, const string& functionality_learning_path, ContextPartGraphNodeGenerator* cpgng, ContextPartGraphMatch* cpgm, pair<ContextPartGraph*, ContextPartGraph*>& cpg)
{
	// prepare folder/file names
	string first_training_mesh_name = training_mesh_pairs.first;
	string second_training_mesh_name = training_mesh_pairs.second;
	string pair_folder = functionality_learning_path + "/" + first_training_mesh_name + "--" + second_training_mesh_name;
	std::cout << " Loading data from folder: " << pair_folder << "...";
	string first_training_mesh_folder = pair_folder + "/" + first_training_mesh_name + "/";
	string second_training_mesh_folder = pair_folder + "/" + second_training_mesh_name + "/";
	string first_training_mesh_filename = first_training_mesh_folder + "sourceMesh.ply";
	string second_training_mesh_filename = second_training_mesh_folder + "targetMesh.ply";
	string first_training_mesh_segments_filename = first_training_mesh_folder + "sourceSegment.txt";
	string second_training_mesh_segments_filename = second_training_mesh_folder + "targetSegment.txt";
	string first_training_mesh_graph_filename = first_training_mesh_folder + "sourceGraphHierarchy.txt";
	string second_training_mesh_graph_filename = second_training_mesh_folder + "targetGraphHierarchy.txt";
	string first_training_mesh_graph_descriptor_filename = first_training_mesh_folder + "sourceGraphDescriptor.txt";
	string second_training_mesh_graph_descriptor_filename = second_training_mesh_folder + "targetGraphDescriptor.txt";
	string first_training_mesh_graph_context_filename = first_training_mesh_folder + "sourceGraphContext.txt";
	string second_training_mesh_graph_context_filename = second_training_mesh_folder + "targetGraphContext.txt";
	string graph_matching_initial_weights_folder = functionality_learning_path + "/";

	// prepare data structures to store segments
	TTriangleMesh first_mesh;
	TTriangleMesh second_mesh;
	vector<TTriangleMesh> first_mesh_groups;
	vector<TTriangleMesh> second_mesh_groups;
	vector< vector<vector<int>> > first_mesh_segments;
	vector< vector<vector<int>> > second_mesh_segments;

	// load first mesh groups (largest parts)
	std::cout << " Loading mesh parts from: " << first_training_mesh_folder << "...";
	int num_first_mesh_parts = 0;
	while (true)
	{
		string part_name = first_training_mesh_folder + to_string(num_first_mesh_parts) + ".ply";
		if (!FileUtil::existsfile(part_name)) break;
		first_mesh_groups.push_back(TTriangleMesh());
		if (!MeshUtil::loadMesh(part_name, first_mesh_groups.back()))
		{
			std::cerr << "Failed to load mesh part " << part_name << std::endl;
			return false;
		}
		num_first_mesh_parts++;
	}

	// load second mesh groups (largest parts)
	std::cout << " Loading mesh parts from: " << second_training_mesh_folder << "...";
	int num_second_mesh_parts = 0;
	while (true)
	{
		string part_name = second_training_mesh_folder + to_string(num_second_mesh_parts) + ".ply";
		if (!FileUtil::existsfile(part_name)) break;
		second_mesh_groups.push_back(TTriangleMesh());
		if (!MeshUtil::loadMesh(part_name, second_mesh_groups.back()))
		{
			std::cerr << "Failed to load mesh part " << part_name << std::endl;
			return false;
		}
		num_second_mesh_parts++;
	}

	// load/generate first mesh segments
	std::cout << " Segmenting mesh " << first_training_mesh_name << " into components" << "...";
	if (FileUtil::existsfile(first_training_mesh_segments_filename)) {
		if (!MeshUtil::loadMesh(first_training_mesh_filename, first_mesh))
		{
			std::cerr << "Failed to load mesh " << first_training_mesh_filename << std::endl;
			return false;
		}
		if (!SegmentGroupApxCvx::loadSegments(first_training_mesh_segments_filename, first_mesh_segments))
		{
			std::cerr << "Failed to load mesh segments from " << first_training_mesh_segments_filename << std::endl;
			return false;
		}
	}
	else {
		if (!SegmentGroupApxCvx::runSegmentation(first_mesh_groups, first_mesh, first_mesh_segments))
		{
			std::cerr << "Failed to generate mesh segments for " << first_training_mesh_filename << std::endl;
			return false;
		}
		if (!SegmentGroupApxCvx::saveSegments(first_training_mesh_segments_filename, first_mesh_segments))
		{
			std::cerr << "Failed to save mesh segments in " << first_training_mesh_segments_filename << std::endl;
			return false;
		}
		if (!MeshUtil::saveMesh(first_training_mesh_filename, first_mesh))
		{
			std::cerr << "Failed to save mesh in " << first_training_mesh_filename << std::endl;
			return false;
		}
	}

	// load/generate second mesh segments 
	std::cout << " Segmenting mesh " << second_training_mesh_name << " into components" << "...";
	if (FileUtil::existsfile(second_training_mesh_segments_filename)) {
		if (!MeshUtil::loadMesh(second_training_mesh_filename, second_mesh))
		{
			std::cerr << "Failed to load mesh " << second_training_mesh_filename << std::endl;
			return false;
		}
		if (!SegmentGroupApxCvx::loadSegments(second_training_mesh_segments_filename, second_mesh_segments))
		{
			std::cerr << "Failed to load mesh segments from " << second_training_mesh_segments_filename << std::endl;
			return false;
		}
	}
	else {
		if (!SegmentGroupApxCvx::runSegmentation(second_mesh_groups, second_mesh, second_mesh_segments))
		{
			std::cerr << "Failed to generate mesh segments for " << second_training_mesh_filename << std::endl;
			return false;
		}
		if (!SegmentGroupApxCvx::saveSegments(second_training_mesh_segments_filename, second_mesh_segments))
		{
			std::cerr << "Failed to save mesh segments in " << second_training_mesh_segments_filename << std::endl;
			return false;
		}
		if (!MeshUtil::saveMesh(second_training_mesh_filename, second_mesh))
		{
			std::cerr << "Failed to save mesh in " << second_training_mesh_filename << std::endl;
			return false;
		}
	}

	// generate functionality graph for first mesh
	std::cout << " Generating graph for mesh " << first_training_mesh_name << "...";
	ContextPartGraph* first_training_mesh_graph = new ContextPartGraph();
	if (FileUtil::existsfile(first_training_mesh_graph_filename)) {
		if (!first_training_mesh_graph->loadGraphHierarchy(first_training_mesh_graph_filename, *cpgng, &first_mesh, &first_mesh_segments))
		{
			std::cerr << "Failed to load graph from " << first_training_mesh_graph_filename << std::endl;
			delete first_training_mesh_graph;
			return false;
		}
	}
	else {
		if (!first_training_mesh_graph->buildGraphHierarchy(*cpgng, &first_mesh, &first_mesh_segments))
		{
			std::cerr << "Failed to generate graph for mesh " << first_training_mesh_name << std::endl;
			delete first_training_mesh_graph;
			return false;
		}
		if (!first_training_mesh_graph->saveGraphHierarchy(first_training_mesh_graph_filename))
		{
			std::cerr << "Failed to save graph for mesh " << first_training_mesh_name << std::endl;
			delete first_training_mesh_graph;
			return false;
		}
	}
	if (FileUtil::existsfile(first_training_mesh_graph_descriptor_filename)) {
		if (!first_training_mesh_graph->loadGraphDescriptor(first_training_mesh_graph_descriptor_filename))
		{
			std::cerr << "Failed to load graph descriptors from " << first_training_mesh_graph_descriptor_filename << std::endl;
			delete first_training_mesh_graph;
			return false;
		}
	}
	else {
		if (!first_training_mesh_graph->buildGraphDescriptor())
		{
			std::cerr << "Failed to compute graph descriptors for mesh " << first_training_mesh_name << std::endl;
			delete first_training_mesh_graph;
			return false;
		}
		if (!first_training_mesh_graph->saveGraphDescriptor(first_training_mesh_graph_descriptor_filename))
		{
			std::cerr << "Failed to save graph descriptors in " << first_training_mesh_graph_descriptor_filename << std::endl;
			delete first_training_mesh_graph;
			return false;
		}
	}
	if (FileUtil::existsfile(first_training_mesh_graph_context_filename))
	{
		if (!first_training_mesh_graph->loadGraphContext(first_training_mesh_graph_context_filename))
		{
			std::cerr << "Failed to load graph context for mesh " << first_training_mesh_name << std::endl;
			delete first_training_mesh_graph;
			return false;
		}
	}
	else {
		if (!first_training_mesh_graph->buildGraphContext())
		{
			std::cerr << "Failed to generate graph context for mesh " << first_training_mesh_name << std::endl;
			delete first_training_mesh_graph;
			return false;
		}
		if (!first_training_mesh_graph->saveGraphContext(first_training_mesh_graph_context_filename))
		{
			std::cerr << "Failed to save graph context in " << first_training_mesh_graph_context_filename << std::endl;
			delete first_training_mesh_graph;
			return false;
		}
	}

	// generate functionality graph for second mesh
	std::cout << " Generating graph for mesh " << second_training_mesh_name << "...";
	ContextPartGraph* second_training_mesh_graph = new ContextPartGraph();
	if (FileUtil::existsfile(second_training_mesh_graph_filename)) {
		if (!second_training_mesh_graph->loadGraphHierarchy(second_training_mesh_graph_filename, *cpgng, &second_mesh, &second_mesh_segments))
		{
			std::cerr << "Failed to load graph from " << second_training_mesh_graph_filename << std::endl;
			delete second_training_mesh_graph;
			return false;
		}
	}
	else {
		if (!second_training_mesh_graph->buildGraphHierarchy(*cpgng, &second_mesh, &second_mesh_segments))
		{
			std::cerr << "Failed to generate graph for mesh " << second_training_mesh_name << std::endl;
			delete second_training_mesh_graph;
			return false;
		}
		if (!second_training_mesh_graph->saveGraphHierarchy(second_training_mesh_graph_filename))
		{
			std::cerr << "Failed to save graph for mesh " << second_training_mesh_name << std::endl;
			delete second_training_mesh_graph;
			return false;
		}
	}
	if (FileUtil::existsfile(second_training_mesh_graph_descriptor_filename)) {
		if (!second_training_mesh_graph->loadGraphDescriptor(second_training_mesh_graph_descriptor_filename))
		{
			std::cerr << "Failed to load graph descriptors from " << second_training_mesh_graph_descriptor_filename << std::endl;
			delete second_training_mesh_graph;
			return false;
		}
	}
	else {
		if (!second_training_mesh_graph->buildGraphDescriptor())
		{
			std::cerr << "Failed to compute graph descriptors for mesh " << second_training_mesh_name << std::endl;
			delete second_training_mesh_graph;
			return false;
		}
		if (!second_training_mesh_graph->saveGraphDescriptor(second_training_mesh_graph_descriptor_filename))
		{
			std::cerr << "Failed to save graph descriptors in " << second_training_mesh_graph_descriptor_filename << std::endl;
			delete second_training_mesh_graph;
			return false;
		}
	}
	if (FileUtil::existsfile(second_training_mesh_graph_context_filename))
	{
		if (!second_training_mesh_graph->loadGraphContext(second_training_mesh_graph_context_filename))
		{
			std::cerr << "Failed to load graph context for mesh " << second_training_mesh_name << std::endl;
			delete second_training_mesh_graph;
			return false;
		}
	}
	else {
		if (!second_training_mesh_graph->buildGraphContext())
		{
			std::cerr << "Failed to generate graph context for mesh " << second_training_mesh_name << std::endl;
			delete second_training_mesh_graph;
			return false;
		}
		if (!second_training_mesh_graph->saveGraphContext(second_training_mesh_graph_context_filename))
		{
			std::cerr << "Failed to save graph context in " << second_training_mesh_graph_context_filename << std::endl;
			delete second_training_mesh_graph;
			return false;
		}
	}


	std::cout << " Graph matching...";
	if (!cpgm->loadWeights(graph_matching_initial_weights_folder)) return false;
	if (!cpgm->loadGraph(*first_training_mesh_graph, *second_training_mesh_graph))
	{
		std::cout << "Failed to load graphs for this pair" << std::endl;
		delete first_training_mesh_graph;
		delete second_training_mesh_graph;
		return false;
	}
	if (!cpgm->process())
	{
		std::cout << "Failed to process graphs for this pair" << std::endl;
		delete first_training_mesh_graph;
		delete second_training_mesh_graph;
		return false;
	}

	cpg.first = first_training_mesh_graph;
	cpg.second = second_training_mesh_graph;

	return true;
}