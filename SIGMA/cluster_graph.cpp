#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

#include "cluster_graph.h"

#include "sigma.h"

ClusterGraph::ClusterGraph(ContigMap* contigs, EdgeQueue* edges) {
	num_contigs_ = (int) contigs->size();
	num_windows_ = 0;

	for (auto it = contigs->begin(); it != contigs->end(); ++it) {
		Contig* contig = (*it).second;

		//num_windows_ += contig->num_windows();

		roots_.insert(new Cluster(contig));
	}
	num_windows_ = (double)Sigma::total_assembly_size / (double)Sigma::contig_window_len;

	while (!edges->empty()) {
		Edge edge = edges->top();

		Cluster* cluster1 = edge.contig1()->cluster();
		Cluster* cluster2 = edge.contig2()->cluster();

		if (cluster1 != cluster2) {
			roots_.insert(new Cluster(cluster1, cluster2));

			roots_.erase(cluster1);
			roots_.erase(cluster2);
		}

		edges->pop();

        for(int i = 0; i < (int) edges->size(); ++i){

            if(edges->size() < num_edges_to_resort){
                break;
            }

            Edge resort_edge = edges->top();
            edges->pop();
            resort_edge.computeDistanceCluster();

            edges->push(resort_edge);
        }
	}

	updateClusters();
}


ClusterGraph::ClusterGraph(ContigMap* contigs, std::string tree_file) {
	//Contruct the hash of each cluster ID to their pointer
	std::unordered_map<std::string, Cluster*> cluster_map;
	std::unordered_map<Cluster*, bool> father_map;

	num_contigs_ = (int) contigs->size(); 
	num_windows_ = 0;
	
	FILE* t_f = fopen(tree_file.c_str(), "r");
	
	char f_id[256];
	char s_id_1[256];
	char s_id_2[256];
	char c_id[256];
	
	double s_1, s_2;
	Cluster* father;
	fprintf(stderr, " **** Reading tree file %s\n", tree_file.c_str());
	while (!feof(t_f)) {
		//320468544	168318128	167063312	-24.914767	-31.573449	-
		if (fscanf(t_f, "%s\t%s\t%s\t%lf\t%lf\t%s\n", f_id, s_id_1, s_id_2, &s_1, &s_2, c_id) == 6) {
			
			if(s_1 == -1){
				continue;
			}
			
			//fprintf(stderr, "Reading tree file cluster: %s\n", f_id);

			//it is a contig
			if(c_id[0] != '-'){
				//fprintf(stderr, "Updating contig: %s -> %s\n", f_id, c_id);std::cin.get();
				auto it = contigs->find(c_id);
				father = new Cluster(it->second);
				father->set_ID(new std::string(f_id));
				cluster_map[f_id] = father;
			}
			//It is an internal node
			else{
				//Search for cluster1
				Cluster* cluster1 = get_cluster_in_hash(s_id_1, &cluster_map);
				
				//Search for cluster1
				Cluster* cluster2 = get_cluster_in_hash(s_id_2, &cluster_map);
				
				//Search for the father
				father = get_cluster_in_hash(f_id, &cluster_map);
				father->set_cluster(cluster1, cluster2);
				
				//To keep track of the roots
				father_map[cluster1] = false;
				father_map[cluster2] = false;
			}
			
			//update the root value for the father node
			auto it = father_map.find(father);
			if(it ==  father_map.end()){
				father_map[father] = true;
			}
		}
	}
	
	num_windows_ = (double)Sigma::total_assembly_size / (double)Sigma::contig_window_len;

	//For each cluster roots initialize the clusters and the root contig set
	std::vector<Contig*> *r_contig = new std::vector<Contig*>();
	for (auto it = father_map.begin(); it != father_map.end(); ++it) {
		//This is a root
		if(it->second){
			Cluster* root = it->first;
			
			roots_.insert(root);
			
			//fprintf(stderr, " **** init_cluster next root \n");	
			init_cluster(root, r_contig);
			root->set_contigs(r_contig);
			r_contig->clear();
		}
	}
	delete r_contig;
	updateClusters();
}


void ClusterGraph::init_cluster(Cluster* cluster, std::vector<Contig*> *v_contig){
	//It is a leaf
	//fprintf(stderr, "\t **** init_cluster %s\n", (cluster->get_ID())->c_str());	
	if(cluster->num_contigs() == 1){
		v_contig->push_back(cluster->contigs()[0]);
	}
	else{
		Cluster* c1 = cluster->child1();
		if(c1 != NULL){
			init_cluster(c1, v_contig);	
		}
		else{
			fprintf(stderr, "\t\t **** init_cluster null c1 %s\n", (cluster->get_ID())->c_str());	
		}
		
		Cluster* c2 = cluster->child2();
		if(c2 != NULL){
			init_cluster(c2, v_contig);	
		}
		else{
			fprintf(stderr, "\t\t **** init_cluster null c2 %s\n", (cluster->get_ID())->c_str());	
		}
		
		cluster->update_cluster_length_read_count(c1, c2);
	}
}

Cluster* ClusterGraph::get_cluster_in_hash(std::string id, std::unordered_map<std::string, Cluster*>* c_map){
	Cluster* cluster;
	auto it = c_map->find(id);
	if(it == c_map->end()){
		cluster = new Cluster();
		cluster->set_ID(new std::string(id));
		(*c_map)[id] = cluster;
	}
	else{
		cluster = it->second;
	}
	return cluster;
}


void ClusterGraph::updateClusters() {
	ClusterStack clusters;

	for (auto it = roots_.begin(); it != roots_.end(); ++it) {
		clusters.push(*it);
	}

	while (!clusters.empty()) {
		Cluster* cluster = clusters.top();
		clusters.pop();

		if (cluster->num_contigs() > 1) {
			cluster->child1()->set_contigs(cluster->contigs());
			cluster->child2()->set_contigs(cluster->contigs() + cluster->child1()->num_contigs());

			clusters.push(cluster->child1());
			clusters.push(cluster->child2());
		}
	}
}

ClusterGraph::~ClusterGraph() {
	ClusterStack clusters;

	for (auto it = roots_.begin(); it != roots_.end(); ++it) {
		Cluster* cluster = *it;

		for (int contig_index = 0; contig_index < cluster->num_contigs(); ++contig_index) {
			delete cluster->contigs()[contig_index];
		}

		delete[] cluster->contigs();

		clusters.push(cluster);
	}

	while (!clusters.empty()) {
		Cluster* cluster = clusters.top();
		clusters.pop();

		if (cluster->num_contigs() > 1) { 
			clusters.push(cluster->child1());
			clusters.push(cluster->child2());
		}

		delete cluster;
	}
}

ClusterSet* ClusterGraph::roots() { return &roots_; }

void ClusterGraph::computeScores(const ProbabilityDistribution* prob_dist) {
	ClusterStack clusters;


	std::string output_clusters_file_path = Sigma::output_dir + "/scored_clusters.dat";
	fprintf(stderr, "Outputing scored clusters in %s\n", output_clusters_file_path.c_str());
	output_file_debug = fopen(output_clusters_file_path.c_str(), "w");
	//fclose(output_file_debug);

	for (auto it = roots_.begin(); it != roots_.end(); ++it) {
		clusters.push(*it);
	}

	while (!clusters.empty()) {
		Cluster* cluster = clusters.top();
		clusters.pop();

		if (cluster->num_contigs() > 1) {
			clusters.push(cluster->child1());
			clusters.push(cluster->child2());
		}

		computeClusterScore(cluster, prob_dist);
	}
}


void ClusterGraph::output_clusters(Cluster* cluster, double cs, double un_cs) {
	
	//std::string output_clusters_file_path = Sigma::output_dir + "/scored_clusters.dat";
	//FILE *output_file = fopen(output_clusters_file_path.c_str(), "a");
	
	//fprintf(output_file, "%d\t%d\t%d\t%f\t", cluster, cluster->child1(), cluster->child2(), cluster->score());
	fprintf(output_file_debug, "%s\t%s\t%s\t%f\t%f\t", (*(cluster->get_ID())).c_str(), (*(cluster->child1()->get_ID())).c_str(), (*(cluster->child2()->get_ID())).c_str(), cs, un_cs);

	//To write the contig idea in case of a leaf
	if(cluster->num_contigs() == 1){
		Contig* contig = cluster->contigs()[0];
		fprintf(output_file_debug, "%s\n", contig->id().c_str());
	}
	else{
		fprintf(output_file_debug, "-\n");
	}
	//	for (int contig_index = 0; contig_index < cluster->num_contigs(); ++contig_index) {
	//		Contig* contig = cluster->contigs()[contig_index];
	//		fprintf(output_file, "%s\t", contig->id().c_str());		
	//	}

	//	fprintf(output_file, "\n");

	//fclose(output_file);
	
	return;
}

void ClusterGraph::computeClusterScore(Cluster* cluster, const ProbabilityDistribution* prob_dist) {
	double score = 0;

	for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
		double cluster_read_count = 0.0;
		double cluster_dispersion = 0.0;
		if (Sigma::USE_WINDOW == 1) {
			cluster_read_count = cluster->arrival_rates()[sample_index] * Sigma::contig_window_len;
			cluster_dispersion = Sigma::R_VALUE;
		}

		for (int contig_index = 0; contig_index < cluster->num_contigs(); ++contig_index) {
			Contig* contig = cluster->contigs()[contig_index];

			if (Sigma::USE_WINDOW == 0) {
				cluster_read_count = cluster->arrival_rates()[sample_index] * contig->modified_length();
				cluster_dispersion = Sigma::R_VALUE*contig->modified_length();
				//fprintf(stdout, " *** global R: %f, cluster_dispersion: %f\n", Sigma::R, cluster_dispersion);
				score += prob_dist->logpf(cluster_read_count, cluster_dispersion, contig->sum_read_counts()[sample_index]);
			} else {
				for (int window_index = 0; window_index < contig->num_windows(); ++window_index) {
					score += prob_dist->logpf(cluster_read_count, cluster_dispersion, contig->read_counts()[sample_index][window_index]);
				}
			}
		}
	}

	//Can add a penalty factor X: score * X
	if(Sigma::USE_WINDOW == 0){
		score -= 0.5 * Sigma::num_samples * log((double)Sigma::total_assembly_nb_contig);
	}
	else{
		score -= 0.5 * Sigma::num_samples * log(num_windows_);
	}

	cluster->set_score(score);
}

void ClusterGraph::computeModels() {
	ClusterStack clusters;
	
	fprintf(stderr, " *** computeModels: Init Root ...\n");
	for (auto it = roots_.begin(); it != roots_.end(); ++it) {
		clusters.push(*it);
	}
	
	fprintf(stderr, " *** computeModels: computeClusterModel ...\n");
	while (!clusters.empty()) {
		Cluster* cluster = clusters.top();

		if (cluster->num_contigs() == 1 || (cluster->child1()->modeled() && cluster->child2()->modeled())) {
			clusters.pop();
			computeClusterModel(cluster);
		} else {
			clusters.push(cluster->child1());
			clusters.push(cluster->child2());
		}
	}
	
	fprintf(stderr, " *** computeModels: Init Root (2) ...\n");
	for (auto it = roots_.begin(); it != roots_.end(); ++it) {
		clusters.push(*it);
	}

	fprintf(stderr, " *** computeModels: Associate contig to cluster ...\n");
	while (!clusters.empty()) {
		Cluster* cluster = clusters.top();
		clusters.pop();

		if (cluster->connected()) {
			for (int contig_index = 0; contig_index < cluster->num_contigs(); ++contig_index) {
				cluster->contigs()[contig_index]->set_cluster(cluster);
			}
		} else {
			clusters.push(cluster->child1());
			clusters.push(cluster->child2());
		}
	}
}

void ClusterGraph::computeClusterModel(Cluster* cluster) {
	if (cluster->num_contigs() == 1) {
		cluster->set_model_score(cluster->score());
		cluster->set_connected(true);
		//To correct to get the full true
		//output_clusters(cluster, cluster->score(), cluster->score());
	} else {
		double connected_score = cluster->score();
		double disconnected_score = cluster->child1()->model_score() + cluster->child2()->model_score();

		//output_clusters(cluster, connected_score, disconnected_score);

		if (connected_score >= disconnected_score){// - Sigma::SPLITTING_PENALTY) {
			cluster->set_model_score(connected_score);
			cluster->set_connected(true);
		} else {
			cluster->set_model_score(disconnected_score);
			cluster->set_connected(false);
		}
	}

	cluster->set_modeled(true);
}

void ClusterGraph::saveClusters(const char* clusters_file_path) {
	FILE* clusters_fp = fopen(clusters_file_path, "w");

	if (clusters_fp != NULL) {
		int cluster_id = 1;

		ClusterStack clusters;
		
		for (auto it = roots_.begin(); it != roots_.end(); ++it) {
			if(Sigma::COMPUTE_SCORE == 0){
				(*it)->set_connected(true);
			}
			clusters.push(*it);
		}
		

		while (!clusters.empty()) {
			Cluster* cluster = clusters.top();
			clusters.pop();

			if (cluster->connected()) {
				//Output the cluster ID were the cut have been made
				//output_clusters(cluster, -1, -1);
				for (int contig_index = 0; contig_index < cluster->num_contigs(); ++contig_index) {
					Contig* contig = cluster->contigs()[contig_index];

					fprintf(clusters_fp, "%s\t%d\t%d\t%f\n",
							contig->id().c_str(), cluster_id, contig->sum_read_counts()[0], cluster->arrival_rates()[0]);
				}

				cluster_id++;
			} else {
				clusters.push(cluster->child1());
				clusters.push(cluster->child2());
			}
		}

		fclose(clusters_fp);
		fclose(output_file_debug);
	} else {
		fprintf(stderr, "Error opening file: %s\n", clusters_file_path);
		exit(EXIT_FAILURE);
	}
}
