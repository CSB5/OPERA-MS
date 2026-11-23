#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <utility>
#include <limits>
#include <vector>

#include "cluster_graph.h"

#include "sigma.h"

// Smallest contig ID (string order) of a cluster: the key that makes the
// output order independent of heap addresses.
static std::string smallest_contig_id(const Cluster* cluster) {
	std::string smallest = cluster->contigs()[0]->id();
	for (int contig_index = 1; contig_index < cluster->num_contigs(); ++contig_index) {
		if (cluster->contigs()[contig_index]->id() < smallest) {
			smallest = cluster->contigs()[contig_index]->id();
		}
	}
	return smallest;
}

static std::vector<Cluster*> sort_by_smallest_contig_id(const std::vector<Cluster*>& clusters) {
	std::vector<std::pair<std::string, Cluster*> > keyed;
	for (auto it = clusters.begin(); it != clusters.end(); ++it) {
		keyed.push_back(std::make_pair(smallest_contig_id(*it), *it));
	}
	std::sort(keyed.begin(), keyed.end(),
			[](const std::pair<std::string, Cluster*>& a, const std::pair<std::string, Cluster*>& b) { return a.first < b.first; });
	std::vector<Cluster*> sorted;
	for (auto it = keyed.begin(); it != keyed.end(); ++it) {
		sorted.push_back(it->second);
	}
	return sorted;
}

std::vector<Cluster*> ClusterGraph::sortedRoots() const {
	return sort_by_smallest_contig_id(std::vector<Cluster*>(roots_.begin(), roots_.end()));
}

ClusterGraph::ClusterGraph(ContigMap* contigs, EdgeQueue* edges) {
	num_contigs_ = (int) contigs->size();
	num_windows_ = 0;

	for (auto it = contigs->begin(); it != contigs->end(); ++it) {
		sContig* contig = (*it).second;

		//num_windows_ += contig->num_windows();

		roots_.insert(new Cluster(contig));
	}
	num_windows_ = (double)Sigma::total_assembly_size / (double)Sigma::contig_window_len;

	// Single linkage along the edges. The current root of a contig is found
	// with a union-find (path compression on Cluster::uf_parent); contig
	// arrays and contig -> cluster links are set once at the end, which
	// keeps the construction near linear even for trees of millions of
	// contigs.
	while (!edges->empty()) {
		Edge edge = edges->top();

		Cluster* cluster1 = find_root(edge.contig1()->cluster());
		Cluster* cluster2 = find_root(edge.contig2()->cluster());

		if (cluster1 != cluster2) {
			Cluster* merged = new Cluster(cluster1, cluster2);
			cluster1->set_uf_parent(merged);
			cluster2->set_uf_parent(merged);
			roots_.insert(merged);

			roots_.erase(cluster1);
			roots_.erase(cluster2);
		}

		edges->pop();
		/*
		//Average linkage !!!
		for(int i = 0; i < (int) edges->size(); ++i){
			
			if(edges->size() < num_edges_to_resort){
				break;
			}

			Edge resort_edge = edges->top();
			edges->pop();
			resort_edge.computeDistanceCluster();
			
			edges->push(resort_edge);
		}
		*/
	}

	buildRootContigs();
	updateClusters();
}

Cluster* ClusterGraph::find_root(Cluster* cluster) {
	Cluster* root = cluster;
	while (root->uf_parent() != NULL) {
		root = root->uf_parent();
	}
	while (cluster != root) {
		Cluster* next = cluster->uf_parent();
		cluster->set_uf_parent(root);
		cluster = next;
	}
	return root;
}

void ClusterGraph::buildRootContigs() {
	// The contigs of a root in the order of the former copying merges:
	// child1's contigs, then child2's, recursively.
	std::vector<Cluster*> stack;
	for (auto it = roots_.begin(); it != roots_.end(); ++it) {
		Cluster* root = *it;
		if (root->num_contigs() == 1) {
			continue;
		}

		sContig** contigs = new sContig*[root->num_contigs()];
		int contig_index = 0;
		stack.push_back(root);
		while (!stack.empty()) {
			Cluster* cluster = stack.back();
			stack.pop_back();
			if (cluster->num_contigs() == 1) {
				// a leaf: take its contig, free its one-element array
				contigs[contig_index++] = cluster->contigs()[0];
				delete[] cluster->contigs();
				cluster->set_contigs((sContig**) NULL);
			} else {
				stack.push_back(cluster->child2());
				stack.push_back(cluster->child1());
			}
		}
		root->set_contigs(contigs);

		for (int i = 0; i < root->num_contigs(); ++i) {
			contigs[i]->set_cluster(root);
		}
	}
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
	std::vector<sContig*> *r_contig = new std::vector<sContig*>();
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


void ClusterGraph::init_cluster(Cluster* cluster, std::vector<sContig*> *v_contig){
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

	if (Sigma::USE_WINDOW == 1 && Sigma::num_samples == 1 && Sigma::pdist_type == "NegativeBinomial") {
		computeScoresFromWindowStats(prob_dist);
		return;
	}

	std::vector<Cluster*> roots = sortedRoots();
	for (auto it = roots.begin(); it != roots.end(); ++it) {
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


// With windows, one sample and the negative binomial, computeClusterScore
// adds for each window w of the cluster (k_w its read count, mu = arrival
// rate * window length and R the dispersion, the same for all windows):
//   R log(R / (R + mu)) + lgamma(R + k_w) + k_w log(mu / (R + mu))
//   - lgamma(k_w + 1) - lgamma(R)
// which sums to
//   W (R log(R / (R + mu)) - lgamma(R)) + L + K log(mu / (R + mu))
// with W the number of windows, K the sum of the k_w and L the sum of
// lgamma(R + k_w) - lgamma(k_w + 1). W, K and L of a node are the sums of
// its children's, so each score takes constant time instead of time
// proportional to the size of the cluster (quadratic over a deep tree).
// The value only differs from the per-window sum by rounding.
void ClusterGraph::computeScoresFromWindowStats(const ProbabilityDistribution* prob_dist) {
	const double R = Sigma::R_VALUE;

	// post-order: children before their parent
	std::vector<std::pair<Cluster*, bool> > stack;
	std::vector<Cluster*> roots = sortedRoots();
	for (auto it = roots.begin(); it != roots.end(); ++it) {
		stack.push_back(std::make_pair(*it, false));
	}

	while (!stack.empty()) {
		Cluster* cluster = stack.back().first;
		const bool children_done = stack.back().second;
		stack.pop_back();

		if (cluster->num_contigs() == 1) {
			sContig* contig = cluster->contigs()[0];
			double lgamma_sum = 0.0;
			long long reads = 0;
			for (int window_index = 0; window_index < contig->num_windows(); ++window_index) {
				const double k = round((double) contig->read_counts()[0][window_index]);
				lgamma_sum += lgamma(R + k) - lgamma(k + 1);
				reads += (long long) k;
			}
			cluster->set_window_stats(lgamma_sum, contig->num_windows(), reads);
		} else if (!children_done) {
			stack.push_back(std::make_pair(cluster, true));
			stack.push_back(std::make_pair(cluster->child2(), false));
			stack.push_back(std::make_pair(cluster->child1(), false));
			continue;
		} else {
			const Cluster* c1 = cluster->child1();
			const Cluster* c2 = cluster->child2();
			cluster->set_window_stats(c1->window_lgamma_sum() + c2->window_lgamma_sum(),
					c1->window_num() + c2->window_num(),
					c1->window_reads() + c2->window_reads());
		}

		computeClusterScoreFromWindowStats(cluster, prob_dist);
	}
}

void ClusterGraph::computeClusterScoreFromWindowStats(Cluster* cluster, const ProbabilityDistribution* prob_dist) {
	const double R = Sigma::R_VALUE;
	const double mu = cluster->arrival_rates()[0] * Sigma::contig_window_len;
	const long long W = cluster->window_num();
	const long long K = cluster->window_reads();

	double score;
	if (W == 0) {
		score = 0.0;
	} else if (mu > 0.0 && std::isfinite(mu)) {
		score = (double) W * (R * log(R / (R + mu)) - lgamma(R))
				+ cluster->window_lgamma_sum()
				+ (double) K * log(mu / (R + mu));
	} else if (mu == 0.0 && K > 0) {
		// not expected (a zero mean with reads in the windows): the
		// per-window sum (-inf or NaN depending on the windows)
		computeClusterScore(cluster, prob_dist);
		return;
	} else {
		// as in computeClusterScore, every window adds NaN: 0 * log(0) for
		// a zero mean, the log of a negative number or of NaN otherwise
		score = std::numeric_limits<double>::quiet_NaN();
	}

	score -= 0.5 * Sigma::num_samples * log(num_windows_);
	cluster->set_score(score);
}

void ClusterGraph::output_clusters(Cluster* cluster, double cs, double un_cs) {
	
	//std::string output_clusters_file_path = Sigma::output_dir + "/scored_clusters.dat";
	//FILE *output_file = fopen(output_clusters_file_path.c_str(), "a");
	
	//fprintf(output_file, "%d\t%d\t%d\t%f\t", cluster, cluster->child1(), cluster->child2(), cluster->score());
	fprintf(output_file_debug, "%s\t%s\t%s\t%f\t%f\t", (*(cluster->get_ID())).c_str(), (*(cluster->child1()->get_ID())).c_str(), (*(cluster->child2()->get_ID())).c_str(), cs, un_cs);

	//To write the contig idea in case of a leaf
	if(cluster->num_contigs() == 1){
		sContig* contig = cluster->contigs()[0];
		fprintf(output_file_debug, "%s\n", contig->id().c_str());
	}
	else{
		fprintf(output_file_debug, "-\n");
	}
	//	for (int contig_index = 0; contig_index < cluster->num_contigs(); ++contig_index) {
	//		sContig* contig = cluster->contigs()[contig_index];
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
			sContig* contig = cluster->contigs()[contig_index];

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
		ClusterStack clusters;

		for (auto it = roots_.begin(); it != roots_.end(); ++it) {
			if(Sigma::COMPUTE_SCORE == 0){
				(*it)->set_connected(true);
			}
			clusters.push(*it);
		}

		// collect the clusters where the trees are cut ...
		std::vector<Cluster*> connected;
		while (!clusters.empty()) {
			Cluster* cluster = clusters.top();
			clusters.pop();

			if (cluster->connected()) {
				connected.push_back(cluster);
			} else {
				clusters.push(cluster->child1());
				clusters.push(cluster->child2());
			}
		}

		// ... and number them by their smallest contig ID
		connected = sort_by_smallest_contig_id(connected);
		int cluster_id = 1;
		for (auto it = connected.begin(); it != connected.end(); ++it) {
			Cluster* cluster = *it;
			//Output the cluster ID were the cut have been made
			//output_clusters(cluster, -1, -1);
			for (int contig_index = 0; contig_index < cluster->num_contigs(); ++contig_index) {
				sContig* contig = cluster->contigs()[contig_index];

				fprintf(clusters_fp, "%s\t%d\t%d\t%f\n",
						contig->id().c_str(), cluster_id, contig->sum_read_counts()[0], cluster->arrival_rates()[0]);
			}

			cluster_id++;
		}

		fclose(clusters_fp);
		fclose(output_file_debug);
	} else {
		fprintf(stderr, "Error opening file: %s\n", clusters_file_path);
		exit(EXIT_FAILURE);
	}
}
