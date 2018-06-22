#ifndef CLUSTER_GRAPH_H_
#define CLUSTER_GRAPH_H_

#include "contig.h"
#include "edge.h"
#include "cluster.h"
#include "probability_distribution.h"

/**
 * @brief A class for representing a graph of hierarchical clustering trees.
 *
 * Represents a greph of hierarchical clustering trees builit along the scaffold
 * edges based on contig arrival rate information.
 */
class ClusterGraph {
public:
	/**
	 * @brief Constructs a graph of hierarchical clustering trees.
	 *
	 * @param contigs	contigs
	 * @param edges		edges
	 */
	ClusterGraph(ContigMap* contigs, EdgeQueue* edges);

	/**
	 * @brief Constructs a graph of hierarchical clustering trees.
	 *
	 * @param contigs	contigs
	 * @param string	File name that contain a precomputed tree
	 */
	ClusterGraph(ContigMap* contigs, std::string tree_file);

	~ClusterGraph(); /**< Default destructor. */

	/**
	 * @brief Getter for roots of hierarchical clustering trees.
	 *
	 * @return roots of hierarchical clustering trees
	 */
	ClusterSet* roots();

	/**
	 * @brief Computes scores for all clusters based on given probability distribution.
	 *
	 * @param prob_dist		probability distribution
	 */
	void computeScores(const ProbabilityDistribution* prob_dist);


	void output_clusters(Cluster* cluster, double cs, double un_cs);
	/**
	 * @brief Computes models for all clustering trees maximizing BIC.
	 */
	void computeModels();

	/**
	 * @brief Saves final clusters to a file.
	 *
	 * @param clusters_file_path	path to file for saving final clusters
	 */
	void saveClusters(const char* clusters_file_path);

private:
	/**
	 * @brief Updates contig array pointers of all clusters.
	 */
	void updateClusters();

    const static int num_edges_to_resort = 50;

	/**
	 * @brief Computes score for the cluster based on given probability distribution.
	 *
	 * @param cluster		cluster
	 * @param prob_dist		probability distribution
	 */
	void computeClusterScore(Cluster* cluster, const ProbabilityDistribution* prob_dist);

	/**
	 * @brief Computes model for the cluster which maximizes BIC.
	 *
	 * @param cluster	cluster
	 */
	void computeClusterModel(Cluster* cluster);

	int num_contigs_; /**< Number of contigs. */
	double num_windows_; /**< Number of windows. */
	ClusterSet roots_; /**< Roots of hierarchical clustering trees. */
	
	/** Help methods  used during tree construction from a tree file*/
	Cluster* get_cluster_in_hash(std::string id, std::unordered_map<std::string, Cluster*>* c_map);
	void init_cluster(Cluster* cluster, std::vector<Contig*> *r_contig);
	//To dump all the output for debuging
	FILE *output_file_debug;
	
};

#endif // CLUSTER_GRAPH_H_
