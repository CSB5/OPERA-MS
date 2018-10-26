#ifndef CLUSTER_H_
#define CLUSTER_H_

#include <stack>
#include <unordered_set>
#include <vector>
#include <algorithm>

#include "contig.h"

/**
 * @brief A class for representing clusters in hierarchical clustering trees.
 *
 * Represents clusters in hierarchical clustering trees.
 */
class Cluster {
public:
	/**
	 * @brief Constructs an empty cluster.
	 */
	Cluster();
	
	/**
	 * @brief Constructs a singleton cluster.
	 *
	 * @param contig	contig
	 */
	Cluster(Contig* contig);

	/**
	 * @brief Constructs a cluster node.
	 *
	 * @param child1	left child
	 * @param child2	right child
	 */
	Cluster(Cluster* child1, Cluster* child2);

	~Cluster(); /**< Default destructor. */


	void update_cluster_length_read_count(Cluster* child1, Cluster* child2);

	/**
	 * @brief Getter for contigs belonging to this cluster.
	 *
	 * @return contigs belonging to this cluster
	 */
	Contig** contigs() const;

	/**
	 * @brief Getter for number of contigs belonging to this cluster.
	 *
	 * @return number of contigs belonging to this cluster
	 */
	int num_contigs() const;

	/**
	 * @brief Getter for total length of contigs belonging to this cluster.
	 *
	 * @return total length of contigs belonging to this cluster
	 */
	int length() const;

	/**
	 * @brief Getter for sum of read counts for all samples.
	 *
	 * @return sum of read counts for all samples
	 */
	int* sum_read_counts() const;

	/**
	 * @brief Getter for arrival rates for all samples.
	 *
	 * @return arrival rates for all samples
	 */
	double* arrival_rates() const;
	double* mean_arrival_rates()const;
	double* median_arrival_rates()const;

	/**
	 * @brief Compute the median arrival rate base on contig bin
	 *
	 * @return void
	 */
	void compute_median_arrival_rate();

	/**
	 * @brief Getter for left child.
	 *
	 * @return left child
	 */
	Cluster* child1() const;

	/**
	 * @brief Getter for right child.
	 *
	 * @return right child
	 */
	Cluster* child2() const;

	/**
	 * @brief Getter for score.
	 *
	 * @return score
	 */
	double score() const;

	/**
	 * @brief Getter for model score.
	 *
	 * @return model score
	 */
	double model_score() const;

	/**
	 * @brief Tests whether the model is computed.
	 * 
	 * @return true if the model is computed, false otherwise
	 */
	bool modeled() const;

	/**
	 * @brief Tests whether this cluster is connected.
	 *
	 * @return true if this cluster is connected, false otherwise
	 */
	bool connected() const;

	/**
	 * @brief Setter for contigs belonging to this cluster.
	 *
	 * @param contigs	contigs belonging to this cluster
	 */
	void set_contigs(Contig** contigs);
	void set_contigs(std::vector< Contig*>* v_contig);

	/**
	 * @brief Setter for score.
	 *
	 * @param score score
	 */
	void set_score(double score);

	/**
	 * @brief Setter for model score.
	 *
	 * @param model_score model score
	 */
	void set_model_score(double model_score);

	/**
	 * @brief Sets a flag which indicates whether the model is computed.
	 *
	 * @param modeled true/false
	 */
	void set_modeled(bool modeled);

	/**
	 * @brief Sets a flag which indicates whether this cluster is connected.
	 *
	 * @param connected true/false
	 */
	void set_connected(bool connected);
	
	void set_cluster(Cluster* c1, Cluster* c2);
	
	//void set_iterator_un_set(std::unordered_set<Cluster*>::iterator it);
	//std::unordered_set<Cluster*>::iterator get_iterator_un_set()const;
	
	void set_ID(std::string* ID);
	std::string* get_ID()const;

private:
	Contig** contigs_; /**< Contigs belonging to this cluster. */

	int num_contigs_; /**< Number of contigs belonging to this cluster. */
	int length_; /**< Total length of contigs belonging to this cluster. */
	int* sum_read_counts_; /**< Sum of read counts for all samples. */
	double* mean_arrival_rates_; /**< Mean arrival rates over for all samples. */
	double* median_arrival_rates_; /**< Median contig arrival rates for all samples. */

	Cluster* child1_; /**< Left child. */
	Cluster* child2_; /**< Right child. */

	double score_; /**< Score. */
	double model_score_; /**< Model score. */
	bool modeled_; /**< A flag which indicates whether the model is computed. */
	bool connected_; /**< A flag which indicates whether this cluster is connected. */
	//std::unordered_set<Cluster*>::iterator it_un_set_;/**iterator to speed up the tree construction**/
	
	std::string* cluster_ID;

};


/**
 * A stack of clusters used for recursive operations on hierarchical clustering trees.
 */
typedef std::stack<Cluster*> ClusterStack;


/**
 * A set of clusters used for storing current roots of hierarchical clustering trees.
 */
typedef std::unordered_set<Cluster*> ClusterSet;

#endif // CLUSTER_H_
