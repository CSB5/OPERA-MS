#ifndef EDGE_H_
#define EDGE_H_

#include <vector>
#include <queue>
#include <unordered_set>
#include <functional>

#include "contig.h"

/**
 * @brief A class for representing scaffold edges in the assembly graph.
 *
 * Represents a scaffold edge which includes two contigs as well as their
 * distance based on their arrival rates.
 */
class Edge {
public:
	/**
	 * @brief Constructs an edge.
	 *
	 * Constructs an edge between two given contigs.
	 *
	 * @param contig1	first contig
	 * @param contig2	second contig
	 */
	Edge(Contig* contig1, Contig* contig2);

	/**
	 * @brief Getter for first contig.
	 *
	 * @return first contig
	 */
	Contig* contig1() const;

	/**
	 * @brief Getter for second contig.
	 *
	 * @return second contig
	 */
	Contig* contig2() const;

	/**
	 * @brief Getter for distance.
	 *
	 * @return distance
	 */
	double distance() const;

	/**
	 * @brief Computes distance between contigs incident to this edge.
	 */
	void computeDistance();

    void computeDistanceCluster();

	/**
	 * @brief Tests if the two given edges are equal.
	 *
	 * Tests if the two given edges are equal i.e. contigs incident to them
	 * have same ids.
	 *
	 * @param edge1		first edge
	 * @param edge2		second edge
	 * @return true if the edges are equal, otherwise false
	 */
	friend bool operator==(const Edge& edge1, const Edge& edge2);

private:
	Contig* contig1_; /**< First contig. */
	Contig* contig2_; /**< Second contig. */
	double distance_; /**< Distance. */
};


/**
 * @brief Function object class for less-than inequality comparison of Edge objects.
 *
 * Function object class for less-than inequality comparison of scaffold edges based on
 * their distances. It is used as a comparator for the priority queue of edges.
 */
class EdgeComparator {
public:
	/**
	 * @brief Computes if the first edge compares less than the second edge.
	 *
	 * Computes if the first edge compares less than the second edge (i.e.
	 * whether the first edge has larger distance then the second edge
	 * and thus a lower priority).
	 * 
	 * @param edge1		first edge
	 * @param edge2		second edge
	 * @return true if the first edge has larger distance than the second edge, false otherwise
	 */
	bool operator()(const Edge& edge1, const Edge& edge2) const;
};


/**
 * @brief Function object class for computing hash values for Edge objects.
 *
 * Function object class for computing hash values for scaffold edges based
 * on ids of contigs incident to the edges.
 */
class EdgeHash {
public:
	/**
	 * @brief Computes hash value for given edge.
	 *
	 * @param edge	edge
	 * @return hash value for given edge
	 */
    size_t operator()(const Edge& edge) const;

private:
	std::hash<std::string> string_hash_; /**< Default hash function object for strings. */
};


/** A hash set of edges used for storing unique edges when reading edges files. */
typedef std::unordered_set<Edge, EdgeHash> EdgeSet;


/** A priority queue of edges used for building clustering trees. */
typedef std::priority_queue<Edge, std::vector<Edge>, EdgeComparator> EdgeQueue;

#endif // EDGE_H_
