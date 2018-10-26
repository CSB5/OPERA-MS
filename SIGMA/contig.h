#ifndef CONTIG_H_
#define CONTIG_H_

#include <string>
#include <vector>
#include <unordered_map>

class Cluster;


/**
 * @brief A class for representing contigs in the assembly graph.
 *
 * Represents a contig in the assembly graph including additional information
 * obtained from read mapping.
 */
class Contig {
public:
	/**
	 * @brief Constructs a contig.
	 *
	 * @param id		id
	 * @param length	length
	 */
	Contig(std::string id, int length);

	/**
	 * @brief Constructs a contig from previously extracted contig information.
	 *
	 * @param id			id
	 * @param length		length
	 * @param left_edge		starting point of the first window
	 * @param right_edge	ending point of the last window
	 * @param num_windows	number of windows
	 */
	Contig(std::string id, int length, int left_edge, int right_edge, int num_windows);

	~Contig(); /**< Default destructor. */

	/**
	 * @brief Getter for id.
	 * 
	 * @return id
	 */
	std::string id() const;

	/**
	 * @brief Getter for length.
	 *
	 * @return length
	 */	
	int length() const;

	/**
	 * @brief Getter for modifed length.
	 *
	 * @return modifed length
	 */	
	int modified_length() const;

	/**
	 * @brief Getter for starting point of the first window.
	 *
	 * @return starting point of the first window
	 */	
	int left_edge() const;

	/**
	 * @brief Getter for ending point of the last window.
	 *
	 * @return ending point of the last window
	 */	
	int right_edge() const;

	/**
	 * @brief Getter for number of windows.
	 *
	 * @return number of windows
	 */	
	int num_windows() const;

	/**
	 * @brief Getter for sum of read counts for all samples.
	 *
	 * @return sum of read counts for all samples
	 */
	int* sum_read_counts() const;

	/**
	 * @brief Getter for sum of read counts for all windows.
	 *
	 * @return sum of read counts for all windows
	 */
	int** read_counts() const;

	/**
	 * @brief Getter for cluster containing this contig.
	 *
	 * During the construction of clustering trees, this is the root
	 * cluster containing this contig. After computing models, this is
	 * the final cluster this contig is assigned to.
	 *
	 * @return cluster containing this contig
	 */
	Cluster* cluster() const;

	/**
	 * @brief Setter for cluster containing this contig.
	 *
	 * @param cluster	cluster containing this contig
	 */
	void set_cluster(Cluster* cluster);

private:
	std::string id_; /**< Id. */
	int length_; /**< Length. */

	int modified_length_; /**< Modified length. */
	int left_edge_; /**< Starting point of the first window. */
	int right_edge_; /**< Ending point of the last window. */
	int num_windows_; /**< Number of windows. */

	int* sum_read_counts_; /**< Sum of read counts for all samples. */
	int** read_counts_; /**< Sum of read counts for all windows. */

	Cluster* cluster_; /**< Cluster containing this contig. */

	/**
	 * @brief Initializes read counts for all samples and windows to 0.
	 */
	void initReadCounts();
};


/** A map with contig ids as keys and pointers to corresponding objects as values. */
typedef std::unordered_map<std::string, Contig*> ContigMap;


/**
 * @brief Class containing internal contigs IO functions.
 * 
 * This class containes internal IO functions for saving and loading
 * extracted contig and mapping information.
 */
class ContigIO {
public:
	/**
	 * @brief Saves extracted contig and mapping information to a file.
	 *
	 * @param contigs					map with contig information
	 * @param sigma_contigs_file_path	path to file for saving contig information
	 */
	static void save_contigs(const ContigMap* contigs, const char* sigma_contigs_file_path);

	/**
	 * @brief Loads extracted contig and mapping information from a file.
	 *
	 * @param sigma_contigs_file_path	path to file for loading contig information
	 * @param contigs					map with contig information
	 */
	static void load_contigs(const char* sigma_contigs_file_path, ContigMap* contigs);
};




/**
* @brief Computes mode (most frequently occuring value) of a given vector.
*  
* @param Rs vector with values
* @return mode of Rs vector	
*/
double compute_mode(std::vector<double> &Rs);



/**
 * @brief Computes a global read count dispersion parameter R on contig windows.
 *
 * @param contigs	map with contig information
 * @return read count dispersion parameter on contig windows
 */
double compute_R(ContigMap* contigs);



#endif // CONTIG_H_
