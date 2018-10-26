#ifndef EDGE_READER_H_
#define EDGE_READER_H_

#include "contig.h"
#include "edge.h"

/**
 * @brief An interface for edges file readers.
 *
 * Provides an interface for reading paired-end/mate-pair associations from
 * the output of an arbitrary edge bundling tool.
 */
class EdgeReader {
public:
	virtual ~EdgeReader(); /**< A virtual destructor. */

	/**
	 * @brief Reads edge information from given edges file.
	 *
	 * @param edges_file			path to edges file
	 * @param contigs				map with contig information
	 * @param edges					set of edges
	 * @param skipped_edges_file	path to skipped edges file
	 */
	virtual void read(const char* edges_file, const ContigMap* contigs, EdgeSet* edges, const char* skipped_edges_file) = 0;

	/**
	 * @brief Filters edges from given edges file.
	 *
	 * @param edges_file			path to edges file
	 * @param contigs				map with contig information
	 * @param filtered_edges_file	path to filtered edges file
	 */
	virtual void filter(const char* edges_file, const ContigMap* contigs, const char* filtered_edges_file, const char* skipped_edges_file) = 0;
};


/**
 * @brief Opera edges file reader.
 *
 * Enables reading edge information from the output of
 * <a href="http://sourceforge.net/projects/operasf/">Opera's</a>
 * bundling routine.
 */
class OperaBundleReader : public EdgeReader {
public:
	OperaBundleReader(); /**< An empty constructor. */

	/**
	 * @brief Reads edge information from Opera's bundle file.
	 *
	 * @copydetails EdgeReader::read(const char*, const ContigMap*, EdgeSet*, const char*)
	 */
	void read(const char* edges_file, const ContigMap* contigs, EdgeSet* edges, const char* skipped_edges_file);

	/**
	 * @brief Filters edges from Opera's bundle file.
	 *
	 * @copydetails EdgeReader::filter(const char*, const ContigMap*, const char*)
	 */
	void filter(const char* edges_file, const ContigMap* contigs, const char* filtered_edges_file, const char* skipped_edges_file);
};

#endif // EDGE_READER_H_
