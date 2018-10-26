#ifndef MAPPING_READER_H_
#define MAPPING_READER_H_

#include "contig.h"

/**
 * @brief An interface for mapping file readers.
 *
 * Provides an interface for reading read-to-contig mapping information from
 * the output of an arbitrary mapping tool.
 */
class MappingReader {
public:
	virtual ~MappingReader(); /**< A virtual destructor. */

	/**
	 * @brief Reads mapping information from given mapping file.
	 *
	 * @param mapping_file		path to mapping file
	 * @param sample_index		index of sequenced sample
	 * @param contigs			map with contig information
	 */
	virtual void read(const char* mapping_file, int sample_index, ContigMap* contigs) = 0;
};


/**
 * @brief SAM mapping file reader.
 *
 * Enables reading read-to-contig mapping information from
 * <a href="http://samtools.github.io/">SAM file format</a>.
 */
class SAMReader : public MappingReader {
public:
	SAMReader(); /**< An empty constructor.*/

	/**
	 * @brief Reads mapping information from SAM format.
	 *
	 * @copydetails MappingReader::read(const char*, int, ContigMap*)
	 */
	void read(const char* mapping_file, int sample_index, ContigMap* contigs);

private:
	/**
	 * @brief Reads mapping information from .sam file or stdin.
	 *
	 * @param mapping_fp		.sam file/stdin pointer
	 * @param sample_index		index of sequenced sample
	 * @param contigs			map with contig information
	 */
	void readInputStream(FILE* mapping_fp, int sample_index, ContigMap* contigs);
};

#endif // MAPPING_READER_H_