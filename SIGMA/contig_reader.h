#ifndef CONTIG_READER_H_
#define CONTIG_READER_H_

#include "contig.h"

/**
 * @brief An interface for contigs file readers.
 *
 * Provides an interface for reading contig information from the output of an
 * arbitrary assembly tool.
 */
class ContigReader {
public:
    virtual ~ContigReader(); /**< A virtual destructor. */

    /**
     * @brief Reads contig information from given contigs file
     *        and calculates total assembly length.
     *
     * @param contigs_file  path to contigs file
     * @param contigs       map with contig information
     * @return long int     total assembly size
     */
    virtual long int read(const char* contigs_file, ContigMap* contigs) = 0;

    /**
     * @brief Calculates total assembly length.
     *
     * @param contigs_file  path to contigs file
     * @param contigs       map with contig information
     * @return int          total assembly size
     */
    virtual void get_assembly_size(const char* contigs_file) = 0;
};

class AllReader : public ContigReader{
public:
    AllReader();
    
    long int read(const char* contigs_file, ContigMap* contigs);


    void get_assembly_size(const char* contigs_file);
};

/**
 * @brief SOAPdenovo contigs file reader.
 *
 * Enables reading contig information from the output of
 * <a href="http://soap.genomics.org.cn/soapdenovo.html">SOAPdenovo</a>
 * assembler.
 */

class SOAPdenovoReader : public ContigReader {
public:
    SOAPdenovoReader(); /**< An empty constructor. */

    /**
     * @brief Reads contig information from SOAPdenovo contigs file.
     *
     * @copydetails ContigReader::read(const char*, ContigMap*)
     */
    long int read(const char* contigs_file, ContigMap* contigs);

    /**
     * @brief Calculates total assemlby length from SOAPdenovo contigs file.
     *
     * @copydetails ContigReader::get_assembly_size(const char*)
     */
    void get_assembly_size(const char* contigs_file);
};

/**
 * @brief RAY contigs file reader.
 *
 * Enables reading contig information from the output of
 * <a href="http://to_add_corret_url.html">RAY</a>
 * assembler.
 */
class RAYReader : public ContigReader {
public:
    RAYReader(); /**< An empty constructor. */

    /**
     * @brief Reads contig information from RAY contigs file.
     *
     * @copydetails ContigReader::read(const char*, ContigMap*)
     */
    long int read(const char* contigs_file, ContigMap* contigs);

    /**
     * @brief Calculates total assemlby length from SOAPdenovo contigs file.
     *
     * @copydetails ContigReader::get_assembly_size(const char*)
     */
    void get_assembly_size(const char* contigs_file);
};



/**
 * @brief Velvet contigs file reader.
 *
 * Enables reading contig information from the output of
 * <a href="https://www.ebi.ac.uk/~zerbino/velvet/">Velvet</a>
 * assembler.
 */
class VelvetReader : public ContigReader {
public:
    VelvetReader(); /**< An empty constructor. */

    /**
     * @brief Reads contig information from Velvet contigs file.
     * 
     * @copydetails ContigReader::read(const char*, ContigMap*)
     */
    long int read(const char* contigs_file, ContigMap* contigs);

    /**
     * @brief Calculates total assemlby length from Velvet contigs file.
     *
     * @copydetails ContigReader::get_assembly_size(const char*)
     */
    void get_assembly_size(const char* contigs_file);
};

#endif // CONTIG_READER_H_endif // CONTIG_READER_H_
