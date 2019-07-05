#pragma once

#ifndef DEBUG
//#define DEBUG
#endif

#ifndef TIME
#define TIME
#endif

#ifndef LOG
#define LOG
#endif



#ifndef ANCHOR
//#define ANCHOR
#endif

#ifndef SPLIT
//#define SPLIT
#endif

#ifndef CHECK
//#define CHECK
#endif

#ifndef ORIENTATION
//#define ORIENTATION
#endif

#include <string>
#include <list>
#include <vector>
class LibInfo;
#include "LibInfo.h"

enum contig_file_format {FASTA, STATISTIC};
enum contig_file_type {VELVET, SOAP, SOAP_CONTIG, NORMAL, OPERA_SCAF};
enum map_type {BOWTIE, OPERA, SAM};
enum read_ori {IN, OUT, FORWARD1, FORWARD2, ERROR};
enum common_orientation {PLUS = 0, MINUS = 1, NONE, PLUSMINUS};
enum contig_position {START, END};
enum extend_orientation {RIGHT, LEFT, BOTH, NEITHER};
enum start_type {BORDER, ALL};
enum contig_type{UNIQUE, SUPER_SMALL, SMALL, BIG, REPEAT};
enum lib_type{SHORT, LONG};
enum library_size{SMALL_LIB, MEDIUM_LIB, BIG_LIB};

using namespace std;

class Configure
{
public:
	Configure(void);
	virtual ~Configure(void);

	//Attributes
	static const int TOO_MANY_PARTIAL_SCAFFOLDS = -2;
	static int FILE_FORMAT;				// fasta or statistic
	static int FILE_TYPE;				// velvet or soap
	static bool FILTER_REPEAT;			// if need to filter the repeat contigs, if false edges between repeats are removed
	static bool KEEP_REPEAT_FULL;			// if repeat need to be kept, if set to true no contigs are set to repeat
	static double REPEAT_THRESHOLD;		// repeat threshold, default = 1.5
	static string CONTIG_FILE;			// the contig file
	static int CONTIG_SIZE_THERSHOLD;	// the contig size threshold, default = 500

	static string MAP_FILE;				// the mapping file
	static int MAP_TYPE;				// the mapping type: bowtie
	static bool CALCULATE_LIB;			// if need to calculate the library information
	static int LIB_MEAN;				// the mean length of library, default = 10000
	static int LIB_STD;					// the standard deviation of library, default = 1000
	static bool CALCULATE_ORI;			// if need to calculate the orientation of paired end reads
	static int READ_ORI;				// the orientation of reads(in, out or forward)
	static string SAMDIR;				//samtools dir

	static int STD_TIMES;				// the times of stdandard deviation, default is 6
	static int CLUSTER_THRESHOLD;		// the paired end read cluster threshold
	static string OUTPUT_FOLDER;		// the result directory

	static string DIRECTORY_SIGN;		// the sign of directory, /(linux) or \(windows)

	static string SCAFFOLD_PREFEX;		// name of obtained scaffolds

	static int MIN_LIBRARY_MEAN;		// the minimum library size, default = 250

	static bool ABORT;					// if abort due to long time running

	// heuristic parameters
	static bool HEURISTIC;				// record if use heuristic method to speed up
	static int STEP;					// the max number of step to find unassigned nodes

	// threshold to consider an edge as invalid
	static int MIN_DIS;

	// overlap graph related variables
	static string OVERLAP_GRAPH_FILE;	// the file name of the overlap graph
	static int KMER;					// the size of the kmer

	// multiple libraries information
	static vector<LibInfo*> *MULTI_LIB_INFO;		// the information of multiple libraries

	// heuristic method
	static int MAX_NUMBER_OF_PARTIAL_SCAFFOLDS;		// the threshold of maximum visited number of partial scaffolds
	static int UNIT_OF_PARTIAL_SCAFFOLDS;                   // the unit of maximum visited number of partial scaffolds

	// only do bundling step
	static bool ONLY_BUNDLE;

	// record if user specified the value of kmer
	static bool USER_SPECIFY_KMER;

	// the increased value of cluster threshold
	static int INCREASED_DELTA;

	// the way to handle multiple libraries
	static bool SIMULTANEOUSLY_HANDLE_LIBRARIES;

	// the file containing the real position of unique contigs
	static string REAL_POSITION_FILE;

	// the upperbound of the library
	static int UPPERBOUND;

	// label if turn on quick mode (only use top N unassigned nodes
	static bool QUICK_MODE;
	static int TOP_VALUE;

	// int library threshold
	static int SUPER_LIBRARY_THRESHOLD_1;
	static int SUPER_LIBRARY_THRESHOLD_2;

	// ################# Statistic Usage ##################
	// the top N of unassigned nodes

	// small library threshold
	static int SMALL_LIBRARY_SIZE;
	
	// super small contig threshold
	static int SUPER_SMALL_CONTIG_THRESHOLD;

	// percentage of subgraphs increasing threshold
	static int PERCENTAGE_OF_INCREASING_THRESHOLD_GRAPHS;

	// haploid coverage, default value = -1
	static int HAPLOID_COVERAGE;

	// ploidy value, default value = 1
	static int PLOIDY;

	// function to destroy allocated MULTI_LIB_INFO
	static void destroy();
};
