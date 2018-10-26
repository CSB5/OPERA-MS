#include "Configure.h"

Configure::Configure(void)
{
}

Configure::~Configure(void)
{
	
}

// contig file related parameters
int Configure::FILE_FORMAT = FASTA;
int Configure::FILE_TYPE = VELVET;
bool Configure::FILTER_REPEAT = true;
bool Configure::KEEP_REPEAT_FULL = false;
double Configure::REPEAT_THRESHOLD = 1.5;
string Configure::CONTIG_FILE = "";
int Configure::CONTIG_SIZE_THERSHOLD = 500;

// mapping file related parameters
string Configure::MAP_FILE = "";
int Configure::MAP_TYPE = OPERA;//This variable is only to check if all the mapping type are edge. In this case the computation of the coverage and so repeat detection are disabled (as allcontigs will have a coverage of 1)
bool Configure::CALCULATE_LIB = true;
int Configure::LIB_MEAN = 10000;
int Configure::LIB_STD = 1000;
bool Configure::CALCULATE_ORI = true;
int Configure::READ_ORI = IN;

// bundle related parameters
int Configure::STD_TIMES = 6;
int Configure::CLUSTER_THRESHOLD = 5;

// result directory
string Configure::OUTPUT_FOLDER = "";
string Configure::DIRECTORY_SIGN = "/";  //"\\";

string Configure::SCAFFOLD_PREFEX = "opera_scaffold_";

bool Configure::ABORT = false;

// heuristic
bool Configure::HEURISTIC = false;
int Configure::STEP = 5;

// minimum valid distance for an edge
int Configure::MIN_DIS = 0;

// overlap graph realated variables
string Configure::OVERLAP_GRAPH_FILE = "";
int Configure::KMER = 49;

// multiple libraries variables
vector<LibInfo*> *Configure::MULTI_LIB_INFO = new vector<LibInfo*>;

//FAST
//int Configure::MAX_NUMBER_OF_PARTIAL_SCAFFOLDS = 100000;
//int Configure::PERCENTAGE_OF_INCREASING_THRESHOLD_GRAPHS = 5;
//SLOW
int Configure::MAX_NUMBER_OF_PARTIAL_SCAFFOLDS = 10000000;
int Configure::PERCENTAGE_OF_INCREASING_THRESHOLD_GRAPHS = 1;

int Configure::UNIT_OF_PARTIAL_SCAFFOLDS = 20000;

bool Configure::ONLY_BUNDLE = false;

bool Configure::USER_SPECIFY_KMER = false;

int Configure::INCREASED_DELTA = 5;

bool Configure::SIMULTANEOUSLY_HANDLE_LIBRARIES = true;

string Configure::REAL_POSITION_FILE = "";

string Configure::SAMDIR = "";

int Configure::UPPERBOUND = 0;

bool Configure::QUICK_MODE = false;

int Configure::TOP_VALUE = 10;

int Configure::SUPER_LIBRARY_THRESHOLD_1 = 1000;

int Configure::SUPER_LIBRARY_THRESHOLD_2 = 10000;

int Configure::SMALL_LIBRARY_SIZE = 1000;
	
int Configure::SUPER_SMALL_CONTIG_THRESHOLD = 500;

int Configure::MIN_LIBRARY_MEAN = 250;

int Configure::PLOIDY = 1;

int Configure::HAPLOID_COVERAGE = -1;

void Configure::destroy(){
	vector<LibInfo*>::iterator iter = MULTI_LIB_INFO->begin();
	while( iter != MULTI_LIB_INFO->end() )
	{
		delete *iter;
		iter = MULTI_LIB_INFO->erase( iter );
	}
	delete MULTI_LIB_INFO;
}