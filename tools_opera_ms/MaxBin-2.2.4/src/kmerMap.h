#ifndef __KMERMAP_H__
#define __KMERMAP_H__

#include "global_inc.h"
#include <math.h>

#ifndef WIN32
	#define _int64 long long
#endif

class kmerMap
{
	public:
		kmerMap(int inputlen, bool is_symmetric);
		~kmerMap();
		int getMapping(char *inputKmer);
		int getReverseMapping(char *inputKmer);
		int getReverseMapping(int index);
		int getEntryNum();
	private:
		bool err;
		bool symmetric;
		int len_mer;
		int entry_num;
		int kmer_num;
		int *kmer_small;
		int *kmer_large;
		int *kmer_map;
		void buildKmerMappingTable();
		unsigned int calNum(char *str);
		void revSeq(char *seq, char *result);
		void numToString(int num, char *str);
};

#endif

