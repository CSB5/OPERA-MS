#ifndef __ABSTRACTDIST_H__
#define __ABSTRACTDIST_H__

#include "global_inc.h"
#include "kmerMap.h"
#include "quickSort.h"
#include "Profiler.h"

class AbstractDist
{
	public:
		// Functions
		AbstractDist(int input_kmerlen);
		~AbstractDist();
		virtual double getDist(const char *inputseq1, const char *inputseq2) {return 0;};
		virtual double getDist(double *input_pro1, double *input_pro2) {return 0;};
		double getNpercentage();
		void setNormalization(bool input);
	protected:
		// Variables
		int kmerlen;
		int entrynum;
		Profiler *pro1;
		Profiler *pro2;
		double *profile1;
		double *profile2;
		int *rank1;
		int *rank2;
		int *tmprank;
		kmerMap *kmap;
		quickSort qs;
		double percent_N;
		bool normalize;
		// Functions
		void init();
		void resetDist();
		void computeRank(double *profile, int *rank);
		int computeDistFromRank();
};

#endif
