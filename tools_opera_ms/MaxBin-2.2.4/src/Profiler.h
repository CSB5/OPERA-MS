#ifndef __PROFILER_H__
#define __PROFILER_H__

#include "kmerMap.h"

class Profiler
{
	public:
		Profiler(int input_kmerlen, const char *seq, kmerMap *input_kmap);
		~Profiler();
		double* getProfile();
		float getPercentN();
		float getTargetProfile(char *ch);
		void reset();
		void addProfile(Profiler *pf, long double weight);
		void calcProfile();
	private:
		// Variables
		int kmerlen;
		double *profile;
		char *tempKmer;
		float percent_N;
		kmerMap *kmap;
		long double total_weight;
		// Functions
		void computeProfile(const char *seq);
};



#endif
