#ifndef __SPEARMANDIST_H__
#define __SPEARMANDIST_H__

#include "global_inc.h"
#include "AbstractDist.h"
#include "kmerMap.h"
#include "quickSort.h"
#include "Profiler.h"

class SpearmanDist : public AbstractDist
{
	public:
		// Functions
		SpearmanDist(int input_kmerlen) : AbstractDist(input_kmerlen) {}
		double getDist(const char *inputseq1, const char *inputseq2);
		double getDist(double *input_pro1, double *input_pro2);
};

#endif
