#ifndef __MANHATTANDIST_H__
#define __MANHATTANDIST_H__

#include "AbstractDist.h"

class ManhattanDist: public AbstractDist
{
	public:
		// Functions
		ManhattanDist(int input_kmerlen): AbstractDist(input_kmerlen) {}
		double getDist(const char *inputseq1, const char *inputseq2);
		double getDist(double *input_pro1, double *input_pro2);
};


#endif
