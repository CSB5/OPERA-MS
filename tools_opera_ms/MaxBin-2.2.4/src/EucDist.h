#ifndef __EUCDIST_H__
#define __EUCDIST_H__

#include "SpearmanDist.h"

class EucDist: public AbstractDist
{
	public:
		// Functions
		EucDist(int input_kmerlen): AbstractDist(input_kmerlen) {}
		double getDist(const char *inputseq1, const char *inputseq2);
		double getDist(double *inputpro1, double *inputpro2);
};


#endif
