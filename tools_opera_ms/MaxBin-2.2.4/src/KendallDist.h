#ifndef __KENDALLDIST_H__
#define __KENDALLDIST_H__

#include "AbstractDist.h"

class KendallDist: public AbstractDist
{
	public:
		// Functions
		KendallDist(int input_kmerlen): AbstractDist(input_kmerlen) {}
		double getDist(const char *inputseq1, const char *inputseq2);

};

#endif
