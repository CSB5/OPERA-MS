#ifndef __RATIOMANHATTANDIST_H__
#define __RATIOMANHATTANDIST_H__

#include "AbstractDist.h"

class RatioManhattanDist: public AbstractDist
{
	public:
		// Functions
		RatioManhattanDist(int input_kmerlen): AbstractDist(input_kmerlen) {init();}
		double getDist(const char *inputseq1, const char *inputseq2);
		double getDist(double *input_pro1, double *input_pro2);
	private:
		void init();
};


#endif
