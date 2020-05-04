#ifndef __MINMANHATTANDIST_H__
#define __MINMANHATTANDIST_H__

#include "AbstractDist.h"

class MinManhattanDist: public AbstractDist
{
	public:
		// Functions
		MinManhattanDist(int input_kmerlen): AbstractDist(input_kmerlen) {init();}
		double getDist(const char *inputseq1, const char *inputseq2);
		double getDist(double *input_pro1, double *input_pro2);
	private:
		void init();
};


#endif

