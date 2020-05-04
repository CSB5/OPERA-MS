#ifndef __NORMALDISTRIBUTION_H__
#define __NORMALDISTRIBUTION_H__

#include "global_inc.h"
#include <math.h>

class NormalDistribution
{
	public:
		NormalDistribution(double input_mean, double input_std);
		double prob(double input);
		//double cumu(double input);
	private:
		// Variables
		double mean;
		double std;
		double pi;
};




#endif