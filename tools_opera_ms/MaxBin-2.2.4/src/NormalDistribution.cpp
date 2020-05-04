#include "NormalDistribution.h"

NormalDistribution::NormalDistribution(double input_mean, double input_std)
{
	pi = 3.1415926;
	mean = input_mean;
	std = input_std;
}

double NormalDistribution::prob(double input)
{
	double ret;

	ret = (1 / (std * sqrt(2 * pi))) * exp(-0.5 * pow(((input - mean) / std), 2));
	return(ret);
}
