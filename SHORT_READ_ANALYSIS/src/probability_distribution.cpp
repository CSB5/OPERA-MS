#include <cstdlib>
#include <cmath>

#include "probability_distribution.h"

/** @{ */
/** Constant for Stirling's series approximation of log(x!). */
static const double LOG_SQRT2PI = 0.5 * log(2.0 * 3.14159265358979323846);
static const double LC1 = 1.0 / 12.0;
static const double LC2 = -1.0 / 360.0;
static const double LC3 = 1.0 / 1260.0;
static const double LC4 = -1.0 / 1680.0;
/** @} */

/**
 * @brief Computes Stirling's series approximation of log(x!).
 *
 * Computes <a href="http://en.wikipedia.org/wiki/Stirling%27s_approximation>">
 * Stirling's series approximation</a> of log(x!) using
 * log(x!) = log(x) + log((x-1)!) = log(x) + log(gamma(x)).
 *
 * @param x		value for which the approximation is computed
 * @return Stirling's series approximation of log(x!)
 */
static inline double stirling_log_factorial(double x) {
	const double r1 = 1.0 / x;
	const double r2 = r1 * r1;
	const double r3 = r1 * r2;
	const double r5 = r2 * r3;
	const double r7 = r2 * r5;

	return LC4 * r7 + LC3 * r5 + LC2 * r3 + LC1 * r1
			+ LOG_SQRT2PI + 0.5 * log(x) + x * (log(x) - 1.0);
}


ProbabilityDistribution::~ProbabilityDistribution() {}


PoissonDistribution::PoissonDistribution() {}

double PoissonDistribution::logpf(double mean, double dispersion, double value) const {
	const double lambda = round(mean);
	const double k = round(value);

	return k * log(lambda) - lambda - stirling_log_factorial(k);
}


NegativeBinomialDistribution::NegativeBinomialDistribution(double R) : R_(R){}
	/*log_p(log(1.0 - 1.0 / R)),
	log_1mp(log(1.0 / R)),
	xo1mx((1.0 / R) / (1.0 - 1.0 / R)) {}*/
	


double NegativeBinomialDistribution::logpf(double mean, double dispersion, double value) const {
	/*const double r = round(xo1mx * mean);
	const double k = round(value);

	return r * log_1mp + k * log_p + stirling_log_factorial(k + r - 1)
	- stirling_log_factorial(k) - stirling_log_factorial(r - 1);*/
	const double k = round(value);

	//TO DO
	//This computation can be splited to improve running time
	return dispersion * log(dispersion / (dispersion + mean)) + //Compute once
		lgamma(dispersion + k) + //Sum of the 2 sons
		k * log(mean / (dispersion + mean)) - //k * compute once
		lgamma(k + 1) -	//Sum of the 2 sons
		lgamma(dispersion); //Compute once globaly

	/*return R_ * log(R_ / (R_ + mean)) + stirling_log_factorial(R_ + k - 1) + 
		   k * log(mean / (R_ + mean)) - stirling_log_factorial(k) -
		   stirling_log_factorial(R_ - 1);*/
}
