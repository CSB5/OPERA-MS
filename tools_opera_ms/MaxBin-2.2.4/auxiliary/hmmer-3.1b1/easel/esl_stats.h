/* Foundation, miscellenea for the statistics modules.
 */
#ifndef eslSTATS_INCLUDED
#define eslSTATS_INCLUDED

/* 1. Summary statistics calculations */
extern int esl_stats_DMean(const double *x, int n, double *opt_mean, double *opt_var);
extern int esl_stats_FMean(const float  *x, int n, double *opt_mean, double *opt_var);
extern int esl_stats_IMean(const int    *x, int n, double *opt_mean, double *opt_var);

/* 2. Special functions */
extern int esl_stats_LogGamma(double x, double *ret_answer);
extern int esl_stats_Psi(double x, double *ret_answer);
extern int esl_stats_IncompleteGamma(double a, double x, double *ret_pax, double *ret_qax);

/* 3. Standard statistical tests */
extern int esl_stats_GTest(int ca, int na, int cb, int nb, double *ret_G, double *ret_P);
extern int esl_stats_ChiSquaredTest(int v, double x, double *ret_answer);

/* 4. Data fitting */
extern int esl_stats_LinearRegression(const double *x, const double *y, const double *sigma, int n,
				      double *opt_a,       double *opt_b,
				      double *opt_sigma_a, double *opt_sigma_b, double *opt_cov_ab,
				      double *opt_cc,      double *opt_Q);
#endif /*eslSTATS_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.1b1; May 2013
 * Copyright (C) 2013 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *
 * SVN $Id: esl_stats.h 771 2012-06-11 12:53:13Z eddys $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/easel/branches/hmmer/3.1/esl_stats.h $
 *****************************************************************/
