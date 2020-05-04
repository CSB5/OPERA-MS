/* Hyperexponential (mixture exponential) distributions.
 * 
 */
#ifndef eslHYPEREXP_INCLUDED
#define eslHYPEREXP_INCLUDED

#ifdef eslAUGMENT_RANDOM
#include <esl_random.h>
#endif
#ifdef eslAUGMENT_HISTOGRAM
#include <esl_histogram.h>
#endif
#ifdef eslAUGMENT_FILEPARSER
#include <esl_fileparser.h>
#endif

typedef struct {
  double *q;			/* mixture coefficients   [0..K-1]*/
  double *lambda;		/* scale params           [0..K-1]*/
  double *wrk;			/* tmp K-vector, for logpdf calc  */
  double  mu;			/* location (x offset) parameter  */
  int     K;			/* # of components                */
  char   *fixlambda;		/* TRUE to constrain a lambda val */
  int     fixmix;		/* TRUE to constrain the q's      */
} ESL_HYPEREXP;


extern ESL_HYPEREXP *esl_hyperexp_Create(int K);
extern void          esl_hyperexp_Destroy(ESL_HYPEREXP *h);
extern int           esl_hyperexp_Copy(ESL_HYPEREXP *src, ESL_HYPEREXP *dest);
extern int           esl_hyperexp_FixedUniformMixture(ESL_HYPEREXP *h);
extern int           esl_hyperexp_SortComponents(ESL_HYPEREXP *h);
extern int           esl_hyperexp_Write(FILE *fp, ESL_HYPEREXP *hxp);
extern int           esl_hyperexp_Dump(FILE *fp, ESL_HYPEREXP *hxp);
#ifdef eslAUGMENT_FILEPARSER
extern int           esl_hyperexp_Read(ESL_FILEPARSER *ef, ESL_HYPEREXP **ret_hxp);
extern int           esl_hyperexp_ReadFile(char *filename, ESL_HYPEREXP **ret_hxp);
#endif

extern double  esl_hxp_pdf    (double x, ESL_HYPEREXP *h);
extern double  esl_hxp_logpdf (double x, ESL_HYPEREXP *h);
extern double  esl_hxp_cdf    (double x, ESL_HYPEREXP *h);
extern double  esl_hxp_logcdf (double x, ESL_HYPEREXP *h);
extern double  esl_hxp_surv   (double x, ESL_HYPEREXP *h);
extern double  esl_hxp_logsurv(double x, ESL_HYPEREXP *h);
extern double  esl_hxp_invcdf (double p, ESL_HYPEREXP *h);

extern double  esl_hxp_generic_pdf   (double x, void *params);
extern double  esl_hxp_generic_cdf   (double x, void *params);
extern double  esl_hxp_generic_surv  (double x, void *params);
extern double  esl_hxp_generic_invcdf(double x, void *params);

extern int esl_hxp_Plot(FILE *fp, ESL_HYPEREXP *h,
			double (*func)(double x, ESL_HYPEREXP *h), 
			double xmin, double xmax, double xstep);

#ifdef eslAUGMENT_RANDOM
extern double esl_hxp_Sample(ESL_RANDOMNESS *r, ESL_HYPEREXP *h);
#endif
#ifdef eslAUGMENT_MINIMIZER
extern int esl_hxp_FitGuess   (double *x, int n, ESL_HYPEREXP *h);
extern int esl_hxp_FitComplete(double *x, int n, ESL_HYPEREXP *h);
#ifdef eslAUGMENT_HISTOGRAM
extern int esl_hxp_FitGuessBinned   (ESL_HISTOGRAM *g, ESL_HYPEREXP *h);
extern int esl_hxp_FitCompleteBinned(ESL_HISTOGRAM *g, ESL_HYPEREXP *h);
#endif
#endif

#endif /*eslHYPEREXP_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.1b1; May 2013
 * Copyright (C) 2013 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *
 * SVN $Id: esl_hyperexp.h 727 2011-10-24 17:17:32Z eddys $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/easel/branches/hmmer/3.1/esl_hyperexp.h $
 *****************************************************************/
