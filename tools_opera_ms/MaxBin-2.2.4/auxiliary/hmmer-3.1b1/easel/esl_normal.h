/* Statistical routines for normal distributions
 * 
 * SRE, Tue Nov 21 14:29:02 2006 [Janelia]
 * SVN $Id: esl_normal.h 664 2011-02-27 17:08:36Z eddys $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/easel/branches/hmmer/3.1/esl_normal.h $
 */
#ifndef eslNORMAL_INCLUDED
#define eslNORMAL_INCLUDED

extern double esl_normal_pdf   (double x, double mu, double sigma);
extern double esl_normal_logpdf(double x, double mu, double sigma);
extern double esl_normal_cdf   (double x, double mu, double sigma);
extern double esl_normal_surv  (double x, double mu, double sigma);

#endif /*eslNORMAL_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.1b1; May 2013
 * Copyright (C) 2013 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
