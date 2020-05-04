/* RNA secondary structure markup in WUSS notation.
 * 
 * SRE, Tue Feb 15 10:15:28 2005
 * SVN $Id: esl_wuss.h 791 2012-08-22 18:03:16Z rivase $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/easel/branches/hmmer/3.1/esl_wuss.h $
 */
#ifndef eslWUSS_INCLUDED
#define eslWUSS_INCLUDED

extern int esl_wuss2ct(char *ss, int len, int *ct);
extern int esl_ct2wuss(int *ct, int n, char *ss);
extern int esl_ct2simplewuss(int *ct, int n, char *ss);
extern int esl_wuss2kh(char *ss, char *kh);
extern int esl_kh2wuss(char *kh, char *ss);
extern int esl_wuss_full(char *oldss, char *newss);
extern int esl_wuss_nopseudo(char *ss1, char *ss2);


#endif /*eslWUSS_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.1b1; May 2013
 * Copyright (C) 2013 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
