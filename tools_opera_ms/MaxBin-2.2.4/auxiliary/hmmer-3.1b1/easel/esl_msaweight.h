/* Sequence weighting algorithms.
 * 
 * SRE, Sun Nov  5 09:11:13 2006 [Janelia]
 * SVN $Id: esl_msaweight.h 664 2011-02-27 17:08:36Z eddys $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/easel/branches/hmmer/3.1/esl_msaweight.h $
 */
#ifndef eslMSAWEIGHT_INCLUDED
#define eslMSAWEIGHT_INCLUDED

#include "esl_msa.h"

extern int esl_msaweight_GSC(ESL_MSA *msa);
extern int esl_msaweight_PB(ESL_MSA *msa);
extern int esl_msaweight_BLOSUM(ESL_MSA *msa, double maxid);
extern int esl_msaweight_IDFilter(const ESL_MSA *msa, double maxid, ESL_MSA **ret_newmsa);


#endif /*eslMSAWEIGHT_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.1b1; May 2013
 * Copyright (C) 2013 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
