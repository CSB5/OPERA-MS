/* Clustering sequences in an MSA by % identity.
 * 
 */
#ifndef eslMSACLUSTER_INCLUDED
#define eslMSACLUSTER_INCLUDED

#include "esl_msa.h"

extern int esl_msacluster_SingleLinkage(const ESL_MSA *msa, double maxid, 
					int **opt_c, int **opt_nin, int *opt_nc);

#endif /*eslMSACLUSTER_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.1b1; May 2013
 * Copyright (C) 2013 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *
 * SVN $Id: esl_msacluster.h 840 2012-12-20 18:16:14Z eddys $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/easel/branches/hmmer/3.1/esl_msacluster.h $
 *****************************************************************/
