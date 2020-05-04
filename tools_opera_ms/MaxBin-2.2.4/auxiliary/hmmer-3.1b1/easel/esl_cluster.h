/* Generalized single linkage clustering.
 * 
 * SRE, Mon Jan  7 09:40:06 2008 [Janelia]
 * SVN $Id: esl_cluster.h 664 2011-02-27 17:08:36Z eddys $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/easel/branches/hmmer/3.1/esl_cluster.h $
 */
#ifndef eslCLUSTER_INCLUDED
#define eslCLUSTER_INCLUDED

extern int esl_cluster_SingleLinkage(void *base, size_t n, size_t size, 
				     int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
				     int *workspace, int *assignments, int *ret_C);
#endif /*eslCLUSTER_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.1b1; May 2013
 * Copyright (C) 2013 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
