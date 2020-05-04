/* PAML interface
 *
 *   "Phylogenetic Analysis by Maximum Likelihood"
 *   Ziheng Yang
 *   http://abacus.gene.ucl.ac.uk/software/paml.html
 *   [Yang97]
 * 
 *           incept: SRE, Tue Jul 13 13:20:08 2004 [St. Louis]
 * upgrade to Easel: SRE, Thu Mar  8 13:26:20 2007 [Janelia]
 * SVN $Id: esl_paml.h 664 2011-02-27 17:08:36Z eddys $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/easel/branches/hmmer/3.1/esl_paml.h $
 */
#ifndef eslPAML_INCLUDED
#define eslPAML_INCLUDED

#include <stdio.h>
#include <esl_dmatrix.h>

extern int esl_paml_ReadE(FILE *fp, ESL_DMATRIX *E, double *pi);


#endif /*eslPAML_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.1b1; May 2013
 * Copyright (C) 2013 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
