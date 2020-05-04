/* Calculates pairwise %id for all aligned sequence pairs in MSA
 */
#include "esl_config.h"

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "esl_distance.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"

static ESL_OPTIONS options[] = {
  /* name                type   default   env  range  togs  reqs incomp    help                                                   docgroup */
  { "-h",         eslARG_NONE,    FALSE, NULL,  NULL, NULL, NULL,  NULL,  "help; show brief info on version and usage",              1 },
  { "--informat",  eslARG_STRING,      NULL,  NULL, NULL,  NULL,  NULL, NULL, "specify the input MSA file is in format <s>", 0 }, 
  { "--outformat", eslARG_STRING, "Clustal",  NULL, NULL,  NULL,  NULL, NULL, "write the output MSA in format <s>",          0 }, 
  { "--dna",       eslARG_NONE,       FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use DNA alphabet",                            0 },
  { "--rna",       eslARG_NONE,       FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use RNA alphabet",                            0 },
  { "--amino",     eslARG_NONE,       FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use protein alphabet",                        0 },
 { 0,0,0,0,0,0,0,0,0,0 },
};

static char banner[] = "calculate pairwise %id for each seq pair in an MSA";
static char usage[]  = "[options] <msafile>";

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char         *msafile = esl_opt_GetArg(go, 1);
  ESL_ALPHABET *abc     = NULL;
  int           infmt   = eslMSAFILE_UNKNOWN;
  ESLX_MSAFILE *afp     = NULL;
  ESL_MSA      *msa     = NULL;
  FILE         *ofp     = stdout;
  int           nali    = 0;
  int           namewidth;
  double        pid;
  int           nid, n;
  int           i,j;
  int           status;

  /* allow user to assert the input MSA alphabet */
  if      (esl_opt_GetBoolean(go, "--rna"))   abc = esl_alphabet_Create(eslRNA);
  else if (esl_opt_GetBoolean(go, "--dna"))   abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, "--amino")) abc = esl_alphabet_Create(eslAMINO); 

  /* allow user to assert the input MSA format */
  if (esl_opt_IsOn(go, "--informat") &&
      (infmt = eslx_msafile_EncodeFormat(esl_opt_GetString(go, "--informat"))) == eslMSAFILE_UNKNOWN)
    esl_fatal("%s is not a valid MSA file format for --informat", esl_opt_GetString(go, "--informat"));

  /* digital open */
  if ( ( status = eslx_msafile_Open(&abc, msafile, NULL, infmt, NULL, &afp)) != eslOK)
    eslx_msafile_OpenFailure(afp, status);

  while ((status = eslx_msafile_Read(afp, &msa)) == eslOK)
    {	
      nali++;

      namewidth = esl_str_GetMaxWidth(msa->sqname, msa->nseq);

      for (i = 0; i < msa->nseq; i++)
	for (j = i+1; j < msa->nseq; j++)
	  {
	    esl_dst_XPairId(abc, msa->ax[i], msa->ax[j], &pid, &nid, &n);
	    fprintf(ofp, "%-*s %-*s %6.2f %6d %6d\n", namewidth, msa->sqname[i], namewidth, msa->sqname[j], pid*100.0, nid, n);
	  }

      esl_msa_Destroy(msa);
    }
  if (nali == 0 || status != eslEOF) eslx_msafile_ReadFailure(afp, status); 

  eslx_msafile_Close(afp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
  
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.1b1; May 2013
 * Copyright (C) 2013 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *
 * SVN $Id: esl-alipid.c 766 2012-06-04 13:11:51Z eddys $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/easel/branches/hmmer/3.1/miniapps/esl-alipid.c $
 *****************************************************************/
