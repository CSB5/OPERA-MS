/* Generate the 6 protein sequences from the reading frames of a nucleotide sequence
 *
 * WMA, Fri Mar 29 10:12 am 2013 [Janelia]
 */
#include "esl_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "esl_translate.h"
#include "easel.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             1 },
  { "-r",        eslARG_INT,    "15",     NULL, "n>0", NULL,  NULL, NULL, "minimum size of orf to report",                    1 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <file>";
static char banner[] = "translate reading frames of a sequence";

static void
cmdline_failure(char *argv0, char *format, ...)
{
  va_list argp;

  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  esl_usage(stdout, argv0, usage);
  printf("\nTo see more help on available options, do %s -h\n\n", argv0);
  exit(1);
}

static void
cmdline_help(char *argv0, ESL_GETOPTS *go) 
{
  esl_banner(stdout, argv0, banner);
  esl_usage (stdout, argv0, usage);
  puts("\n where general options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
  exit(0);
}

int main(int argc, char **argv)
{
  ESL_GETOPTS    *go       = NULL; 
  char           *filename = NULL;
  ESL_SQFILE   *sqfp       = NULL;
  ESL_SQ        *sq        = NULL;
  int             i;
  ESL_SQ        **prot;
  int             p_count;
  int min; 

  /* Parse command line */
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) cmdline_failure(argv[0], "Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) cmdline_failure(argv[0], "Error in app configuration: %s\n",   go->errbuf);
  if (esl_opt_GetBoolean(go, "-h") )                   cmdline_help(argv[0], go);
  if (esl_opt_ArgNumber(go) != 1)                      cmdline_failure(argv[0], "Incorrect number of command line arguments.\n");
  min = esl_opt_GetInteger(go, "-r");
    
  filename = esl_opt_GetArg(go, 1);

  if(esl_sqfile_Open(filename, eslSQFILE_FASTA, NULL, &sqfp) != eslOK) esl_fatal("could not open sequence file");
  if((sq = esl_sq_Create()) == NULL) esl_fatal("allocation failed");
  if(esl_sqio_Read(sqfp, sq) != eslOK) esl_fatal("could not read sequence file");
      
  if(esl_trans_orf(sq, &prot, &p_count, min) != eslOK) esl_fatal("translation failed");
      
  for(i = 0; i < p_count; i++)
  {
    esl_sqio_Write(stdout, prot[i], eslSQFILE_FASTA, 0);
    esl_sq_Destroy(prot[i]);
  }
  free(prot);
  
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
 *****************************************************************/ 

