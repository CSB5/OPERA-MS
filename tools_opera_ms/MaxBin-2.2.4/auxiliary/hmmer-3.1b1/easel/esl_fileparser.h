/* A simple token-based file parsing system.
 */
#ifndef eslFILEPARSER_INCLUDED
#define eslFILEPARSER_INCLUDED

#include <stdio.h>
#include "easel.h"

typedef struct {
  FILE *fp;			/* open file pointer, for reading                  */
  char *buf;			/* current line; will be modified by esl_strtok(). */
  int   buflen;			/* current allocated length of buf                 */
  char *s;			/* used by esl_strtok(); current position in buf.  */
  char  commentchar;		/* often '#'                                       */

  char *tok;			/* _NextLine() may remember a token...             */
  int   toklen;			/* ... and its length                              */
  char  tokchar;		/* ... and char that got overwritten by \0, if any */

  char *filename;		/* name of opened file; or NULL (if just a stream) */
  int   linenumber;		/* what line is loaded into buf; 1..nlines         */
  char  errbuf[eslERRBUFSIZE];  /* for holding error diagnostics                   */

  int   is_buffer;              /* the file has been buffered into memory          */
  char *mem_buffer;             /* pointer to the buffered file                    */
  int   mem_size;               /* size of the buffered file                       */
  int   mem_pos;                /* current position in the buffer                  */
} ESL_FILEPARSER;

extern int  esl_fileparser_Open(const char *filename, const char *envvar, ESL_FILEPARSER **ret_efp);
extern ESL_FILEPARSER *esl_fileparser_Create(FILE *fp);
extern ESL_FILEPARSER *esl_fileparser_CreateMapped(void *buffer, int size);
extern int  esl_fileparser_SetCommentChar  (ESL_FILEPARSER *efp, char c);
extern int  esl_fileparser_GetToken        (ESL_FILEPARSER *efp, char **opt_tok, int *opt_toklen);
extern int  esl_fileparser_NextLine        (ESL_FILEPARSER *efp);
extern int  esl_fileparser_NextLinePeeked  (ESL_FILEPARSER *efp, char *prefix, int plen);
extern int  esl_fileparser_GetTokenOnLine  (ESL_FILEPARSER *efp, char **opt_tok, int *opt_toklen);
extern int  esl_fileparser_GetRemainingLine(ESL_FILEPARSER *efp, char **ret_s);
extern void esl_fileparser_Destroy         (ESL_FILEPARSER *efp);
extern void esl_fileparser_Close           (ESL_FILEPARSER *efp);

#endif /*eslFILEPARSER_INCLUDED */
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.1b1; May 2013
 * Copyright (C) 2013 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
