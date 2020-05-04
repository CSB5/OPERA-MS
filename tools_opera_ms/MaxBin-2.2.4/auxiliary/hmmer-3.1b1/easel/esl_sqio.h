/* Unaligned sequence file i/o.
 * 
 * SVN $Id: esl_sqio.h 863 2013-04-18 13:40:02Z wheelert $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/easel/branches/hmmer/3.1/esl_sqio.h $
 */
#ifndef eslSQIO_INCLUDED
#define eslSQIO_INCLUDED

#include <stdio.h>

#include "esl_sqio_ascii.h"
#ifdef eslAUGMENT_NCBI
#include "esl_sqio_ncbi.h"
#endif
#include "esl_sq.h"

#ifdef eslAUGMENT_ALPHABET
#include "esl_alphabet.h"
#endif
#ifdef eslAUGMENT_MSA
#include "esl_msa.h"
#include "esl_msafile.h"
#endif

/*::cexcerpt::sq_sqio_data::begin::*/
/* ESL_SQDATA:
 * Data for different sequence formats.
 */
typedef union {
  ESL_SQASCII_DATA ascii;
#ifdef eslAUGMENT_NCBI
  ESL_SQNCBI_DATA  ncbi;
#endif
} ESL_SQDATA;
/*::cexcerpt::sq_sqio_data::end::*/

/* ESL_SQFILE:
 * An open sequence file for reading.
 */
typedef struct esl_sqio_s {
  char *filename;	      /* Name of file (for diagnostics)           */

  /* In digital mode, we have an alphabet ptr                             */
  int   do_digital;	      /* TRUE if we're reading in digital mode    */
#if defined(eslAUGMENT_ALPHABET)  
  const ESL_ALPHABET *abc;
#else
  void               *abc;
#endif

  /* Format-specific configuration                                        */
  int     format;	      /* Format code of this file                 */
  ESL_DSQ inmap[128];	      /* an input map, 0..127                     */

  /* function pointers to format specific routines                        */
  int   (*position)        (struct esl_sqio_s *sqfp, off_t offset);
  void  (*close)           (struct esl_sqio_s *sqfp);

  int   (*set_digital)     (struct esl_sqio_s *sqfp, const ESL_ALPHABET *abc);
  int   (*guess_alphabet)  (struct esl_sqio_s *sqfp, int *ret_type);

  int   (*read)            (struct esl_sqio_s *sqfp, ESL_SQ *sq);
  int   (*read_info)       (struct esl_sqio_s *sqfp, ESL_SQ *sq);
  int   (*read_seq)        (struct esl_sqio_s *sqfp, ESL_SQ *sq);
  int   (*read_window)     (struct esl_sqio_s *sqfp, int C, int W, ESL_SQ *sq);
  int   (*echo)            (struct esl_sqio_s *sqfp, const ESL_SQ *sq, FILE *ofp);

  int   (*read_block)      (struct esl_sqio_s *sqfp, ESL_SQ_BLOCK *sqBlock, int max_residues, int max_sequences, int long_target);

#ifdef eslAUGMENT_SSI
  int   (*open_ssi)        (struct esl_sqio_s *sqfp, const char *ssifile_hint);
  int   (*pos_by_key)      (struct esl_sqio_s *sqfp, const char *key);
  int   (*pos_by_number)   (struct esl_sqio_s *sqfp, int which);

  int   (*fetch)           (struct esl_sqio_s *sqfp, const char *key, ESL_SQ *sq);
  int   (*fetch_info)      (struct esl_sqio_s *sqfp, const char *key, ESL_SQ *sq);
  int   (*fetch_subseq)    (struct esl_sqio_s *sqfp, const char *source, int64_t start, int64_t end, ESL_SQ *sq);
#endif

  int   (*is_rewindable)   (const struct esl_sqio_s *sqfp);
  const char *(*get_error) (const struct esl_sqio_s *sqfp);

  ESL_SQDATA data;            /* format specific data                     */
} ESL_SQFILE;

#if defined(eslAUGMENT_ALPHABET)  
/* ESL_SQCACHE:
 * A entire database cached into memory.
 */
typedef struct esl_sqcache_s {
  char               *filename;    /* Name of file (for diagnostics)              */
  int                 format;      /* Format code of this file                    */

  const ESL_ALPHABET *abc;         /* alphabet for database                       */

  uint32_t            seq_count;   /* number of sequences                         */
  uint64_t            res_count;   /* number of residues                          */
  uint32_t            max_seq;     /* longest sequence                            */

  ESL_SQ             *sq_list;     /* list of cached sequences [0 .. seq_count-1] */

  void               *residue_mem; /* memory holding the residues                 */
  void               *header_mem;  /* memory holding the header strings           */

  uint64_t            res_size;    /* size of residue memory allocation           */
  uint64_t            hdr_size;    /* size of header memory allocation            */
} ESL_SQCACHE;
#endif

/*::cexcerpt::sq_sqio_format::begin::*/
/* Unaligned file format codes
 * These codes are coordinated with the msa module.
 *   - 0 is an unknown/unassigned format (eslSQFILE_UNKNOWN, eslMSAFILE_UNKNOWN)
 *   - <=100 is reserved for sqio, for unaligned formats
 *   - >100  is reserved for msa, for aligned formats
 */
#define eslSQFILE_UNKNOWN      0
#define eslSQFILE_FASTA        1
#define eslSQFILE_EMBL         2     /* EMBL/Swiss-Prot/TrEMBL */
#define eslSQFILE_GENBANK      3     /* GenBank */
#define eslSQFILE_DDBJ         4     /* DDBJ (currently passed to GenBank parser */
#define eslSQFILE_UNIPROT      5     /* UniProt (passed to EMBL parser) */
#define eslSQFILE_NCBI         6     /* NCBI (blast db) */
#define eslSQFILE_DAEMON       7     /* Farrar's "daemon" format for hmmpgmd queries: fasta with // end-of-record terminator */
#define eslSQFILE_HMMPGMD      8     /* Farrar's hmmpgmd database format: fasta w/ extra header line starting in '#' */
#define eslSQFILE_FMINDEX      9     /* the pressed FM-index format used for the FM-index-based MSV acceleration */
/*::cexcerpt::sq_sqio_format::end::*/


/* eslREADBUFSIZE is the fixed size of a block to bring in at one time,
 * in character-based (fread()) parsers (like the FASTA parser).
 */
#define eslREADBUFSIZE  4096

extern int  esl_sqfile_Open(const char *seqfile, int fmt, const char *env, ESL_SQFILE **ret_sqfp);
extern int  esl_sqfile_Position(ESL_SQFILE *sqfp, off_t offset);
extern void esl_sqfile_Close(ESL_SQFILE *sqfp);

#ifdef eslAUGMENT_ALPHABET
extern int  esl_sqfile_OpenDigital(const ESL_ALPHABET *abc, const char *filename, int format, const char *env, ESL_SQFILE **ret_sqfp);
extern int  esl_sqfile_SetDigital(ESL_SQFILE *sqfp, const ESL_ALPHABET *abc);
extern int  esl_sqfile_GuessAlphabet(ESL_SQFILE *sqfp, int *ret_type);

extern int  esl_sqfile_Cache(const ESL_ALPHABET *abc, const char *seqfile, int fmt, const char *env, ESL_SQCACHE **ret_sqcache);
extern void esl_sqfile_Free(ESL_SQCACHE *sqcache);
#endif

const char  *esl_sqfile_GetErrorBuf(const ESL_SQFILE *sqfp);
extern int   esl_sqfile_IsRewindable(const ESL_SQFILE *sqfp);
extern int   esl_sqio_Ignore(ESL_SQFILE *sqfp, const char *ignoredchars);
extern int   esl_sqio_AcceptAs(ESL_SQFILE *sqfp, char *xchars, char readas);
extern int   esl_sqio_EncodeFormat(char *fmtstring);
extern char *esl_sqio_DecodeFormat(int fmt);
extern int   esl_sqio_IsAlignment(int fmt);

extern int   esl_sqio_Read        (ESL_SQFILE *sqfp, ESL_SQ *sq);
extern int   esl_sqio_ReadInfo    (ESL_SQFILE *sqfp, ESL_SQ *sq);
extern int   esl_sqio_ReadWindow  (ESL_SQFILE *sqfp, int C, int W, ESL_SQ *sq);
extern int   esl_sqio_ReadSequence(ESL_SQFILE *sqfp, ESL_SQ *sq);
extern int   esl_sqio_Echo        (ESL_SQFILE *sqfp, const ESL_SQ *sq, FILE *ofp);
extern int   esl_sqio_ReadBlock   (ESL_SQFILE *sqfp, ESL_SQ_BLOCK *sqBlock, int max_residues, int max_sequences, int long_target);

#ifdef eslAUGMENT_SSI
extern int   esl_sqfile_OpenSSI         (ESL_SQFILE *sqfp, const char *ssifile_hint);
extern int   esl_sqfile_PositionByKey   (ESL_SQFILE *sqfp, const char *key);
extern int   esl_sqfile_PositionByNumber(ESL_SQFILE *sqfp, int which);

extern int   esl_sqio_Fetch      (ESL_SQFILE *sqfp, const char *key, ESL_SQ *sq);
extern int   esl_sqio_FetchInfo  (ESL_SQFILE *sqfp, const char *key, ESL_SQ *sq);
extern int   esl_sqio_FetchSubseq(ESL_SQFILE *sqfp, const char *source, int64_t start, int64_t end, ESL_SQ *sq);
#endif

extern int   esl_sqio_Write(FILE *fp, ESL_SQ *s, int format, int update);
extern int   esl_sqio_Parse(char *buffer, int size, ESL_SQ *s, int format);

#endif /*eslSQIO_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.1b1; May 2013
 * Copyright (C) 2013 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
