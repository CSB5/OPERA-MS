/* Tracking cpu/system/elapsed time used by a process.
 * 
 * SRE, Wed Feb 22 19:30:36 2006 [St. Louis] [moved to Easel]
 * SRE, Thu Aug  3 08:00:35 2000 [St. Louis] [moved to SQUID]
 * SRE, Fri Nov 26 14:54:21 1999 [St. Louis] [HMMER]
 * SVN $Id: esl_stopwatch.h 664 2011-02-27 17:08:36Z eddys $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/easel/branches/hmmer/3.1/esl_stopwatch.h $
 */
#ifndef eslSTOPWATCH_INCLUDED
#define eslSTOPWATCH_INCLUDED

#include <time.h>
#ifdef HAVE_TIMES
#include <sys/times.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* need for sysconf() */
#endif

typedef struct {
  /* t0 and cpu0 keep base, when the watch was Start()'ed */
#ifdef HAVE_TIMES
  clock_t    t0;		/* Wall time, POSIX times()      */
  struct tms cpu0;		/* CPU/system time, POSIX times()*/
#else
  time_t  t0;			/* Wall time, fallback to ANSI time()  */
  clock_t cpu0;			/* CPU time, fallback to ANSI clock()  */
#endif

  /* elapsed/user/sys are t-t0 results for the last time the
   * watch was Stop()'ed.
   */
  double elapsed;               /* elapsed time, seconds */
  double user;                  /* CPU time, seconds     */
  double sys;                   /* system time, seconds  */
} ESL_STOPWATCH;


extern ESL_STOPWATCH *esl_stopwatch_Create(void);
extern void           esl_stopwatch_Destroy(ESL_STOPWATCH *w);

extern int esl_stopwatch_Start(ESL_STOPWATCH *w);
extern int esl_stopwatch_Stop(ESL_STOPWATCH *w);
extern int esl_stopwatch_Display(FILE *fp, ESL_STOPWATCH *w, char *prefix);

extern int esl_stopwatch_Include(ESL_STOPWATCH *master, ESL_STOPWATCH *w);


#endif /*eslSTOPWATCH_INCLUDED*/ 
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.1b1; May 2013
 * Copyright (C) 2013 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/


