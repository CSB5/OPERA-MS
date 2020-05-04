/* Simple master/worker data parallelization using POSIX threads.
 * 
 * SVN $Id: esl_threads.h 753 2012-03-09 04:09:05Z eddys $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/easel/branches/hmmer/3.1/esl_threads.h $
 */
#ifndef eslTHREADS_INCLUDED
#define eslTHREADS_INCLUDED

#include <pthread.h>

typedef struct {
  int             threadCount;      /* number of active worker threads                           */
  pthread_t      *threadId;	    /* threadId for each worker thread; [0..threadCount-1]       */
  void          **data;		    /* data pointer for each worker thread; [0..threadCount-1]   */

  int             startThread;      /* number of worker threads currently blocked at start mutex */
  pthread_mutex_t startMutex;	    /* the starting gate                                         */
  pthread_cond_t  startCond;	    /* the signal that workers are synchronized and may start    */

  void           (*func)(void *);   /* each worker thread runs this function; arg is to data[]   */
} ESL_THREADS;


extern ESL_THREADS *esl_threads_Create(void (*func)(void *));
extern void         esl_threads_Destroy(ESL_THREADS *obj);

extern int esl_threads_AddThread     (ESL_THREADS *obj, void *data);
extern int esl_threads_GetWorkerCount(ESL_THREADS *obj);
extern int esl_threads_WaitForStart  (ESL_THREADS *obj);
extern int esl_threads_WaitForFinish (ESL_THREADS *obj);

extern int   esl_threads_Started (ESL_THREADS *obj, int *ret_workeridx);
extern void *esl_threads_GetData (ESL_THREADS *obj, int workeridx);
extern int   esl_threads_Finished(ESL_THREADS *obj, int workeridx);

extern int esl_threads_CPUCount(int *ret_ncpu);

#endif /*eslTHREADS_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.1b1; May 2013
 * Copyright (C) 2013 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
