/*
 * bdmprun.h
 *
 * This file contains the various header inclusions
 *
 * Started 3/31/13
 * George
 */

#define _BSD_SOURCE
#define _GNU_SOURCE

#include <GKlib.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <malloc.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <signal.h>
#include <setjmp.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/wait.h>
#include <unistd.h>
#include <mqueue.h>
#include <semaphore.h>
#include <pthread.h>
#include <mpi.h>


#if defined(ENABLE_OPENMP)
  #include <omp.h>
#endif

#include <bdmpi.h>
#include <bdipc.h>
#include <bdcommon.h>
#include "defs.h"
#include "macros.h"
#include "struct.h"
#include "proto.h"


#if defined(COMPILER_MSC)
#if defined(rint)
  #undef rint
#endif
#define rint(x) ((idx_t)((x)+0.5))  /* MSC does not have rint() function */
#define __func__ "dummy-function"
#endif


/* This is a global variable for accessing all information related to the master
   node's job. It is defined in main.c and as extern in all other files.  */
#if !defined(BDMPRUN_MAIN)
extern
#endif
mjob_t *job;

