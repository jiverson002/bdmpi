/*
 * common.h
 *
 * This file contains the various header inclusions
 *
 * Started 8/28/13
 * George
 */

#define _BSD_SOURCE
#define _GNU_SOURCE
#define _DEFAULT_SOURCE

#include <GKlib.h>
#include <malloc.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <signal.h>
#include <setjmp.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <mqueue.h>
#include <semaphore.h>
#include <pthread.h>


#if defined(ENABLE_OPENMP)
  #include <omp.h>
#endif

#include <bdmpi.h>
#include <bdipc.h>
#include <bdcommon.h>
