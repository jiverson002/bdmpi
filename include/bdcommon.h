/*!
\file
\brief The various data structures/constants shared by the master and the slaves
\date Started 4/3/13
\author George
*/

#ifndef _BDCOMMON_H_
#define _BDCOMMON_H_

#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <semaphore.h>
#include <mqueue.h>


/*************************************************************************/
/* Internal constant definitions  */
/*************************************************************************/
/* The size of the global shared memory region.
   TODO: This needs to be determined automatically, as it now
   limits the maximum number of running processes to roughly
   BDMPI_GLOBALSIZE/sizeof(int) */
#define BDMPI_GLOBALSMSIZE       12576

/* The initial size of the # of communicators array */
#define BDMPI_INIT_MAXNCOMM      2048

/* The maximum length of the complete working directory path */
#define BDMPI_WDIR_LEN           1024

/****************************************************************************/
/*!
 *  \details  Enable the use of sb_discard() throughout the BDMPI library.
 */
/****************************************************************************/
#define BDMPI_SB_DISCARD   1

/****************************************************************************/
/*!
 *  \details  Enable the use of sb_saveall() throughout the BDMPI library.
 */
/****************************************************************************/
#define BDMPI_SB_SAVEALL   2

/****************************************************************************/
/*!
 *  \details  Enable the "lazy-write" strategy in the sbmalloc library.  This
 *            means that memory allocations controlled by the sbmalloc library
 *            will not be written to disk until there is ``sufficient''
 *            pressure on the total DRAM to warrant such an action.  In this
 *            case, ``sufficient'' is determined by the resident memory
 *            command line parameter `-rm='.
 *
 *  \note     While compatible, it is not recommended to use this option with
 *            the #BDMPI_SB_SAVEALL option, since the latter will essentially
 *            negate the advantages of this strategy.
 */
/****************************************************************************/
#define BDMPI_SB_LAZYWRITE 4

/****************************************************************************/
/*!
 *  \details  Enable the "lazy-read" strategy in the sbmalloc library.  This
 *            means that memory allocations controlled by the sbmalloc library
 *            will not be read from disk and read protected until the
 *            application makes a read / write attempt to the memory location
 *            corresponding to the allocation.  Furthermore, rather than read
 *            the entire allocation chunk, the first time that any system page
 *            within it is accessed, memory is read and protected at a
 *            resolution of an sbpage, which can be any multiple of a system
 *            page.
 */
/****************************************************************************/
#define BDMPI_SB_LAZYREAD  8

/****************************************************************************/
/*!
 *  \details  Enable the "asynchronous I/O" strategy in the sbmalloc library.
 *            This means that when a memory request is made for a page, the
 *            page will be read if it is not already and returned immediately.
 *            Then, in the background, an `I/O thread' will continue reading
 *            the rest of the chunk that the page was from.
 */
/****************************************************************************/
#define BDMPI_SB_ASIO      16


/*************************************************************************/
/* Common macros */
/*************************************************************************/
#define S_IFSET(a,b)      IFSET(job->jdesc->dbglvl, (a), (b))
#define M_IFSET(a,b)      IFSET(job->dbglvl, (a), (b))
#define S_SB_IFSET(FLAG)  if ((FLAG) == (job->jdesc->sbopts&(FLAG)))
#define SB_SB_IFSET(FLAG) if ((FLAG) == (sbinfo->opts&(FLAG)))

#define BDWARN(expr)\
  do {\
    char hostname[9]; \
    if (!(expr)) {\
      gethostname(hostname, 9);\
      fprintf(stderr, "[%8s:%6d][WARN]Test failed on line %d of file %s.\n", \
         hostname, (int)getpid(), __LINE__, __FILE__);\
    }\
  } while(0)

#define BDASSERT(expr)\
  do {\
    char hostname[9]; \
    if (!(expr)) {\
      gethostname(hostname, 9);\
      fprintf(stderr, "[%8s:%6d] ASSERTION failed on line %d of file %s. [strerr: %s]\n", \
         hostname, (int)getpid(), __LINE__, __FILE__, strerror(errno));\
      abort();\
    }\
  } while(0)

#if 1
#define BD_GET_LOCK(lock)\
  {\
    int retval;\
    char hostname[9]; \
    if ((retval = pthread_mutex_lock(lock)) != 0) {\
      gethostname(hostname, 9);\
      fprintf(stderr, "[%8s:%5d] Error: Mutex lock failed on line %d of file %s. [retval: %d %s]\n", \
         hostname, (int)getpid(), __LINE__, __FILE__, retval, strerror(retval));\
      abort();\
    }\
  }
#else
#define BD_GET_LOCK(lock)\
  {\
    int retval;\
    char hostname[9];\
    struct timespec ts;\
    clock_gettime(CLOCK_REALTIME, &ts);\
    ts.tv_sec += 10;\
    if (ETIMEDOUT == (retval=pthread_mutex_timedlock(lock, &ts))) {\
      fprintf(stderr, "[%6ld:%s,%4d]: timed out waiting for mutex (%p)\n",\
        syscall(SYS_gettid), basename(__FILE__), __LINE__, (void*)(lock));\
      pthread_mutex_lock(lock);\
      fprintf(stderr, "[%6ld:%s,%4d]: locked mutex\n", syscall(SYS_gettid),\
        __FILE__, __LINE__);\
    }\
    else if (0 != retval) {\
      gethostname(hostname, 9);\
      fprintf(stderr, "[%8s:%5d] Error: Mutex lock failed on line %d of file %s. [retval: %d %s]\n", \
         hostname, (int)getpid(), __LINE__, __FILE__, retval, strerror(retval));\
      abort();\
    }\
    fprintf(stderr, "[%6ld:%s,%4d]: get lock (%p)\n",\
      syscall(SYS_gettid), basename(__FILE__), __LINE__,\
      (void*)(lock));\
  }
#endif

#define BD_LET_LOCK(lock)\
  {\
    int retval;\
    char hostname[9]; \
    if ((retval = pthread_mutex_unlock(lock)) != 0) {\
      gethostname(hostname, 9);\
      fprintf(stderr, "[%8s:%5d] Error. Mutex unlock failed on line %d of file %s. [retval: %d %s]\n", \
         hostname, (int)getpid(), __LINE__, __FILE__, retval, strerror(retval));\
      abort();\
    }\
    /*fprintf(stderr, "[%6ld:%s,%4d]: let lock (%p)\n",\
      syscall(SYS_gettid), basename(__FILE__), __LINE__,\
      (void*)(lock));*/\
  }

#define BD_TRY_LOCK(lock, haslock)\
  {\
    int retval;\
    char hostname[9];\
    if (0 == (retval=pthread_mutex_trylock(lock))) {\
      haslock = 1;\
      /*fprintf(stderr, "[%6ld:%s,%4d]: get lock (%p)\n",\
        syscall(SYS_gettid), basename(__FILE__), __LINE__,\
        (void*)(lock));*/\
    }\
    else if (EBUSY == retval) {\
      haslock = 0;\
    }\
    else {\
      gethostname(hostname, 9);\
      fprintf(stderr, "[%8s:%5d] Error. Mutex unlock failed on line %d of file %s. [retval: %d %s]\n", \
         hostname, (int)getpid(), __LINE__, __FILE__, retval, strerror(retval));\
      abort();\
    }\
  }

#if 1
#define BD_GET_SEM(sem)                                                   \
{                                                                         \
  char hostname[9];                                                       \
  if (-1 == sem_wait(sem)) {                                              \
    gethostname(hostname, 9);                                             \
    fprintf(stderr, "[%8s:%5d] Error: Semaphore wait failed on line %d "  \
      "of file %s. [errno: %d %s]\n", hostname, (int)getpid(), __LINE__,  \
      __FILE__, errno, strerror(errno));                                  \
    abort();                                                              \
  }                                                                       \
}
#else
#define BD_GET_SEM(sem)\
  {\
    char hostname[9];\
    struct timespec ts;\
    clock_gettime(CLOCK_REALTIME, &ts);\
    ts.tv_sec += 10;\
    if (0 != sem_timedwait(sem, &ts)) {\
      if (ETIMEDOUT == errno) {\
        fprintf(stderr, "[%6ld:%s,%4d]: timed out waiting for semahore (%p)\n",\
          syscall(SYS_gettid), basename(__FILE__), __LINE__, (void*)(sem));\
        sem_wait(sem);\
        fprintf(stderr, "[%6ld:%s,%4d]: get semahore (%p)\n",\
          syscall(SYS_gettid), basename(__FILE__), __LINE__, (void*)(sem));\
      }\
      else {\
        gethostname(hostname, 9);\
        fprintf(stderr, "[%8s:%5d] Error: Mutex lock failed on line %d of file %s. [errno: %d %s]\n", \
           hostname, (int)getpid(), __LINE__, __FILE__, errno, strerror(errno));\
        abort();\
      }\
    }\
    else {\
      fprintf(stderr, "[%6ld:%s,%4d]: get semahore (%p)\n",\
        syscall(SYS_gettid), basename(__FILE__), __LINE__,\
        (void*)(sem));\
    }\
  }
#endif

#define BD_LET_SEM(sem)                                                   \
{                                                                         \
  char hostname[9];                                                       \
  if (-1 == sem_post(sem)) {                                              \
    gethostname(hostname, 9);                                             \
    fprintf(stderr, "[%8s:%5d] Error: Semaphore post failed on line %d "  \
      "of file %s. [errno: %d %s]\n", hostname, (int)getpid(), __LINE__,  \
      __FILE__, errno, strerror(errno));                                  \
    abort();                                                              \
  }                                                                       \
  /*fprintf(stderr, "[%6ld:%s,%4d]: let semahore (%p)\n",\
    syscall(SYS_gettid), basename(__FILE__), __LINE__,\
    (void*)(sem));*/\
}

#define BD_TRY_SEM(SEM, BOOL)                                             \
do {                                                                      \
  char hostname[9];                                                       \
  if (0 == sem_trywait(SEM)) {                                            \
    (BOOL) = 1;                                                           \
    /*fprintf(stderr, "[%6ld:%s,%4d]: get semahore (%p)\n",\
      syscall(SYS_gettid), basename(__FILE__), __LINE__,\
      (void*)(SEM));*/\
  }                                                                       \
  else if (EAGAIN == errno) {                                             \
    (BOOL) = 0;                                                           \
  }                                                                       \
  else {                                                                  \
    gethostname(hostname, 9);                                             \
    fprintf(stderr, "[%8s:%5d] Error: Semaphore try lock failed on line " \
      "%d of file %s. [errno: %d %s]\n", hostname, (int)getpid(),         \
      __LINE__, __FILE__, errno, strerror(errno));                        \
    abort();                                                              \
  }                                                                       \
} while (0)

#define BD_COND_WAIT(cond, mtx)                                             \
{                                                                         \
  int retval;                                                             \
  char hostname[9];                                                       \
  if (0 != (retval=pthread_cond_wait(cond, mtx))) {                       \
    gethostname(hostname, 9);                                             \
    fprintf(stderr, "[%8s:%5d] Error: Condition wait failed on line %d "  \
      "of file %s. [retval: %d %s]\n", hostname, (int)getpid(), __LINE__, \
      __FILE__, retval, strerror(retval));                                \
    abort();                                                              \
  }                                                                       \
}

#define BD_GET_RDLOCK(lock)\
  do {\
    int i, retval;\
    char hostname[9]; \
    for (i=0; i<100; i++) {\
      if ((retval = pthread_rwlock_rdlock(lock)) == 0) break;\
      gethostname(hostname, 9);\
      fprintf(stderr, "[%8s:%5d] Error. RWlock rdlock failed on line %d of file %s. [retval: %d EAGAIN:%d EINVAL:%d EDEADLK:%d %s]\n", \
         hostname, (int)getpid(), __LINE__, __FILE__, retval, EAGAIN, EINVAL, EDEADLK, strerror(retval));\
    }\
    if (i==100) abort();\
  } while(0)

#define BD_GET_WRLOCK(lock)\
  do {\
    int i, retval;\
    char hostname[9]; \
    for (i=0; i<100; i++) {\
      if ((retval = pthread_rwlock_wrlock(lock)) == 0) break;\
      gethostname(hostname, 9);\
      fprintf(stderr, "[%8s:%5d] Error. RWlock wrlock failed on line %d of file %s. [retval: %d EAGAIN:%d EINVAL:%d EDEADLK:%d %s]\n", \
         hostname, (int)getpid(), __LINE__, __FILE__, retval, EAGAIN, EINVAL, EDEADLK, strerror(retval));\
    }\
    if (i==100) abort();\
  } while(0)

#define BD_LET_RDLOCK(lock)\
  {\
    int retval;\
    char hostname[9]; \
    if ((retval = pthread_rwlock_unlock(lock)) != 0) {\
      printf("***RWlock unlock failed on line %d of file %s. [retval: %d %s]\n", \
         __LINE__, __FILE__, retval, strerror(retval));\
      gethostname(hostname, 9);\
      fprintf(stderr, "[%8s:%5d] Error. RWlock unlock failed on line %d of file %s. [retval: %d %s]\n", \
         hostname, (int)getpid(), __LINE__, __FILE__, retval, strerror(retval));\
      abort();\
    }\
  }

#define BD_LET_WRLOCK(lock)\
  BD_LET_RDLOCK(lock)

#define BD_LET_RWLOCK(lock)\
  BD_LET_RDLOCK(lock)



/*************************************************************************/
/*! Debug Levels */
/*************************************************************************/
typedef enum {
  BDMPI_DBG_IPCM       = 1,    /*!< Shows detailed information about the master protocol */
  BDMPI_DBG_IPCS       = 2,    /*!< Shows detailed information about the slave protocol */
  BDMPI_DBG_SLVPOOL    = 4,    /*!< Shows detailed information about the slave pool */
  BDMPI_DBG_IPCM2      = 8,    /*!< Shows detailed information about the master protocol */
} mdbglvl_et;


/*************************************************************************/
/*! The job description information that will reside on the 'global'
    shared memory */
/*************************************************************************/
typedef struct {
  int nnodes;               /*!< The number of masters (i.e., MPI nodes) */
  int mynode;               /*!< The MPI rank of the master node */
  int np;                   /*!< The total number of slave processes over all nodes */
  int soffset;              /*!< The number of slaves in lower-ranked nodes */
  int ns;                   /*!< The number of slave processes for this node/master */
  int nr;                   /*!< The maximum number of slaves allowed to run */
  mdbglvl_et dbglvl;        /*!< The dbglvl of the execution */
  int sbopts;               /*!< The sb library options */
  size_t smsize;            /*!< The size of the shared memory comm buffer */
  size_t imsize;            /*!< The maximum size for in-memory buffering */
  size_t sbsize;            /*!< The minimum size for storage-backed allocation */
  size_t pgsize;            /*!< Number of system pages making up a single sb page */
  pid_t mpid;               /*!< The pid of the master */
  char wdir[BDMPI_WDIR_LEN]; /*!< The working directory */
} bdjdesc_t;


/*************************************************************************/
/*! The different types of messages between slaves and master */
/*************************************************************************/
typedef enum {
  BDMPI_MSGTYPE_INIT         =1,   /*!< a slave finished with init */
  BDMPI_MSGTYPE_FINALIZE     =2,   /*!< a slave finalize operation */
  BDMPI_MSGTYPE_BARRIER      =3,   /*!< a barrier operation */

  BDMPI_MSGTYPE_SEND         =4,   /*!< a data send operation */

  BDMPI_MSGTYPE_RECV         =50,  /*!< a data recv operation */
  BDMPI_MSGTYPE_IRECV        =51,  /*!< a data recv operation */
  BDMPI_MSGTYPE_PROBE        =52,  /*!< a probe operation */
  BDMPI_MSGTYPE_IPROBE       =53,  /*!< an iprobe operation */

  BDMPI_MSGTYPE_BCASTI       =60,  /*!< a broadcast init operation */
  BDMPI_MSGTYPE_BCASTF       =61,  /*!< a broadcast finish operation */
  BDMPI_MSGTYPE_ALLGATHERI   =62,  /*!< an all gather init operation */
  BDMPI_MSGTYPE_ALLGATHERF   =63,  /*!< an all gather finish operation */

  BDMPI_MSGTYPE_REDUCEI      =70,  /*!< a reduction init operation */
  BDMPI_MSGTYPE_REDUCEF      =71,  /*!< a reduction finish operation */
  BDMPI_MSGTYPE_ALLREDUCEI   =72,  /*!< an all reduction init operation */
  BDMPI_MSGTYPE_ALLREDUCEF   =73,  /*!< an all reduction finish operation */
  BDMPI_MSGTYPE_MERGEI       =74,  /*!< a merge init operation */
  BDMPI_MSGTYPE_MERGEF       =75,  /*!< a merge finish operation */

  BDMPI_MSGTYPE_GATHERI      =80,  /*!< a gather init operation */
  BDMPI_MSGTYPE_GATHERF      =81,  /*!< a gather finish operation */
  BDMPI_MSGTYPE_SCATTERI     =82,  /*!< a scatter init operation */
  BDMPI_MSGTYPE_SCATTERF     =83,  /*!< a scatter finish operation */
  BDMPI_MSGTYPE_ALLTOALLI    =84,  /*!< an all-to-all init operation */
  BDMPI_MSGTYPE_ALLTOALLF    =85,  /*!< an all-to-all finish operation */

  BDMPI_MSGTYPE_COMMDUP      =90, /*!< a comm_dup operation */
  BDMPI_MSGTYPE_COMMFREE     =91, /*!< a comm_free operation */
  BDMPI_MSGTYPE_COMMSPLIT    =92, /*!< a comm_split operation */

  BDMPI_MSGTYPE_CID          =100, /*!< a to master-node request for next mpi_commid */

  BDMPI_MSGTYPE_MEMLOAD      =200, /*!< a memory load operation */
  BDMPI_MSGTYPE_MEMSAVE      =201, /*!< a memory save operation */

  BDMPI_MSGTYPE_PROCEED      =210,  /*!< slave should proceed with execution */
  BDMPI_MSGTYPE_MEMFREE      =211,  /*!< slave should free its memory */

  BDMPI_MSGTYPE_NOOP         =999  /*!< a dummy message type */
} bdmsgtype_et;


/*************************************************************************/
/*! The message header communicated between slaves and master */
/*************************************************************************/
typedef struct {
  bdmsgtype_et msgtype;    /*!< The type of the message */
  int mcomm;               /*!< The master's ID of the communicator */
  int myrank;              /*!< The rank of the slave sending the msg */
  int copid;               /*!< The ID of collective operation */
  ssize_t fnum;            /*!< The globally-unique file number */
  int rmfile;              /*!< Indicates that any associated disk buffer can be rmed */
  int tag;                 /*!< The tag of the message */
  int source;              /*!< The rank of the source */
  int dest;                /*!< The rank of the destination */
  size_t count;            /*!< The length of the message in BDMPI_Datatype */
  int mpi_tag;             /*!< The mpi tag to be used for the data transfers */
  BDMPI_Datatype datatype;  /*!< The datatype of the message */
  BDMPI_Op op;              /*!< The reduction operation */
} bdmsg_t;


/*************************************************************************/
/*! Structures for the val-loc datatypes */
/*************************************************************************/
typedef struct {float val; int loc;} bdvlp_fi_t;
typedef struct {double val; int loc;} bdvlp_di_t;
typedef struct {long val; int loc;} bdvlp_li_t;
typedef struct {short val; int loc;} bdvlp_si_t;
typedef struct {int val; int loc;} bdvlp_ii_t;


/*************************************************************************/
/*! Prototypes shared across master and slave */
/*************************************************************************/
/* datatypes */
size_t bdmp_sizeof(BDMPI_Datatype datatype);
size_t bdmp_msize(size_t count, BDMPI_Datatype datatype);
int datatype_isvalid(BDMPI_Datatype datatype);

/* reduce.c */
void reduce_op(void *a, void *b, size_t count, BDMPI_Datatype datatype,
         BDMPI_Op op);
int op_isvalid(BDMPI_Op op);

/* scan.c */
void scan_init_op(void *a, size_t count, BDMPI_Datatype datatype,
         BDMPI_Op op);
void scan_op(void *a, void *b, void *r, int a2r, size_t count,
         BDMPI_Datatype datatype, BDMPI_Op op);

/* util.c */
int bdmp_madvise(char *ptr, size_t size, int advise);

/* xfer.c */
void xfer_setwdir(char *wdir);
ssize_t xfer_getfnum();
void xfer_unlink(ssize_t fnum);
void xfer_in_scb(bdscb_t *scb, void *buf, size_t count, BDMPI_Datatype datatype);
void xfer_out_scb(bdscb_t *scb, void *buf, size_t count, BDMPI_Datatype datatype);
void xfer_in_disk(ssize_t fnum, char *buf, size_t count, BDMPI_Datatype datatype, int rmfile);
void xfer_out_disk(ssize_t fnum, char *buf, size_t count, BDMPI_Datatype datatype);

/* debug.c */
void bdprintf(char *f_str,...);

#endif
