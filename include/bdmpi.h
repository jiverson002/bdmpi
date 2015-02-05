/*!
\file bdmpi.h
\brief This file contains function prototypes and constant definitions for BDMPI
 *
\author George
\date   Started 3/31/13
\version\verbatim $Id$\endverbatim
*/

#ifndef _BDMPI_H_
#define _BDMPI_H_

#ifndef _GNU_SOURCE
# define _GNU_SOURCE
#endif
#include <unistd.h>
#include <sys/syscall.h>
#include <sys/types.h>


/* Uniform definitions for various compilers */
#if defined(_MSC_VER)
  #define COMPILER_MSC
#endif
#if defined(__ICC)
  #define COMPILER_ICC
#endif
#if defined(__GNUC__)
  #define COMPILER_GCC
#endif

/* Include c99 int definitions and need constants. When building the library,
 * these are already defined by GKlib; hence the test for _GKLIB_H_ */
#ifndef _GKLIB_H_
#ifdef COMPILER_MSC
#include <limits.h>

typedef __int32 int32_t;
typedef __int64 int64_t;
#define PRId32       "I32d"
#define PRId64       "I64d"
#define SCNd32       "ld"
#define SCNd64       "I64d"
#define INT32_MIN    ((int32_t)_I32_MIN)
#define INT32_MAX    _I32_MAX
#define INT64_MIN    ((int64_t)_I64_MIN)
#define INT64_MAX    _I64_MAX
#else
#include <inttypes.h>
#endif
#endif



/* BDMPI's version number */
#define BDMPI_VER_MAJOR         0
#define BDMPI_VER_MINOR         1
#define BDMPI_VER_SUBMINOR      0


/* Wild-cards */
#define BDMPI_ANY_SOURCE         ((int)-1)
#define BDMPI_ANY_TAG            ((int)-1)

#define BDMPI_PROC_NULL          ((int)-2)


/* Return codes  */
#define BDMPI_SUCCESS        ((int)1001)  /*!< Returned normally */
#define BDMPI_INPROGRESS     ((int)1002)  /*!< Indicates an active async operation */
#define BDMPI_ERR_UNKNOWN    ((int)1003)  /*!< A generic error */
#define BDMPI_ERR_COMM       ((int)1004)  /*!< Invalid communicator */
#define BDMPI_ERR_ARG        ((int)1005)  /*!< Invalid argument */
#define BDMPI_ERR_BUFFER     ((int)1006)  /*!< Invalid buffer */
#define BDMPI_ERR_COUNT      ((int)1007)  /*!< Invalid count */
#define BDMPI_ERR_TYPE       ((int)1008)  /*!< Invalid type */
#define BDMPI_ERR_TAG        ((int)1009)  /*!< Invalid tag */
#define BDMPI_ERR_ROOT       ((int)1010)  /*!< Invalid root */
#define BDMPI_ERR_RANK       ((int)1011)  /*!< Invalid rank */
#define BDMPI_ERR_TRUNCATE   ((int)1012)  /*!< Message truncated on receive */
#define BDMPI_ERR_OP         ((int)1013)  /*!< Invalid reduction operation */



/* Reduction operation types */
typedef int BDMPI_Op;

#define BDMPI_OP_NULL   ((BDMPI_Op) -1)

#define BDMPI_MAX       ((BDMPI_Op)1001)    /*!< max reduction */
#define BDMPI_MIN       ((BDMPI_Op)1002)    /*!< min reduction */
#define BDMPI_SUM       ((BDMPI_Op)1003)    /*!< sum reduction */
#define BDMPI_PROD      ((BDMPI_Op)1004)    /*!< prod reduction */
#define BDMPI_LAND      ((BDMPI_Op)1005)    /*!< logical and reduction */
#define BDMPI_BAND      ((BDMPI_Op)1006)    /*!< boolean and reduction */
#define BDMPI_LOR       ((BDMPI_Op)1007)    /*!< logical or reduction */
#define BDMPI_BOR       ((BDMPI_Op)1008)    /*!< boolean or reduction */
#define BDMPI_LXOR      ((BDMPI_Op)1009)    /*!< logical xor reduction */
#define BDMPI_BXOR      ((BDMPI_Op)1010)    /*!< boolean xor reduction */
#define BDMPI_MAXLOC    ((BDMPI_Op)1011)    /*!< max value and location reduction */
#define BDMPI_MINLOC    ((BDMPI_Op)1012)    /*!< min value and location reduction */



/* Data types */
typedef int BDMPI_Datatype;

#define BDMPI_DATATYPE_NULL     ((BDMPI_Datatype)-1)

#define BDMPI_CHAR                ((BDMPI_Datatype)1001)
#define BDMPI_SIGNED_CHAR         ((BDMPI_Datatype)1002)
#define BDMPI_UNSIGNED_CHAR       ((BDMPI_Datatype)1003)
#define BDMPI_BYTE                ((BDMPI_Datatype)1004)
#define BDMPI_WCHAR               ((BDMPI_Datatype)1005)
#define BDMPI_SHORT               ((BDMPI_Datatype)1006)
#define BDMPI_UNSIGNED_SHORT      ((BDMPI_Datatype)1007)
#define BDMPI_INT                 ((BDMPI_Datatype)1008)
#define BDMPI_UNSIGNED            ((BDMPI_Datatype)1009)
#define BDMPI_LONG                ((BDMPI_Datatype)1010)
#define BDMPI_UNSIGNED_LONG       ((BDMPI_Datatype)1011)
#define BDMPI_LONG_LONG_INT       ((BDMPI_Datatype)1012)
#define BDMPI_UNSIGNED_LONG_LONG  ((BDMPI_Datatype)1013)
#define BDMPI_INT8_T              ((BDMPI_Datatype)1014)
#define BDMPI_UINT8_T             ((BDMPI_Datatype)1015)
#define BDMPI_INT16_T             ((BDMPI_Datatype)1016)
#define BDMPI_UINT16_T            ((BDMPI_Datatype)1017)
#define BDMPI_INT32_T             ((BDMPI_Datatype)1018)
#define BDMPI_UINT32_T            ((BDMPI_Datatype)1019)
#define BDMPI_INT64_T             ((BDMPI_Datatype)1020)
#define BDMPI_UINT64_T            ((BDMPI_Datatype)1021)
#define BDMPI_SIZE_T              ((BDMPI_Datatype)1022)
#define BDMPI_SSIZE_T             ((BDMPI_Datatype)1023)
#define BDMPI_FLOAT               ((BDMPI_Datatype)1024)
#define BDMPI_DOUBLE              ((BDMPI_Datatype)1025)
#define BDMPI_FLOAT_INT           ((BDMPI_Datatype)1026)
#define BDMPI_DOUBLE_INT          ((BDMPI_Datatype)1027)
#define BDMPI_LONG_INT            ((BDMPI_Datatype)1028)
#define BDMPI_SHORT_INT           ((BDMPI_Datatype)1029)
#define BDMPI_2INT                ((BDMPI_Datatype)1030)


/*! Communicators */
typedef struct {
  int mcomm;   /*!< The master's ID of the communicator */
  int size;    /*!< The total number of processes in the communicator */
  int rank;    /*!< The global rank of the process in the communicator */
  int lsize;   /*!< The number of intra processes in the communicator */
  int lrank;   /*!< The intra-rank of the slave process */
  int nsize;   /*!< The number of master nodes in the communicator */
  int nrank;   /*!< The node rank of the slaves */
  int rrank;   /*!< The global rank of the lrank==0 slave (i.e., the local "root") */
  int copid;   /*!< A counter forming the ID of each successive collective op */
} bdscomm_t;

typedef bdscomm_t *BDMPI_Comm;

#define BDMPI_COMM_NULL          NULL

#ifdef BDMPLIB_INIT_C
BDMPI_Comm BDMPI_COMM_WORLD  = BDMPI_COMM_NULL;
BDMPI_Comm BDMPI_COMM_SELF   = BDMPI_COMM_NULL;
BDMPI_Comm BDMPI_COMM_CWORLD = BDMPI_COMM_NULL;
BDMPI_Comm BDMPI_COMM_NODE   = BDMPI_COMM_NULL;
#else
extern BDMPI_Comm BDMPI_COMM_WORLD;
extern BDMPI_Comm BDMPI_COMM_SELF;
extern BDMPI_Comm BDMPI_COMM_CWORLD;
extern BDMPI_Comm BDMPI_COMM_NODE;
#endif


/*! Status */
typedef struct {
  int MPI_SOURCE;
  int MPI_TAG;
  int MPI_ERROR;
  int BDMPI_SOURCE;
  int BDMPI_TAG;
  int BDMPI_ERROR;
  int dest;
  BDMPI_Datatype datatype;
  size_t count;
  BDMPI_Comm comm;
} BDMPI_Status;

#define BDMPI_STATUS_IGNORE ((BDMPI_Status *)1)


/*! Requests types */
typedef enum {
  BDMPI_REQUEST_ISEND,           /*!< Request created due to an Isend */
  BDMPI_REQUEST_IRECV,           /*!< Request created due to an Irecv */
  BDMPI_REQUEST_REDUCEI,         /*!< Request created due to a Reduce_init */
} bdmp_request_et;


/*! Request */
typedef struct {
  bdmp_request_et type;
  int state;
  void *buf;
  BDMPI_Status status;
  int copid;
} bdrequest_t;

typedef bdrequest_t * BDMPI_Request;

#define BDMPI_REQUEST_NULL NULL


/*------------------------------------------------------------------------
* Prototypes
*-------------------------------------------------------------------------*/
#ifdef _WINDLL
#define BDMPI_API(type) __declspec(dllexport) type __cdecl
#elif defined(__cdecl)
#define BDMPI_API(type) type __cdecl
#else
#define BDMPI_API(type) type
#endif


#ifdef __cplusplus
extern "C" {
#endif

BDMPI_API(int)    BDMPI_Print_vm_info(char *hdr);
BDMPI_API(int)    BDMPI_mlockall(int flags);
BDMPI_API(int)    BDMPI_munlockall(void);
BDMPI_API(int)    BDMPI_mlock(const void *addr, size_t len);
BDMPI_API(int)    BDMPI_munlock(const void *addr, size_t len);
BDMPI_API(int)    BDMPI_Entercritical(void);
BDMPI_API(int)    BDMPI_Exitcritical(void);

BDMPI_API(void *) BDMPI_sbmalloc(size_t size);
BDMPI_API(void *) BDMPI_sbrealloc(void *oldptr, size_t size);
BDMPI_API(void)   BDMPI_sbfree(void *ptr);
BDMPI_API(void)   BDMPI_sbload(void *ptr);
BDMPI_API(void)   BDMPI_sbloadall(void);
BDMPI_API(void)   BDMPI_sbunload(void *ptr);
BDMPI_API(void)   BDMPI_sbunloadall(void);
BDMPI_API(void)   BDMPI_sbdiscard(void *ptr, ssize_t size);

BDMPI_API(int) BDMPI_Init(int *argc, char **argv[]);
BDMPI_API(int) BDMPI_Finalize();
BDMPI_API(int) BDMPI_Comm_size(BDMPI_Comm comm, int *size);
BDMPI_API(int) BDMPI_Comm_rank(BDMPI_Comm comm, int *rank);
BDMPI_API(int) BDMPI_Comm_lsize(BDMPI_Comm comm, int *lsize);
BDMPI_API(int) BDMPI_Comm_lrank(BDMPI_Comm comm, int *lrank);
BDMPI_API(int) BDMPI_Comm_nsize(BDMPI_Comm comm, int *nsize);
BDMPI_API(int) BDMPI_Comm_nrank(BDMPI_Comm comm, int *nrank);
BDMPI_API(int) BDMPI_Comm_rrank(BDMPI_Comm comm, int *rrank);
BDMPI_API(int) BDMPI_Comm_dup(BDMPI_Comm comm, BDMPI_Comm *newcomm);
BDMPI_API(int) BDMPI_Comm_free(BDMPI_Comm *comm);
BDMPI_API(int) BDMPI_Comm_split(BDMPI_Comm comm, int color, int key, BDMPI_Comm *newcomm);
BDMPI_API(int) BDMPI_Send(void *buf, size_t count, BDMPI_Datatype datatype,
                   int dest, int tag, BDMPI_Comm comm);
BDMPI_API(int) BDMPI_Isend(void *buf, size_t count, BDMPI_Datatype datatype,
                   int dest, int tag, BDMPI_Comm comm, BDMPI_Request *request);
BDMPI_API(int) BDMPI_Recv(void *buf, size_t count, BDMPI_Datatype datatype,
                   int source, int tag, BDMPI_Comm comm, BDMPI_Status *status);
BDMPI_API(int) BDMPI_Irecv(void *buf, size_t count, BDMPI_Datatype datatype,
                   int source, int tag, BDMPI_Comm comm, BDMPI_Request *request);
BDMPI_API(int) BDMPI_Sendrecv(void *sendbuf, size_t sendcount, BDMPI_Datatype sendtype,
                   int dest, int sendtag, void *recvbuf, size_t recvcount,
                   BDMPI_Datatype recvtype, int source, int recvtag, BDMPI_Comm comm,
                   BDMPI_Status *status);
BDMPI_API(int) BDMPI_Test(BDMPI_Request *request, int *flag, BDMPI_Status *status);
BDMPI_API(int) BDMPI_Wait(BDMPI_Request *request, BDMPI_Status *status);
BDMPI_API(int) BDMPI_Waitall(int count, BDMPI_Request *request, BDMPI_Status *status);
BDMPI_API(int) BDMPI_Get_count(BDMPI_Status *status, BDMPI_Datatype datatype,
                   size_t *count);
BDMPI_API(int) BDMPI_Barrier(BDMPI_Comm comm);
BDMPI_API(int) BDMPI_Bcast(void *buf, size_t count, BDMPI_Datatype datatype,
                   int root, BDMPI_Comm comm);
BDMPI_API(int) BDMPI_Allgather(void *sendbuf, size_t sendcount, BDMPI_Datatype sendtype,
                   void *recvbuf, size_t recvcount, BDMPI_Datatype recvtype,
                   BDMPI_Comm comm);
BDMPI_API(int) BDMPI_Allgatherv(void *sendbuf, size_t sendcount, BDMPI_Datatype sendtype,
                   void *recvbuf, size_t *recvcounts, size_t *displs,
                   BDMPI_Datatype recvtype, BDMPI_Comm comm);
BDMPI_API(int) BDMPI_Gather(void *sendbuf, size_t sendcount, BDMPI_Datatype sendtype,
                   void *recvbuf, size_t recvcount, BDMPI_Datatype recvtype,
                   int root, BDMPI_Comm comm);
BDMPI_API(int) BDMPI_Gatherv(void *sendbuf, size_t sendcount, BDMPI_Datatype sendtype,
                   void *recvbuf, size_t *recvcounts, size_t *displs,
                   BDMPI_Datatype recvtype, int root, BDMPI_Comm comm);
BDMPI_API(int) BDMPI_Scatter(void *sendbuf, size_t sendcount, BDMPI_Datatype sendtype,
                   void *recvbuf, size_t recvcount, BDMPI_Datatype recvtype,
                   int root, BDMPI_Comm comm);
BDMPI_API(int) BDMPI_Scatterv(void *sendbuf, size_t *sendcounts, size_t *displs,
                   BDMPI_Datatype sendtype, void *recvbuf, size_t recvcount,
                   BDMPI_Datatype recvtype, int root, BDMPI_Comm comm);
BDMPI_API(int) BDMPI_Alltoall(void *sendbuf, size_t sendcount, BDMPI_Datatype sendtype,
                   void *recvbuf, size_t recvcount, BDMPI_Datatype recvtype,
                   BDMPI_Comm comm);
BDMPI_API(int) BDMPI_Alltoallv(void *sendbuf, size_t *sendcounts, size_t *sdipls,
                   BDMPI_Datatype sendtype, void *recvbuf, size_t *recvcounts,
                   size_t *rdispls, BDMPI_Datatype recvtype, BDMPI_Comm comm);
BDMPI_API(int) BDMPI_Reduce(void *sendbuf, void *recvbuf, size_t count,
                   BDMPI_Datatype datatype, BDMPI_Op op, int root, BDMPI_Comm comm);
BDMPI_API(int) BDMPI_Reduce_init(void *sendbuf, void *recvbuf, size_t count,
                   BDMPI_Datatype datatype, BDMPI_Op op, int root, BDMPI_Comm comm,
                   BDMPI_Request *request);
BDMPI_API(int) BDMPI_Reduce_fine(void *recvbuf, size_t count, BDMPI_Datatype datatype,
                   BDMPI_Op op, int root, BDMPI_Comm comm, BDMPI_Request *request);
BDMPI_API(int) BDMPI_Allreduce(void *sendbuf, void *recvbuf, size_t count,
                   BDMPI_Datatype datatype, BDMPI_Op op, BDMPI_Comm comm);
BDMPI_API(int) BDMPI_Merge(void *sendbuf, int *sendids, size_t sendcount,
                   void *recvbuf, int *recvids, size_t *r_recvcount,
                   BDMPI_Datatype datatype, BDMPI_Op op, int root, BDMPI_Comm comm);
BDMPI_API(int) BDMPI_Probe(int source, int tag, BDMPI_Comm comm, BDMPI_Status *status);
BDMPI_API(int) BDMPI_Iprobe(int source, int tag, BDMPI_Comm comm, int *flag,
                   BDMPI_Status *status);
BDMPI_API(int) BDMPI_Scan(void *sendbuf, void *recvbuf, size_t count,
                   BDMPI_Datatype datatype, BDMPI_Op op, BDMPI_Comm comm);
BDMPI_API(int) BDMPI_Exscan(void *sendbuf, void *recvbuf, size_t count,
                   BDMPI_Datatype datatype, BDMPI_Op op, BDMPI_Comm comm);
BDMPI_API(double) BDMPI_Wtime(void);

#ifdef __cplusplus
}
#endif


#ifndef BDMPRUN
/*------------------------------------------------------------------------
* The MPI-compliant interface
*-------------------------------------------------------------------------*/
#define MPI_UNDEFINED    (-32766)

/* Wild-cards */
#define MPI_ANY_SOURCE  BDMPI_ANY_SOURCE
#define MPI_ANY_TAG     BDMPI_ANY_TAG


/* Return codes  */
#define MPI_SUCCESS       BDMPI_SUCCESS
#define MPI_INPROGRESS    BDMPI_INPROGRESS
#define MPI_ERR_UNKNOWN   BDMPI_ERR_UNKNOWN
#define MPI_ERR_COMM      BDMPI_ERR_COMM
#define MPI_ERR_ARG       BDMPI_ERR_ARG
#define MPI_ERR_BUFFER    BDMPI_ERR_BUFFER
#define MPI_ERR_COUNT     BDMPI_ERR_COUNT
#define MPI_ERR_TYPE      BDMPI_ERR_TYPE
#define MPI_ERR_TAG       BDMPI_ERR_TAG
#define MPI_ERR_ROOT      BDMPI_ERR_ROOT
#define MPI_ERR_RANK      BDMPI_ERR_RANK
#define MPI_ERR_TRUNCATE  BDMPI_ERR_TRUNCATE
#define MPI_ERR_OP        BDMPI_ERR_OP


/* Reduction operation types */
typedef BDMPI_Op MPI_Op;

#define MPI_OP_NULL BDMPI_OP_NULL

#define MPI_MAX       BDMPI_MAX
#define MPI_MIN       BDMPI_MIN
#define MPI_SUM       BDMPI_SUM
#define MPI_PROD      BDMPI_PROD
#define MPI_LAND      BDMPI_LAND
#define MPI_BAND      BDMPI_BAND
#define MPI_LOR       BDMPI_LOR
#define MPI_BOR       BDMPI_BOR
#define MPI_LXOR      BDMPI_LXOR
#define MPI_BXOR      BDMPI_BXOR
#define MPI_MAXLOC    BDMPI_MAXLOC
#define MPI_MINLOC    BDMPI_MINLOC



/* Data types */
typedef BDMPI_Datatype MPI_Datatype;

#define MPI_DATATYPE_NULL BDMPI_DATATYPE_NULL

#define MPI_CHAR                BDMPI_CHAR
#define MPI_SIGNED_CHAR         BDMPI_SIGNED_CHAR
#define MPI_UNSIGNED_CHAR       BDMPI_UNSIGNED_CHAR
#define MPI_BYTE                BDMPI_BYTE
#define MPI_WCHAR               BDMPI_WCHAR
#define MPI_SHORT               BDMPI_SHORT
#define MPI_UNSIGNED_SHORT      BDMPI_UNSIGNED_SHORT
#define MPI_INT                 BDMPI_INT
#define MPI_UNSIGNED            BDMPI_UNSIGNED
#define MPI_LONG                BDMPI_LONG
#define MPI_UNSIGNED_LONG       BDMPI_UNSIGNED_LONG
#define MPI_LONG_LONG_INT       BDMPI_LONG_LONG_INT
#define MPI_UNSIGNED_LONG_LONG  BDMPI_UNSIGNED_LONG_LONG
#define MPI_INT8_T              BDMPI_INT8_T
#define MPI_UINT8_T             BDMPI_UINT8_T
#define MPI_INT16_T             BDMPI_INT16_T
#define MPI_UINT16_T            BDMPI_UINT16_T
#define MPI_INT32_T             BDMPI_INT32_T
#define MPI_UINT32_T            BDMPI_UINT32_T
#define MPI_INT64_T             BDMPI_INT64_T
#define MPI_UINT64_T            BDMPI_UINT64_T
#define MPI_SIZE_T              BDMPI_SIZE_T
#define MPI_SSIZE_T             BDMPI_SSIZE_T
#define MPI_FLOAT               BDMPI_FLOAT
#define MPI_DOUBLE              BDMPI_DOUBLE
#define MPI_FLOAT_INT           BDMPI_FLOAT_INT
#define MPI_DOUBLE_INT          BDMPI_DOUBLE_INT
#define MPI_LONG_INT            BDMPI_LONG_INT
#define MPI_SHORT_INT           BDMPI_SHORT_INT
#define MPI_2INT                BDMPI_2INT


/* Communicators */
typedef BDMPI_Comm MPI_Comm;

#define MPI_COMM_NULL     BDMPI_COMM_NULL
#define MPI_COMM_WORLD    BDMPI_COMM_WORLD
#define MPI_COMM_SELF     BDMPI_COMM_SELF

/* Processes */
#define MPI_PROC_NULL BDMPI_PROC_NULL

/* Status */
typedef BDMPI_Status MPI_Status;

#define MPI_STATUS_IGNORE BDMPI_STATUS_IGNORE


/*! Request */
typedef BDMPI_Request MPI_Request;

#define MPI_REQUEST_NULL BDMPI_REQUEST_NULL


/*------------------------------------------------------------------------
* Prototypes
*-------------------------------------------------------------------------*/
#ifdef _WINDLL
#define MPI_API(type) __declspec(dllexport) type __cdecl
#elif defined(__cdecl)
#define MPI_API(type) type __cdecl
#else
#define MPI_API(type) type
#endif


#ifdef __cplusplus
extern "C" {
#endif

MPI_API(int) MPI_Print_vm_info(char *hdr);
MPI_API(int) MPI_Comm_lsize(MPI_Comm comm, int *lsize);
MPI_API(int) MPI_Comm_lrank(MPI_Comm comm, int *lrank);
MPI_API(int) MPI_Comm_nsize(MPI_Comm comm, int *nsize);
MPI_API(int) MPI_Comm_nrank(MPI_Comm comm, int *nrank);
MPI_API(int) MPI_Comm_rrank(MPI_Comm comm, int *rrank);
MPI_API(int) MPI_Reduce_init(void *sendbuf, void *recvbuf, int count,
                   MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm,
                   MPI_Request *request);
MPI_API(int) MPI_Reduce_fine(void *recvbuf, int count, MPI_Datatype datatype,
                   MPI_Op op, int root, MPI_Comm comm, MPI_Request *request);
MPI_API(int) MPI_Merge(void *sendbuf, int *sendids, int sendcount,
                   void *recvbuf, int *recvids, int *r_recvcount,
                   MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);

MPI_API(int) MPI_Init(int *argc, char **argv[]);
MPI_API(int) MPI_Finalize();
MPI_API(int) MPI_Comm_size(MPI_Comm comm, int *size);
MPI_API(int) MPI_Comm_rank(MPI_Comm comm, int *rank);
MPI_API(int) MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm);
MPI_API(int) MPI_Comm_free(MPI_Comm *comm);
MPI_API(int) MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm);
MPI_API(int) MPI_Send(void *buf, int count, MPI_Datatype datatype,
                   int dest, int tag, MPI_Comm comm);
MPI_API(int) MPI_Isend(void *buf, int count, MPI_Datatype datatype,
                   int dest, int tag, MPI_Comm comm, MPI_Request *request);
MPI_API(int) MPI_Recv(void *buf, int count, MPI_Datatype datatype,
                   int source, int tag, MPI_Comm comm, MPI_Status *status);
MPI_API(int) MPI_Irecv(void *buf, int count, MPI_Datatype datatype,
                   int source, int tag, MPI_Comm comm, MPI_Request *request);
MPI_API(int) MPI_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                   int dest, int sendtag, void *recvbuf, int recvcount,
                   MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm,
                   MPI_Status *status);
MPI_API(int) MPI_Test(MPI_Request *request, int *flag, MPI_Status *status);
MPI_API(int) MPI_Wait(MPI_Request *request, MPI_Status *status);
MPI_API(int) MPI_Waitall(int count, MPI_Request *request, MPI_Status *status);
MPI_API(int) MPI_Get_count(MPI_Status *status, MPI_Datatype datatype,
                   int *count);
MPI_API(int) MPI_Barrier(MPI_Comm comm);
MPI_API(int) MPI_Bcast(void *buf, int count, MPI_Datatype datatype,
                   int root, MPI_Comm comm);
MPI_API(int) MPI_Allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                   void *recvbuf, int recvcount, MPI_Datatype recvtype,
                   MPI_Comm comm);
MPI_API(int) MPI_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                   void *recvbuf, int *recvcounts, int *displs,
                   MPI_Datatype recvtype, MPI_Comm comm);
MPI_API(int) MPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                   void *recvbuf, int recvcount, MPI_Datatype recvtype,
                   int root, MPI_Comm comm);
MPI_API(int) MPI_Gatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                   void *recvbuf, int *recvcounts, int *displs,
                   MPI_Datatype recvtype, int root, MPI_Comm comm);
MPI_API(int) MPI_Scatter(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                   void *recvbuf, int recvcount, MPI_Datatype recvtype,
                   int root, MPI_Comm comm);
MPI_API(int) MPI_Scatterv(void *sendbuf, int *sendcounts, int *displs,
                   MPI_Datatype sendtype, void *recvbuf, int recvcount,
                   MPI_Datatype recvtype, int root, MPI_Comm comm);
MPI_API(int) MPI_Alltoall(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                   void *recvbuf, int recvcount, MPI_Datatype recvtype,
                   MPI_Comm comm);
MPI_API(int) MPI_Alltoallv(void *sendbuf, int *sendcounts, int *sdipls,
                   MPI_Datatype sendtype, void *recvbuf, int *recvcounts,
                   int *rdispls, MPI_Datatype recvtype, MPI_Comm comm);
MPI_API(int) MPI_Reduce(void *sendbuf, void *recvbuf, int count,
                   MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
MPI_API(int) MPI_Allreduce(void *sendbuf, void *recvbuf, int count,
                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
MPI_API(int) MPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status);
MPI_API(int) MPI_Iprobe(int source, int tag, MPI_Comm comm, int *flag,
                   MPI_Status *status);
MPI_API(int) MPI_Scan(void *sendbuf, void *recvbuf, int count,
                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
MPI_API(int) MPI_Exscan(void *sendbuf, void *recvbuf, int count,
                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
MPI_API(double) MPI_Wtime(void);

#ifdef __cplusplus
}
#endif

//void* _malloc(size_t, char*, int);
//#define malloc(SIZE) _malloc(SIZE, __FILE__, __LINE__)


#endif


#endif  /* _BDMPI_H_ */
