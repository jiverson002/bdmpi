/*
 * fortran_api.c
 *
 * Contains the top-level version of the BDMP API wrapped for consumption by
 * fotran. 
 *
 *
 * Started 11/02/2013
 * Dominique
 *
 */

#include "bdmplib.h"


/* size of tables for communicators and requests */
#define POOLSIZE 256

/* handle fortran's different int sizes */
#ifndef BDMPI_FORTRAN_INTEGER_4
typedef int64_t _fint_t;
#else
typedef int32_t _fint_t;
#endif


/*************************************************************************/
/* Static variables */
/*************************************************************************/
static BDMPI_Request reqs[POOLSIZE];
static _fint_t nreqs;
static BDMPI_Comm comms[POOLSIZE];
static _fint_t ncomms;


/*************************************************************************/
/* Functions */
/*************************************************************************/
void bdmpi_init_(_fint_t * err)
{
  int x = 1;
  const char * arg = "BDMPI_FORTRAN_API";
  char ** vargs = (char**)(&arg); 
  *err = (_fint_t)BDMPI_Init(&x,&vargs);
  comms[0] = BDMPI_COMM_WORLD;
  comms[1] = BDMPI_COMM_SELF;
  comms[2] = BDMPI_COMM_CWORLD;
  comms[3] = BDMPI_COMM_NODE;
  ncomms = 4;
  nreqs = 0;
}


void bdmpi_finalize_(_fint_t * err)
{
  _fint_t r;
  r = (_fint_t)BDMPI_Finalize();
  *err = r;
}


void bdmpi_comm_size_(_fint_t * comm, _fint_t *size, 
    _fint_t * err)
{
  int isize;
  _fint_t r;
  r = (_fint_t)BDMPI_Comm_size(comms[*comm],&isize);
  *size = (_fint_t)isize;
  *err = r;
}


void bdmpi_comm_rank_(_fint_t * comm, _fint_t *rank, 
    _fint_t * err)
{
  int irank;
  _fint_t r;
  r = (_fint_t)BDMPI_Comm_rank(comms[*comm],&irank);
  *rank = (_fint_t)irank;
  *err = r;
}


void bdmpi_comm_lsize_(_fint_t *comm, _fint_t *lsize, 
    _fint_t * err)
{
  int ilsize;
  _fint_t r;
  r = (_fint_t)BDMPI_Comm_lsize(comms[*comm], &ilsize);
  *lsize = (_fint_t)ilsize;
  *err = r;
  
}


void bdmpi_comm_lrank_(_fint_t *comm, _fint_t *lrank, 
    _fint_t * err)
{
  int ilrank;
  _fint_t r;
  r = (_fint_t)BDMPI_Comm_lrank(comms[*comm], &ilrank);
  *lrank=(_fint_t)ilrank;
  *err=r;
}


void bdmpi_comm_nsize_(_fint_t *comm, _fint_t *nsize, 
    _fint_t * err)
{
  int insize;
  _fint_t r;
  r = (_fint_t)BDMPI_Comm_nsize(comms[*comm], &insize);
  *nsize = (_fint_t)insize;
  *err=r;
}


void bdmpi_comm_nrank_(_fint_t *comm, _fint_t *nrank, 
    _fint_t * err)
{
  int inrank;
  _fint_t r;
  r = (_fint_t)BDMPI_Comm_nrank(comms[*comm], &inrank);
  *nrank = (_fint_t)inrank;
  *err=r;
}

void bdmpi_comm_rrank_(_fint_t *comm, _fint_t *rrank, 
    _fint_t * err)
{
  int irrank;
  _fint_t r;
  r = (_fint_t)BDMPI_Comm_rrank(comms[*comm], &irrank);
  *rrank = (_fint_t)irrank;
  *err=r;
}


void bdmpi_comm_dup_(_fint_t *comm, _fint_t *newcomm, 
    _fint_t * err)
{
  _fint_t r;
  r = (_fint_t)BDMPI_Comm_dup(comms[*comm], comms+ncomms);
  *newcomm = ncomms++;
  *err=r;
}


void bdmpi_comm_free_(_fint_t ** comm, _fint_t * err)
{
  _fint_t r;
  r = (_fint_t)BDMPI_Comm_free(comms+(**comm));
  *err=r;
}


void bdmpi_comm_split_(_fint_t *comm, _fint_t *color, 
    _fint_t *key, _fint_t *newcomm, _fint_t * err)
{
  _fint_t r;
  r = (_fint_t)BDMPI_Comm_split(comms[*comm], *color, *key, comms+ncomms);
  *newcomm=ncomms++;
  *err=r;
}


void bdmpi_send_(void *buf, _fint_t *count, _fint_t*datatype, 
    _fint_t *dest, _fint_t *tag, _fint_t *comm, 
    _fint_t * err)
{
  _fint_t r;
  r = (_fint_t)BDMPI_Send(buf, *count,*datatype, *dest, *tag, comms[*comm]);
  *err=r;
}


void bdmpi_isend_(void *buf, _fint_t *count, _fint_t*datatype, 
    _fint_t *dest, _fint_t *tag, _fint_t *comm, 
    _fint_t *request, _fint_t * err)
{
  _fint_t r;
  if (nreqs >= POOLSIZE) {
    nreqs = 0;
  }
  r = (_fint_t)BDMPI_Isend(buf, *count,*datatype, *dest, *tag, comms[*comm], 
      reqs+nreqs);
  if (r != BDMPI_SUCCESS) {
    *request = -1;
  } else {
    *request = nreqs++;
  }
  *err=r;
}


void bdmpi_recv_(void *buf, _fint_t *count, _fint_t*datatype, 
    _fint_t *source, _fint_t *tag, _fint_t *comm, 
    _fint_t *status, _fint_t * err)
{
  _fint_t r;
  r = (_fint_t)BDMPI_Recv(buf, *count,*datatype, *source, *tag, comms[*comm],
      (BDMPI_Status*)status);
  *err=r;
}


void bdmpi_irecv_(void *buf, _fint_t *count, _fint_t*datatype, 
    _fint_t *source, _fint_t *tag, _fint_t *comm, 
    _fint_t *request, _fint_t * err)
{
  _fint_t r;
  if (nreqs >= POOLSIZE) {
    nreqs = 0;
  }
  r = (_fint_t)BDMPI_Irecv(buf, *count,*datatype, *source, *tag, comms[*comm],
      reqs+nreqs);
  if (r != BDMPI_SUCCESS) {
    *request = -1;
  } else {
    *request = nreqs++;
  }
  *err=r;
}


void bdmpi_sendrecv_(void *sendbuf, _fint_t *sendcount, 
    _fint_t * sendtype, _fint_t * dest, _fint_t * sendtag, 
    void *recvbuf, _fint_t *recvcount, _fint_t * recvtype, 
    _fint_t * source, _fint_t * recvtag, _fint_t *comm, 
    _fint_t *status, _fint_t * err)
{
  _fint_t r;
  r = (_fint_t)BDMPI_Sendrecv(sendbuf, *sendcount, *sendtype, *dest, *sendtag, 
      recvbuf, *recvcount, *recvtype, *source, *recvtag, comms[*comm], 
      (BDMPI_Status*)status);
  *err = r;
}


void bdmpi_test_(_fint_t *request, _fint_t *flag, _fint_t *status, _fint_t * err)
{
  int iflag;
  _fint_t r;
  r = (_fint_t)BDMPI_Test(reqs+(*request), &iflag, (BDMPI_Status*)status);
  *flag = (_fint_t)iflag;
  *err = r;
}


void bdmpi_wait_(_fint_t *request, _fint_t *status, _fint_t * err)
{
  _fint_t r;
  if (*request != -1) {
    r = (_fint_t)BDMPI_Wait(reqs+(*request), (BDMPI_Status*)status);
  } else {
    r = BDMPI_SUCCESS;
  }
  *err = r;
}


void bdmpi_waitall_(_fint_t *count, _fint_t *request, _fint_t *status, _fint_t * err)
{
  _fint_t r,i,n;
  n = *count;
  BDMPI_Request * myreqs;
  myreqs = (BDMPI_Request*)malloc(sizeof(BDMPI_Request)*n);
  for (i=0;i<n;++i) {
    if (request[i] == -1) {
      myreqs[i] = BDMPI_REQUEST_NULL;
    } else {
      myreqs[i] = reqs[request[i]];
    }
  }
  r = (_fint_t)BDMPI_Waitall(*count, myreqs, (BDMPI_Status*)status);
  free(myreqs);
  *err = r;
}


void bdmpi_get_count_(_fint_t *status, _fint_t*datatype, 
    _fint_t *count, _fint_t * err)
{
  _fint_t r;
  size_t rcount;
  r = (_fint_t)BDMPI_Get_count((BDMPI_Status*)status,*datatype, &rcount);
  *count = (_fint_t)rcount;
  *err = r;
}


void bdmpi_barrier_(_fint_t *comm, _fint_t * err)
{
  _fint_t r;
  r = (_fint_t)BDMPI_Barrier(comms[*comm]);
  *err = r;
}


void bdmpi_bcast_(void *buf, _fint_t *count, _fint_t*datatype, _fint_t *root, 
    _fint_t *comm, _fint_t * err)
{
  _fint_t r;
  r = (_fint_t)BDMPI_Bcast(buf, *count,*datatype, *root, comms[*comm]);
  *err = r;
}


void bdmpi_allgather_(void *sendbuf, _fint_t *sendcount, _fint_t *sendtype, 
    void *recvbuf, _fint_t *recvcount, _fint_t *recvtype, _fint_t *comm,
    _fint_t * err)
{
  _fint_t r;
  r = (_fint_t)BDMPI_Allgather(sendbuf, *sendcount, *sendtype, recvbuf, *recvcount, 
      *recvtype, comms[*comm]);
  *err = r;
}


void bdmpi_allgatherv_(void *sendbuf, _fint_t *sendcount, 
    _fint_t *sendtype, void *recvbuf, _fint_t *recvcounts, _fint_t *displs, 
    _fint_t *recvtype, _fint_t *comm, _fint_t * err)
{
  _fint_t r,i;
  int size;
  BDMPI_Comm_size(comms[*comm],&size);
  size_t * rrecvcounts = (size_t*)malloc(sizeof(size_t)*size);
  size_t * rdispls = (size_t*)malloc(sizeof(size_t)*size);
  r = (_fint_t)BDMPI_Allgatherv(sendbuf, *sendcount, *sendtype, recvbuf, rrecvcounts, 
      rdispls, *recvtype, comms[*comm]);
  for (i=0;i<size;++i) {
    recvcounts[i] = (_fint_t)rrecvcounts[i];
  }
  free(rrecvcounts);
  for (i=0;i<size;++i) {
    displs[i] = (_fint_t)rdispls[i];
  }
  free(rdispls);
  *err = r;
}


void bdmpi_gather_(void *sendbuf, _fint_t *sendcount, _fint_t *sendtype, 
    void *recvbuf, _fint_t *recvcount, _fint_t *recvtype, _fint_t *root, 
    _fint_t *comm, _fint_t * err)
{
  _fint_t r;
  r = (_fint_t)BDMPI_Gather(sendbuf, *sendcount, *sendtype, recvbuf, *recvcount, 
      *recvtype, *root, comms[*comm]);
  *err = r;
}


void bdmpi_gatherv_(void *sendbuf, _fint_t *sendcount, _fint_t *sendtype, 
    void *recvbuf, _fint_t *recvcounts, _fint_t *displs, _fint_t *recvtype, 
    _fint_t *root, _fint_t *comm, _fint_t * err)
{
  _fint_t r,i;
  int size, rank;
  size_t * rrecvcounts = NULL;
  size_t * rdispls = NULL;
  BDMPI_Comm_size(comms[*comm],&size);
  BDMPI_Comm_rank(comms[*comm],&rank);
  
  if (rank == *root) {
    rrecvcounts = (size_t*)malloc(sizeof(size_t)*size);
    rdispls = (size_t*)malloc(sizeof(size_t)*size);
    r = (_fint_t)BDMPI_Gatherv(sendbuf, *sendcount, *sendtype, recvbuf, rrecvcounts, 
        rdispls, *recvtype, *root, comms[*comm]);
    for (i=0;i<size;++i) {
      recvcounts[i] = (_fint_t)rrecvcounts[i];
    }
    free(rrecvcounts);
    for (i=0;i<size;++i) {
      displs[i] = (_fint_t)rdispls[i];
    }
    free(rdispls);
  } else {
    r = (_fint_t)BDMPI_Gatherv(sendbuf, *sendcount, *sendtype, recvbuf, rrecvcounts, 
        rdispls, *recvtype, *root, comms[*comm]);
  }
  *err = r;
}


void bdmpi_scatter_(void *sendbuf, _fint_t *sendcount, _fint_t *sendtype, 
    void *recvbuf, _fint_t *recvcount, _fint_t *recvtype, _fint_t *root, 
    _fint_t *comm, _fint_t * err)
{
  _fint_t r;
  r = (_fint_t)BDMPI_Scatter(sendbuf, *sendcount, *sendtype, recvbuf, *recvcount, 
      *recvtype, *root, comms[*comm]);
  *err = r;
}


/*
void bdmpi_scatterv_(void *sendbuf, _fint_t *sendcounts, _fint_t *displs, 
    _fint_t sendtype, void *recvbuf, _fint_t recvcount, 
    _fint_t recvtype, _fint_t root, _fint_t *comm, _fint_t * err)
{
  return BDMPI_Scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, 
      recvcount, recvtype, root, comms[*comm]);
}

*/

void bdmpi_alltoall_(void *sendbuf, _fint_t *sendcount, _fint_t *sendtype, 
    void *recvbuf, _fint_t *recvcount, _fint_t *recvtype, _fint_t *comm, 
    _fint_t * err)
{
  _fint_t r;
  r = (_fint_t)BDMPI_Alltoall(sendbuf, *sendcount, *sendtype, recvbuf, *recvcount, 
      *recvtype, comms[*comm]);
  *err = r;
}


/*
void bdmpi_alltoallv_(void *sendbuf, _fint_t *sendcounts, _fint_t *sdipls, 
    _fint_t sendtype, void *recvbuf, _fint_t *recvcounts, 
    _fint_t *rdispls, _fint_t recvtype, _fint_t *comm, _fint_t * err)
{
  return BDMPI_Alltoallv(sendbuf, sendcounts, sdipls, sendtype, recvbuf, 
      recvcounts, rdispls, recvtype, comms[*comm]);
}

*/

void bdmpi_reduce_(void *sendbuf, void *recvbuf, _fint_t *count, 
    _fint_t*datatype, _fint_t *op, _fint_t *root, _fint_t *comm, _fint_t * err)
{
  _fint_t r;
  r = (_fint_t)BDMPI_Reduce(sendbuf, recvbuf, *count,*datatype, *op, *root, 
      comms[*comm]);
  *err = r;
}


void bdmpi_reduce_init_(void *sendbuf, void *recvbuf, _fint_t *count, 
    _fint_t*datatype, _fint_t *op, _fint_t *root, _fint_t *comm, 
    BDMPI_Request *request, _fint_t * err)
{
  _fint_t r;
  r = (_fint_t)BDMPI_Reduce_init(sendbuf, recvbuf, *count,*datatype, *op, *root, 
      comms[*comm], request);
  *err = r;
}


void bdmpi_reduce_fine_(void *recvbuf, _fint_t *count, _fint_t*datatype, 
    _fint_t *op, _fint_t *root, _fint_t *comm, BDMPI_Request *request, 
    _fint_t * err)
{
  _fint_t r;
  r = (_fint_t)BDMPI_Reduce_fine(recvbuf, *count,*datatype, *op, *root, comms[*comm], 
      request);
  *err = r;
}


void bdmpi_allreduce_(void *sendbuf, void *recvbuf, _fint_t *count, 
    _fint_t*datatype, _fint_t *op, _fint_t *comm, _fint_t * err)
{
  _fint_t r;
  r = (_fint_t)BDMPI_Allreduce(sendbuf, recvbuf, *count,*datatype, *op, comms[*comm]);
  *err = r;
}


/*
void bdmpi_merge_(void *sendbuf, _fint_t *sendids, _fint_t sendcount, void *recvbuf, 
    _fint_t *recvids, _fint_t *r_recvcount, _fint_t*datatype, _fint_t * op, 
    _fint_t root, _fint_t *comm, _fint_t * err)
{
  return BDMPI_Merge(sendbuf, sendids, sendcount, recvbuf, recvids, 
      r_recvcount,*datatype, op, root, comms[*comm]);
}

*/

void bdmpi_probe_(_fint_t *source, _fint_t *tag, _fint_t *comm, 
    _fint_t *status, _fint_t * err)
{
  _fint_t r;
  r = (_fint_t)BDMPI_Probe(*source, *tag, comms[*comm], (BDMPI_Status*)status);
  *err=r;
}


void bdmpi_iprobe_(_fint_t *source, _fint_t *tag, _fint_t *comm, _fint_t *flag, 
    _fint_t *status, _fint_t * err)
{
  _fint_t r;
  int iflag;
  r = (_fint_t)BDMPI_Iprobe(*source, *tag, comms[*comm], &iflag, 
      (BDMPI_Status*)status);
  *flag = (_fint_t)iflag;
  *err=r;
}


void bdmpi_scan_(void *sendbuf, void *recvbuf, _fint_t *count, 
    _fint_t *datatype, _fint_t * op, _fint_t *comm, _fint_t * err)
{
  _fint_t r;
  r = (_fint_t)BDMPI_Scan(sendbuf, recvbuf, *count,*datatype, *op, comms[*comm]);
  *err=r;
}


void bdmpi_exscan_(void *sendbuf, void *recvbuf, _fint_t *count, 
    _fint_t*datatype, _fint_t * op, _fint_t *comm, _fint_t * err)
{
  _fint_t r;
  r = (_fint_t)BDMPI_Exscan(sendbuf, recvbuf, *count,*datatype, *op, comms[*comm]);
  *err=r;
}


/*************************************************************************/
/* Return functions declared in bdmpif.h */
/*************************************************************************/

double bdmpi_wtime_(void)
{
  return BDMPI_Wtime();
}




