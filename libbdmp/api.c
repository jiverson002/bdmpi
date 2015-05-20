/*!
\file
\brief Public \bdmpi APIs 
\date Started 4/3/2013
\author George
*/

#include "bdmplib.h"


/* Guard for uninitialized BDMPI */
#define EXITIFNOTINIT(job) \
  do {\
    if (job == NULL) {\
      fprintf(stderr, "BDMPI_Init() has not been called yet.\n");\
      exit(EXIT_FAILURE);\
    }\
  } while(0)


/*************************************************************************/
/* Static variables */
/*************************************************************************/
static sjob_t *job=NULL;      /*!< Holds information about the slave job */



/*! 
\defgroup mpiapilist List of implemented \mpi functions 

This is the set of \mpi functions that are currently implemented in \bdmpi. 

For each of these functions, \bdmpi provides a variant that starts with the `MPI_`
prefix and a variant that starts with the `BDMPI_` prefix. The calling sequence of
the first variant is identical to the \mpi specification whereas the calling sequence
of the second variant has been modified to make it 64-bit compliant (e.g., replaced
most of the sizes that \mpi assumed that were `int` to either `size_t` or `ssize_t`).

@{
*/


/*************************************************************************/
/* BDMPI_ variant */
/*************************************************************************/
int BDMPI_Init(int *argc, char **argv[])
{
  return bdmp_Init(&job, argc, argv);
}

int BDMPI_Finalize()
{
  return bdmp_Finalize(job);
}

int BDMPI_Comm_size(BDMPI_Comm comm, int *size)
{
  return bdmp_Comm_size(job, comm, size);
}

int BDMPI_Comm_rank(BDMPI_Comm comm, int *rank)
{
  EXITIFNOTINIT(job);
  return bdmp_Comm_rank(job, comm, rank);
}

int BDMPI_Comm_dup(BDMPI_Comm comm, BDMPI_Comm *newcomm)
{
  EXITIFNOTINIT(job);
  return bdmp_Comm_dup(job, comm, newcomm);
}

int BDMPI_Comm_free(BDMPI_Comm *comm)
{
  EXITIFNOTINIT(job);
  return bdmp_Comm_free(job, comm);
}

int BDMPI_Comm_split(BDMPI_Comm comm, int color, int key, BDMPI_Comm *newcomm)
{
  EXITIFNOTINIT(job);
  return bdmp_Comm_split(job, comm, color, key, newcomm);
}

int BDMPI_Send(void *buf, size_t count, BDMPI_Datatype datatype, int dest, 
         int tag, BDMPI_Comm comm)
{
  EXITIFNOTINIT(job);
  return bdmp_Send(job, buf, count, datatype, dest, tag, comm);
}

int BDMPI_Isend(void *buf, size_t count, BDMPI_Datatype datatype, int dest, 
         int tag, BDMPI_Comm comm, BDMPI_Request *request)
{
  EXITIFNOTINIT(job);
  return bdmp_Isend(job, buf, count, datatype, dest, tag, comm, request);
}

int BDMPI_Recv(void *buf, size_t count, BDMPI_Datatype datatype, int source, 
         int tag, BDMPI_Comm comm, BDMPI_Status *status)
{
  EXITIFNOTINIT(job);
  return bdmp_Recv(job, buf, count, datatype, source, tag, comm, status);
}

int BDMPI_Irecv(void *buf, size_t count, BDMPI_Datatype datatype, int source, 
         int tag, BDMPI_Comm comm, BDMPI_Request *request)
{
  EXITIFNOTINIT(job);
  return bdmp_Irecv(job, buf, count, datatype, source, tag, comm, request);
}

int BDMPI_Probe(int source, int tag, BDMPI_Comm comm, BDMPI_Status *status)
{
  EXITIFNOTINIT(job);
  return bdmp_Probe(job, source, tag, comm, status);
}

int BDMPI_Iprobe(int source, int tag, BDMPI_Comm comm, int *flag, 
         BDMPI_Status *status)
{
  EXITIFNOTINIT(job);
  return bdmp_Iprobe(job, source, tag, comm, flag, status);
}

int BDMPI_Sendrecv(void *sendbuf, size_t sendcount, BDMPI_Datatype sendtype, 
         int dest, int sendtag, void *recvbuf, size_t recvcount, 
         BDMPI_Datatype recvtype, int source, int recvtag, BDMPI_Comm comm, 
         BDMPI_Status *status)
{
  int ierror;
  EXITIFNOTINIT(job);

  ierror = bdmp_Send(job, sendbuf, sendcount, sendtype, dest, sendtag, comm);
  if (ierror == BDMPI_SUCCESS)
    return bdmp_Recv(job, recvbuf, recvcount, recvtype, source, recvtag, comm, status);
  else
    return ierror;
}

int BDMPI_Get_count(BDMPI_Status *status, BDMPI_Datatype datatype, size_t *count)
{
  EXITIFNOTINIT(job);
  return bdmp_Get_count(job, status, datatype, count);
}

int BDMPI_Test(BDMPI_Request *request, int *flag, BDMPI_Status *status)
{
  EXITIFNOTINIT(job);
  return bdmp_Test(job, request, flag, status);
}

int BDMPI_Wait(BDMPI_Request *request, BDMPI_Status *status)
{
  EXITIFNOTINIT(job);
  return bdmp_Wait(job, request, status);
}

int BDMPI_Waitall(int count, BDMPI_Request *requests, BDMPI_Status *statuses)
{
  int i, ierror=BDMPI_SUCCESS;
  EXITIFNOTINIT(job);

  for (i=0; i<count; i++) {
    ierror = bdmp_Wait(job, &requests[i], &statuses[i]);
  }
  return ierror;
}

int BDMPI_Barrier(BDMPI_Comm comm)
{
  EXITIFNOTINIT(job);
  return bdmp_Barrier(job, comm);
}

int BDMPI_Bcast(void *buf, size_t count, BDMPI_Datatype datatype, int root, 
         BDMPI_Comm comm)
{
  EXITIFNOTINIT(job);

  return bdmp_Bcast(job, buf, count, datatype, root, comm);
}

int BDMPI_Allgather(void *sendbuf, size_t sendcount, BDMPI_Datatype sendtype,
         void *recvbuf, size_t recvcount,  BDMPI_Datatype recvtype, 
         BDMPI_Comm comm)
{
  size_t p, recvcounts[comm->size], displs[comm->size];
  EXITIFNOTINIT(job);

  for (p=0; p<comm->size; p++) {
    recvcounts[p] = recvcount;
    displs[p]     = p*recvcount;
  }

  return BDMPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, 
              displs, recvtype, comm);
}

int BDMPI_Allgatherv(void *sendbuf, size_t sendcount, BDMPI_Datatype sendtype,
         void *recvbuf, size_t *recvcounts, size_t *displs, 
         BDMPI_Datatype recvtype, BDMPI_Comm comm)
{
  EXITIFNOTINIT(job);

  return bdmp_Allgatherv(job, sendbuf, sendcount, sendtype, 
              recvbuf, recvcounts, displs, recvtype, comm);
}

int BDMPI_Reduce(void *sendbuf, void *recvbuf, size_t count, 
         BDMPI_Datatype datatype, BDMPI_Op op, int root, BDMPI_Comm comm)
{
  EXITIFNOTINIT(job);

  return bdmp_Reduce(job, sendbuf, recvbuf, count, datatype, op, root, comm);
}

int BDMPI_Allreduce(void *sendbuf, void *recvbuf, size_t count, 
         BDMPI_Datatype datatype, BDMPI_Op op, BDMPI_Comm comm)
{
  EXITIFNOTINIT(job);

  return bdmp_Allreduce(job, sendbuf, recvbuf, count, datatype, op, comm);
}

int BDMPI_Gather(void *sendbuf, size_t sendcount, BDMPI_Datatype sendtype,
         void *recvbuf, size_t recvcount, BDMPI_Datatype recvtype,
         int root, BDMPI_Comm comm)
{
  int p, npes=comm->size;
  size_t recvcounts[npes], rdispls[npes];
  EXITIFNOTINIT(job);

  for (p=0; p<npes; p++) {
    recvcounts[p] = recvcount;
    rdispls[p] = p*recvcount;
  }

  return BDMPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts,
              rdispls, recvtype, root, comm);
}

int BDMPI_Gatherv(void *sendbuf, size_t sendcount, BDMPI_Datatype sendtype,
         void *recvbuf, size_t *recvcounts, size_t *rdispls, 
         BDMPI_Datatype recvtype, int root, BDMPI_Comm comm)
{
  EXITIFNOTINIT(job);

  if (comm->nsize == 1) 
    return bdmp_Gatherv_node(job, sendbuf, sendcount, sendtype, recvbuf,
                 recvcounts, rdispls, recvtype, root, comm);
  else
    return bdmp_Gatherv_p2p(job, sendbuf, sendcount, sendtype, recvbuf,
                 recvcounts, rdispls, recvtype, root, comm);

}

int BDMPI_Scatter(void *sendbuf, size_t sendcount, BDMPI_Datatype sendtype,
         void *recvbuf, size_t recvcount, BDMPI_Datatype recvtype,
         int root, BDMPI_Comm comm)
{
  int p, npes=comm->size;
  size_t sendcounts[npes], sdispls[npes];
  EXITIFNOTINIT(job);

  for (p=0; p<npes; p++) {
    sendcounts[p] = sendcount;
    sdispls[p] = p*sendcount;
  }

  return BDMPI_Scatterv(sendbuf, sendcounts, sdispls, sendtype,
               recvbuf, recvcount, recvtype, root, comm);

}

int BDMPI_Scatterv(void *sendbuf, size_t *sendcounts, size_t *sdispls, 
         BDMPI_Datatype sendtype, void *recvbuf, size_t recvcount, 
         BDMPI_Datatype recvtype, int root, BDMPI_Comm comm)
{
  EXITIFNOTINIT(job);

  if (comm->nsize == 1) 
    return bdmp_Scatterv_node(job, sendbuf, sendcounts, sdispls, sendtype,
                 recvbuf, recvcount, recvtype, root, comm);
  else
    return bdmp_Scatterv_p2p(job, sendbuf, sendcounts, sdispls, sendtype,
                 recvbuf, recvcount, recvtype, root, comm);
}

int BDMPI_Alltoall(void *sendbuf, size_t sendcount, BDMPI_Datatype sendtype,
         void *recvbuf, size_t recvcount, BDMPI_Datatype recvtype,
         BDMPI_Comm comm)
{
  int p, npes=comm->size;
  size_t sendcounts[npes], recvcounts[npes], sdispls[npes], rdispls[npes];
  EXITIFNOTINIT(job);

  for (p=0; p<npes; p++) {
    sendcounts[p] = sendcount;
    recvcounts[p] = recvcount;
    sdispls[p] = p*sendcount;
    rdispls[p] = p*recvcount;
  }

  return BDMPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype,
               recvbuf, recvcounts, rdispls, recvtype, comm);

}

int BDMPI_Alltoallv(void *sendbuf, size_t *sendcounts, size_t *sdispls, 
         BDMPI_Datatype sendtype, void *recvbuf, size_t *recvcounts, 
         size_t *rdispls, BDMPI_Datatype recvtype, BDMPI_Comm comm)
{
  EXITIFNOTINIT(job);

  if (comm->nsize == 1) 
    return bdmp_Alltoallv_node(job, sendbuf, sendcounts, sdispls, sendtype,
                 recvbuf, recvcounts, rdispls, recvtype, comm);
  else
    return bdmp_Alltoallv_p2p(job, sendbuf, sendcounts, sdispls, sendtype,
                 recvbuf, recvcounts, rdispls, recvtype, comm);
}

int BDMPI_Scan(void *sendbuf, void *recvbuf, size_t count, 
         BDMPI_Datatype datatype, BDMPI_Op op, BDMPI_Comm comm)
{
  EXITIFNOTINIT(job);

  return bdmp_Scan(job, sendbuf, recvbuf, count, datatype, op, comm);
}

int BDMPI_Exscan(void *sendbuf, void *recvbuf, size_t count, 
         BDMPI_Datatype datatype, BDMPI_Op op, BDMPI_Comm comm)
{
  EXITIFNOTINIT(job);

  return bdmp_Exscan(job, sendbuf, recvbuf, count, datatype, op, comm);
}

double BDMPI_Wtime(void)
{
  return gk_WClockSeconds();
}



/*************************************************************************/
/* MPI_ variant */
/*************************************************************************/
int MPI_Init(int *argc, char **argv[])
{
  return BDMPI_Init(argc, argv);
}

int MPI_Finalize()
{
  return BDMPI_Finalize();
}

int MPI_Comm_size(MPI_Comm comm, int *size)
{
  return BDMPI_Comm_size(comm, size);
}

int MPI_Comm_rank(MPI_Comm comm, int *rank)
{
  return BDMPI_Comm_rank(comm, rank);
}

int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm)
{
  return BDMPI_Comm_dup(comm, newcomm);
}

int MPI_Comm_free(MPI_Comm *comm)
{
  return BDMPI_Comm_free(comm);
}

int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm)
{
  return BDMPI_Comm_split(comm, color, key, newcomm);
}

int MPI_Send(void const *buf, int count, MPI_Datatype datatype, int dest, 
         int tag, MPI_Comm comm)
{
  return BDMPI_Send((void*)buf, count, datatype, dest, tag, comm);
}

int MPI_Isend(void const *buf, int count, MPI_Datatype datatype, int dest, 
         int tag, MPI_Comm comm, MPI_Request *request)
{
  return BDMPI_Isend((void*)buf, count, datatype, dest, tag, comm, request);
}

int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, 
         int tag, MPI_Comm comm, MPI_Status *status)
{
  return BDMPI_Recv(buf, count, datatype, source, tag, comm, status);
}

int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, 
         int tag, MPI_Comm comm, MPI_Request *request)
{
  return BDMPI_Irecv(buf, count, datatype, source, tag, comm, request);
}

int MPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status)
{
  return BDMPI_Probe(source, tag, comm, status);
}

int MPI_Iprobe(int source, int tag, MPI_Comm comm, int *flag, 
         MPI_Status *status)
{
  return BDMPI_Iprobe(source, tag, comm, flag, status);
}

int MPI_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype, 
         int dest, int sendtag, void *recvbuf, int recvcount, 
         MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, 
         MPI_Status *status)
{
  return BDMPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag,
                        recvbuf, recvcount, recvtype, source, recvtag,
                        comm, status);
}

int MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count)
{
  int rcode;
  size_t rcount=0;

  rcode = BDMPI_Get_count(status, datatype, &rcount);
  *count = rcount;
  return rcode;
}

int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status)
{
  return BDMPI_Test(request, flag, status);
}

int MPI_Wait(MPI_Request *request, MPI_Status *status)
{
  return BDMPI_Wait(request, status);
}

int MPI_Waitall(int count, MPI_Request *requests, MPI_Status *statuses)
{
  return BDMPI_Waitall(count, requests, statuses);
}

int MPI_Barrier(MPI_Comm comm)
{
  return BDMPI_Barrier(comm);
}

int MPI_Bcast(void *buf, int count, MPI_Datatype datatype, int root, 
         MPI_Comm comm)
{
  return BDMPI_Bcast(buf, count, datatype, root, comm);
}

int MPI_Allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
         void *recvbuf, int recvcount,  MPI_Datatype recvtype, 
         MPI_Comm comm)
{

  return BDMPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, 
              recvtype, comm);
}

int MPI_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
         void *recvbuf, int *recvcounts,  int *displs, MPI_Datatype recvtype, 
         MPI_Comm comm)
{
  int i;
  size_t _recvcounts[comm->size], _displs[comm->size];

  for (i=0; i<comm->size; i++) {
    _recvcounts[i] = recvcounts[i];
    _displs[i]     = displs[i];
  }

  return BDMPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, _recvcounts, 
              _displs, recvtype, comm);
}

int MPI_Reduce(void const *sendbuf, void *recvbuf, int count, 
         MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
{
  return BDMPI_Reduce((void*)sendbuf, recvbuf, count, datatype, op, root, comm);
}

int MPI_Allreduce(void *sendbuf, void *recvbuf, int count, 
         MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  return BDMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
}

int MPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
         void *recvbuf, int recvcount, MPI_Datatype recvtype,
         int root, MPI_Comm comm)
{
  return BDMPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount,
              recvtype, root, comm);
}

int MPI_Gatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
         void *recvbuf, int *recvcounts, int *rdispls, 
         MPI_Datatype recvtype, int root, MPI_Comm comm)
{
  EXITIFNOTINIT(job);

  int i;
  size_t _recvcounts[comm->size], _rdispls[comm->size];

  for (i=0; i<comm->size; i++) {
    _recvcounts[i] = recvcounts[i];
    _rdispls[i]    = rdispls[i];
  }

  return BDMPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, _recvcounts, 
              _rdispls, recvtype, root, comm);

}

int MPI_Scatter(void *sendbuf, int sendcount, MPI_Datatype sendtype,
         void *recvbuf, int recvcount, MPI_Datatype recvtype,
         int root, MPI_Comm comm)
{
  return BDMPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype,
              root, comm);
}

int MPI_Scatterv(void *sendbuf, int *sendcounts, int *sdispls, 
         MPI_Datatype sendtype, void *recvbuf, int recvcount, 
         MPI_Datatype recvtype, int root, MPI_Comm comm)
{
  EXITIFNOTINIT(job);

  int i;
  size_t _sendcounts[comm->size], _sdispls[comm->size];

  for (i=0; i<comm->size; i++) {
    _sendcounts[i] = sendcounts[i];
    _sdispls[i]    = sdispls[i];
  }

  return BDMPI_Scatterv(sendbuf, _sendcounts, _sdispls, sendtype, recvbuf,
              recvcount, recvtype, root, comm);
}

int MPI_Alltoall(void *sendbuf, int sendcount, MPI_Datatype sendtype,
         void *recvbuf, int recvcount, MPI_Datatype recvtype,
         MPI_Comm comm)
{
  return BDMPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, 
              recvtype, comm);
}

int MPI_Alltoallv(void const * const sendbuf, int const * const sendcounts,
         int const * const sdispls, MPI_Datatype sendtype, void *recvbuf,
         int const * const recvcounts, int const * const rdispls,
         MPI_Datatype recvtype, MPI_Comm comm)
{
  EXITIFNOTINIT(job);

  int i;
  size_t _sendcounts[comm->size], _sdispls[comm->size],
         _recvcounts[comm->size], _rdispls[comm->size];

  for (i=0; i<comm->size; i++) {
    _sendcounts[i] = sendcounts[i];
    _sdispls[i]    = sdispls[i];
    _recvcounts[i] = recvcounts[i];
    _rdispls[i]    = rdispls[i];
  }

  return BDMPI_Alltoallv((void*)sendbuf, _sendcounts, _sdispls, sendtype, recvbuf,
              _recvcounts, _rdispls, recvtype, comm);

}

int MPI_Scan(void *sendbuf, void *recvbuf, int count, 
         MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  return BDMPI_Scan(sendbuf, recvbuf, count, datatype, op, comm);
}

int MPI_Exscan(void *sendbuf, void *recvbuf, int count, 
         MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  return BDMPI_Exscan(sendbuf, recvbuf, count, datatype, op, comm);
}

double MPI_Wtime(void)
{
  return gk_WClockSeconds();
}

/*! @} */



/*************************************************************************/
/* BDMPI specific APIs */
/*************************************************************************/

/*! 
\if html
\defgroup bdmpiapi List of \bdmpi specific functions 
@{
\endif
*/

/*! \defgroup bdmpicommlist Functions for further communicator interrogation

The description of the functions that follow assume that each compute node runs a
single instance of the `bdmprun` master process. In that case, there is a one-to-one
mapping between the set of slaves spawn by the `bdmprun` process and the set of ranks
of the \bdmpi program that execute on the node. However, when `mpiexec` assigns
multiple instances of `bdmprun` on the same compute node, then the above mapping does
not hold. In such cases, the term *node* does not correspond to a compute node but to
a `bdmprun` process and the set of ranks that execute on that node correspond to the
slave processes spawned by that `bdmprun` instance.

@{
*/


/*************************************************************************/
/*! 

It is used to get the local size of a communicator. The *local size* is the number of
ranks of the supplied communicator that are assigned to the same node as the calling
process. For example, if a communicator has 5, 3, and 8 processes assigned to three
different nodes, then the local size of all the processes in the first node will be
5, for the second node will be 3, and for the third node will be 8.

\param[in] comm is the communicator being interrogated.

\param[out] lsize returns the local size of the calling process in \p comm.
       
*/
/*************************************************************************/
int BDMPI_Comm_lsize(BDMPI_Comm comm, int *lsize)
{
  return bdmp_Comm_lsize(job, comm, lsize);
}


/*************************************************************************/
/*! 

It is used to get the local rank of the calling process within the supplied
communicator. The *local rank* of process is the number of lower-ranked processes
that are assigned to the same node in the communicator. For example, if \p comm has
20 ranks and ranks 4, 5, 8, and 15 are assigned to the same node, then the intra-node
ranks of these four ranks are 0, 1, 2, and 3, respectively.

\param[in] comm is the communicator being interrogated.

\param[out] lrank returns the intra-node rank of the calling process in the \p comm.
       
*/
/*************************************************************************/
int BDMPI_Comm_lrank(BDMPI_Comm comm, int *lrank)
{
  EXITIFNOTINIT(job);
  return bdmp_Comm_lrank(job, comm, lrank);
}


/*************************************************************************/
/*! 

It is used to get the number of nodes used by the processes of the supplied
communicator.

\param[in] comm is the communicator being interrogated.

\param[out] nsize returns the number of nodes used by the union of the 
ranks in \p comm.
       
*/
/*************************************************************************/
int BDMPI_Comm_nsize(BDMPI_Comm comm, int *nsize)
{
  return bdmp_Comm_nsize(job, comm, nsize);
}


/*************************************************************************/
/*! 

It is used to get the rank of the node on which the calling process is assigned. 

\param[in] comm is the communicator being interrogated.

\param[out] nrank returns the rank of the node in \p comm of the calling process. 
Note that all the slaves that are assigned to the same node in \p comm will return 
the same node rank.
       
*/
/*************************************************************************/
int BDMPI_Comm_nrank(BDMPI_Comm comm, int *nrank)
{
  EXITIFNOTINIT(job);
  return bdmp_Comm_nrank(job, comm, nrank);
}


/*************************************************************************/
/*! 

It is used to get the rank of the process on the same node that has a local rank of
0. This is referred to as *root rank*. For example, if a communicator has 10
processes that are distributed among three nodes as follows: {0, 1, 2, 3}, {4, 5},
and {6, 7, 8, 9}; then, the root rank for all processes in the first node will
be 0, the root rank for all processes in the second node will be 4, whereas
for the third node will be 6. 

\param[in] comm is the communicator being interrogated.

\param[out] rrank returns the root rank of calling process in \p comm.
       
*/
/*************************************************************************/
int BDMPI_Comm_rrank(BDMPI_Comm comm, int *rrank)
{
  EXITIFNOTINIT(job);
  return bdmp_Comm_rrank(job, comm, rrank);
}

/*! @} */


/*! \defgroup bdmpimutexlist Functions for critical sections 

These functions and the `-nr` option of \bdmprun are used to implement a generalized
mutex-like synchronization protocol involving the slave processes that were spawned
by the same master process. For the rest of this discussion we will refer to the set
of slaves spawned by the same master process as *related slaves*. 

In the context of \bdmpi, a *critical section* is a part of the code that can be
executed by fewer slaves than the number of related slaves that can be concurrently
running (i.e., as controlled by the `-nr` option of \bdmprun). 

\bdmpi uses POSIX semaphores to implement critical sections.

@{
*/

/**************************************************************************/
/*! 
  
BDMPI_Entercritical() is used to indicate the start of the critical section. If the
number of related slaves that are already executing in their critical sections is
equal to \c -nr, then the calling process blocks until one of the other processes
exits its critical section (by calling BDMPI_Exitcritical()). At that point, if no
other process wants to enter its critical section, the function returns control to
the calling process. If there is another process that wants to enter its critical
section, the one that called BDMPI_Entercritical() will be allowed to proceed. Note
that the blocking of a slave process performed by BDMPI_Entercritical() does not
change its execution state w.r.t. \bdmpi. That is, the slave process is still
considered to be running.

*/
/**************************************************************************/
int BDMPI_Entercritical(void)
{
  bdlock_lock(job->criticalMX);

  return BDMPI_SUCCESS;
}


/**************************************************************************/
/*! 
  
BDMPI_Entercritical() is used to indicate that the calling process is exiting its
critical section. 

*/
/**************************************************************************/
int BDMPI_Exitcritical(void)
{
  bdlock_unlock(job->criticalMX);

  return BDMPI_SUCCESS;
}


/*! @} */

/*! \if html @} \endif */



/**************************************************************************/
/*! Prints information about the resident vss */
/**************************************************************************/
int BDMPI_Print_vm_info(char *hdr)
{
  size_t vmsize, vmrss;

  gk_GetVMInfo(&vmsize, &vmrss);
  printf("[%s:%06d] VmSize: %10zu VmRSS: %10zu\n", hdr, getpid(), vmsize, vmrss);

  return BDMPI_SUCCESS;
}

int MPI_Print_vm_info(char *hdr)
{
  return BDMPI_Print_vm_info(hdr);
}


/**************************************************************************/
/*! BDMPI's custom-wrappers over the mlock/munlock functions. The wrapper
    functions ensure that the locking operations are controlled by an 
    inter-process mutex. */
/**************************************************************************/
int BDMPI_mlock(const void *addr, size_t len)
{
  int rcode;

  bdlock_lock(job->mlockMX);
  rcode = mlock(addr, len);
  bdlock_unlock(job->mlockMX);

  return rcode;
}

int BDMPI_munlock(const void *addr, size_t len)
{
  return munlock(addr, len);
}

int BDMPI_mlockall(int flags)
{
  int rcode;

  bdlock_lock(job->mlockMX);
  rcode = mlockall(flags);
  bdlock_unlock(job->mlockMX);

  return rcode;
}

int BDMPI_munlockall(void)
{
  return munlockall();
}


/*************************************************************************/
/*************************************************************************/
int BDMPI_Reduce_init(void *sendbuf, void *recvbuf, size_t count, 
         BDMPI_Datatype datatype, BDMPI_Op op, int root, BDMPI_Comm comm,
         BDMPI_Request *request)
{
  EXITIFNOTINIT(job);

  return bdmp_Reduce_init(job, sendbuf, recvbuf, count, datatype, op, root, 
               comm, request);
}


/*************************************************************************/
/*************************************************************************/
int BDMPI_Reduce_fine(void *recvbuf, size_t count, BDMPI_Datatype datatype, 
         BDMPI_Op op, int root, BDMPI_Comm comm, BDMPI_Request *request)
{
  EXITIFNOTINIT(job);

  return bdmp_Reduce_fine(job, recvbuf, count, datatype, op, root, comm, request);
}


/*************************************************************************/
/*************************************************************************/
int BDMPI_Merge(void *sendbuf, int *sendids, size_t sendcount, 
               void *recvbuf, int *recvids, size_t *r_recvcount,
               BDMPI_Datatype datatype, BDMPI_Op op, int root, BDMPI_Comm comm)
{
  EXITIFNOTINIT(job);

  return bdmp_Merge(job, sendbuf, sendids, sendcount, recvbuf, recvids, 
              r_recvcount, datatype, op, root, comm);
}

