/*
 * proto.h
 *
 * This file contains function prototypes
 *
 * Started 3/31/2013
 * George
 *
 * $Id: proto.h 10513 2011-07-07 22:06:03Z karypis $
 *
 */

#ifndef _PROTO_H_
#define _PROTO_H_

/* init.c */
int bdmp_Init(sjob_t **r_job, int *argc, char **argv[]);

/* finalize.c */
int bdmp_Finalize(sjob_t *job);

/* comm.c */
int bdmp_Comm_size(sjob_t *job, BDMPI_Comm comm, int *size);
int bdmp_Comm_lsize(sjob_t *job, BDMPI_Comm comm, int *lsize);
int bdmp_Comm_nsize(sjob_t *job, BDMPI_Comm comm, int *nsize);
int bdmp_Comm_rank(sjob_t *job, BDMPI_Comm comm, int *rank);
int bdmp_Comm_lrank(sjob_t *job, BDMPI_Comm comm, int *lrank);
int bdmp_Comm_nrank(sjob_t *job, BDMPI_Comm comm, int *nrank);
int bdmp_Comm_rrank(sjob_t *job, BDMPI_Comm comm, int *rrank);
int bdmp_Comm_dup(sjob_t *job, BDMPI_Comm comm, BDMPI_Comm *newcomm);
int bdmp_Comm_free(sjob_t *job, BDMPI_Comm *comm);
int bdmp_Comm_split(sjob_t *job, BDMPI_Comm comm, int color, int key,
          BDMPI_Comm *newcomm);

/* send.c */
int bdmp_Send(sjob_t *job, void *buf, size_t count, BDMPI_Datatype datatype,
          int dest, int tag, BDMPI_Comm comm);
int bdmp_Isend(sjob_t *job, void *buf, size_t count, BDMPI_Datatype datatype,
          int dest, int tag, BDMPI_Comm comm, BDMPI_Request *request);

/* recv.c */
int bdmp_Recv(sjob_t *job, void *buf, size_t count, BDMPI_Datatype datatype,
          int source, int tag, BDMPI_Comm comm, BDMPI_Status *status);
int bdmp_Irecv(sjob_t *job, void *buf, size_t count, BDMPI_Datatype datatype,
          int source, int tag, BDMPI_Comm comm, BDMPI_Request *request);

/* probe.c */
int bdmp_Probe(sjob_t *job, int source, int tag, BDMPI_Comm comm,
          BDMPI_Status *status);
int bdmp_Iprobe(sjob_t *job, int source, int tag, BDMPI_Comm comm,
          int *flag, BDMPI_Status *status);

/* test.c */
int bdmp_Test(sjob_t *job, BDMPI_Request *request, int *flag, BDMPI_Status *status);

/* wait.c */
int bdmp_Wait(sjob_t *job, BDMPI_Request *request, BDMPI_Status *status);

/* getcount.c */
int bdmp_Get_count(sjob_t *job, BDMPI_Status *status, BDMPI_Datatype datatype,
          size_t *count);

/* barrier.c */
int bdmp_Barrier(sjob_t *job, BDMPI_Comm comm);

/* bcast.c */
int bdmp_Bcast(sjob_t *job, void *buf, size_t count, BDMPI_Datatype datatype,
          int root, BDMPI_Comm comm);

/* allgather.c */
int bdmp_Allgatherv(sjob_t *job,
          void *sendbuf, size_t sendcount, BDMPI_Datatype sendtype,
          void *recvbuf, size_t *recvcounts, size_t *displs,
          BDMPI_Datatype recvtype, BDMPI_Comm comm);

/* reduce.c */
int bdmp_Reduce(sjob_t *job, void *sendbuf, void *recvbuf, size_t count,
          BDMPI_Datatype datatype, BDMPI_Op op, int root, BDMPI_Comm comm);
int bdmp_Reduce_init(sjob_t *job, void *sendbuf, void *recvbuf, size_t count,
          BDMPI_Datatype datatype, BDMPI_Op op, int root, BDMPI_Comm comm,
          BDMPI_Request *request);
int bdmp_Reduce_fine(sjob_t *job, void *recvbuf, size_t count,
          BDMPI_Datatype datatype, BDMPI_Op op, int root, BDMPI_Comm comm,
          BDMPI_Request *request);
int bdmp_Allreduce(sjob_t *job, void *sendbuf, void *recvbuf, size_t count,
          BDMPI_Datatype datatype, BDMPI_Op op, BDMPI_Comm comm);

/* merge.c */
int bdmp_Merge(sjob_t *job, void *sendbuf, int *sendids, int sendcount,
          void *recvbuf, int *recvids, size_t *r_recvcount,
          BDMPI_Datatype datatype, BDMPI_Op op, int root, BDMPI_Comm comm);

/* alltoall.c */
int bdmp_Alltoallv_node(sjob_t *job,
          void *sendbuf, size_t *sendcounts, size_t *sdispls, BDMPI_Datatype sendtype,
          void *recvbuf, size_t *recvcounts, size_t *rdispls, BDMPI_Datatype recvtype,
          BDMPI_Comm comm);
int bdmp_Alltoallv_p2p(sjob_t *job,
          void *sendbuf, size_t *sendcounts, size_t *sdispls, BDMPI_Datatype sendtype,
          void *recvbuf, size_t *recvcounts, size_t *rdispls, BDMPI_Datatype recvtype,
          BDMPI_Comm comm);

/* scatter.c */
int bdmp_Scatterv_node(sjob_t *job,
          void *sendbuf, size_t *sendcounts, size_t *sdispls, BDMPI_Datatype sendtype,
          void *recvbuf, size_t recvcount, BDMPI_Datatype recvtype, int root,
          BDMPI_Comm comm);
int bdmp_Scatterv_p2p(sjob_t *job,
          void *sendbuf, size_t *sendcounts, size_t *sdispls, BDMPI_Datatype sendtype,
          void *recvbuf, size_t recvcount, BDMPI_Datatype recvtype, int root,
          BDMPI_Comm comm);


/* gather.c */
int bdmp_Gatherv_node(sjob_t *job,
          void *sendbuf, size_t sendcount, BDMPI_Datatype sendtype, void *recvbuf,
          size_t *recvcounts, size_t *rdispls, BDMPI_Datatype recvtype, int root,
          BDMPI_Comm comm);
int bdmp_Gatherv_p2p(sjob_t *job,
          void *sendbuf, size_t sendcount, BDMPI_Datatype sendtype, void *recvbuf,
          size_t *recvcounts, size_t *rdispls, BDMPI_Datatype recvtype, int root,
          BDMPI_Comm comm);

/* scan.c */
int bdmp_Scan(sjob_t *job, void *sendbuf, void *recvbuf, size_t count,
        BDMPI_Datatype datatype, BDMPI_Op op, BDMPI_Comm icomm);
int bdmp_Exscan(sjob_t *job, void *sendbuf, void *recvbuf, size_t count,
        BDMPI_Datatype datatype, BDMPI_Op op, BDMPI_Comm icomm);

/* debug.c */
int slv_printf(char *f_str,...);

/* sbmalloc.c */
int sb_init(char *fstem, sjob_t * const job);
int sb_finalize();
void *sb_malloc(size_t nbytes);
void *sb_realloc(void *ptr, size_t nbytes);
void sb_free(void *buf);
int sb_exists(void *ptr);
void sb_save(void *buf);
void sb_saveall();
size_t sb_saveall_internal();
void sb_load(void *buf);
void sb_loadall();
void sb_discard(void *ptr, ssize_t size);

/* route.c */
void slv_route(sjob_t * const job, bdmsg_t const * const gomsg);

#endif
