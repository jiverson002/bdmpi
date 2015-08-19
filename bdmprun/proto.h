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


/* cmdline.c */
mjob_t *parse_cmdline(int argc, char *argv[]);

/* bdmprun.c */
void spawn_slaves(mjob_t *job);
void setup_master_prefork(mjob_t *job);
void setup_master_postfork(mjob_t *job);
void cleanup_master(mjob_t *job);

/* comm.c */
void comm_setup(mjob_t *job);
void comm_cleanup(mjob_t *job);
void comm_free(mjob_t *job, int comm);
void comm_realloc_master(mjob_t *job);
void *mstr_comm_dup(void *arg);
void *mstr_comm_free(void *arg);
void *mstr_comm_split(void *arg);
void comm_sort_splitdata(size_t n, splitdata_t *skeys);


/* slvpool.c */
void slvpool_listen(mjob_t *job);
void slvpool_check_deadlock(mjob_t *job);
void slvpool_showstats(mjob_t *job);
void slvpool_route(mjob_t *job, bdmsg_t *msg);
void slvpool_service_xcomms(mjob_t *job);
void slvpool_wakeup_some(mjob_t *job);
int slvpool_select_task_to_wakeup(mjob_t *job, int type);
void slvpool_remove(mjob_t *job, int rank);
void slvpool_remove_pid(mjob_t *job, pid_t cpid);
void slvpool_abort(int all, char *f_str,...);
void slvpool_kill_all(mjob_t *job);
void slvpool_munblock(mjob_t *job, int rank);
void slvpool_cunblock(mjob_t *job, int rank);
void slvpool_mblock(mjob_t *job, int rank);
void slvpool_cblock(mjob_t *job, int rank);

/* init.c */
void mstr_init(mjob_t *job, bdmsg_t *msg);

/* finalize.c */
void mstr_finalize(mjob_t *job, bdmsg_t *msg);

/* send.c */
void *mstr_send(void *arg);
void *mstr_send_remote(void *arg);

/* recv.c */
void *mstr_recv(void *arg);
void *mstr_irecv(void *arg);
void *mstr_recv_remote(void *arg);
void *mstr_recvd_remote(void *arg);

/* probe.c */
void *mstr_probe(void *arg);
void *mstr_iprobe(void *arg);

/* barrier.c */
void *mstr_barrier(void *arg);

/* bcast.c */
void *mstr_bcast_init(void *arg);
void *mstr_bcast_recv(void *arg);

/* allgather.c */
void *mstr_allgather_send(void *arg);
void *mstr_allgather_recv(void *arg);

/* reduce.c */
void *mstr_reduce_send(void *arg);
void *mstr_reduce_recv(void *arg);
void *mstr_allreduce_send(void *arg);
void *mstr_allreduce_recv(void *arg);

/* alltoall.c */
void *mstr_alltoall_send(void *arg);
void *mstr_alltoall_recv(void *arg);

/* scatter.c */
void *mstr_scatter_send(void *arg);
void *mstr_scatter_recv(void *arg);

/* gather.c */
void *mstr_gather_send(void *arg);
void *mstr_gather_recv(void *arg);


/* pending.c */
void pending_setup(mjob_t *job);
void pending_cleanup(mjob_t *job);
void pending_freeheader(mjob_t *job, header_t **r_hdr);
void      pending_addsend(mjob_t *job, bdmsg_t *msg, void *buf, size_t len);
header_t *pending_getsend(mjob_t *job, bdmsg_t *msg, int rmheader);
void      pending_locksend(mjob_t *job, bdmsg_t *msg);
void      pending_unlocksend(mjob_t *job, bdmsg_t *msg);
void      pending_addbcast(mjob_t *job, bdmsg_t *msg, void *buf, size_t len, int icnt);
header_t *pending_getbcast(mjob_t *job, bdmsg_t *msg, int countdown, int *r_counter);
void      pending_addreduce(mjob_t *job, bdmsg_t *msg, void *buf, size_t len, int icnt);
header_t *pending_getreduce(mjob_t *job, bdmsg_t *msg);
void      pending_addallreduce(mjob_t *job, bdmsg_t *msg, void *buf, size_t len);
header_t *pending_getallreduce(mjob_t *job, bdmsg_t *msg, int *r_counter);
void      pending_addmerge(mjob_t *job, bdmsg_t *msg, void *buf, int icnt);
header_t *pending_getmerge(mjob_t *job, bdmsg_t *msg);
void      pending_addallgather(mjob_t *job, bdmsg_t *msg, void *buf, size_t len);
header_t *pending_getallgather(mjob_t *job, bdmsg_t *msg, int *r_counter);
void      pending_addalltoall(mjob_t *job, bdmsg_t *msg, void *buf, size_t len);
header_t *pending_getalltoall(mjob_t *job, bdmsg_t *msg);
void      pending_addscatter(mjob_t *job, bdmsg_t *msg, void *buf, size_t len);
header_t *pending_getscatter(mjob_t *job, bdmsg_t *msg);
void      pending_addgather(mjob_t *job, bdmsg_t *msg, void *buf, size_t len);
header_t *pending_getgather(mjob_t *job, bdmsg_t *msg);

/* babel.c */
int babel_get_srank(bdmcomm_t *comm, int grank);
int babel_get_lrank(bdmcomm_t *comm, int grank);
int babel_get_grank(bdmcomm_t *comm, int lrank);
int babel_get_node(bdmcomm_t *comm, int grank);
int babel_is_local(bdmcomm_t *comm, int grank);
int babel_get_my_mcomm(mjob_t *job, int mpi_commid);

/* mpi.c */
MPI_Datatype mpi_dt(BDMPI_Datatype datatype);
MPI_Op mpi_op(BDMPI_Op op);

/* mnode.c */
void mnode_service_xcomms(mjob_t *job);
int mnode_get_next_commid(mjob_t *job);

/* merge.c */
void *mstr_merge_send(void *arg);
void *mstr_merge_recv(void *arg);
void merge_infos(mergeinfo_t *a, mergeinfo_t *b, BDMPI_Datatype dtype, BDMPI_Op op);

/* memory.c */
void * mstr_mem_load(void * const arg);
void * mstr_mem_save(void * const arg);
void memory_wakeup_some(mjob_t * const job, int const source, size_t const size);
int memory_select_task_to_wakeup(mjob_t *job, int type);

#endif
