/*!
\file
\brief Various functions for p2p recv operations.
\date Started 4/6/2013
\author George
*/


#include "bdmprun.h"


/*************************************************************************/
/*! Response to a BDMPI_Recv.
    Checks if a previous send has already been posted, if yes it proceeds
    to service the request, otherwise it denies it, which ends up blocking
    the process.
*/
/*************************************************************************/
void *mstr_recv(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  int response, srank, commid, source_node;
  bdmsg_t mmsg;
  header_t *hdr;
  bdmcomm_t *comm;

  BD_GET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_recv: source: %d, dest: %d [entering]\n",
        job->mynode, msg->myrank, msg->source, msg->dest));

  /* hook to the key info */
  comm   = job->comms[msg->mcomm];
  commid = comm->mpi_commid;
  srank  = babel_get_srank(comm, msg->myrank);

  /* see if the send has been posted */
  pending_locksend(job, msg);
  if ((hdr = pending_getsend(job, msg, 1)) == NULL) {
    M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_recv: mblocking srank: %d\n",
          job->mynode, msg->myrank, srank));

    /* this is within a pending_locksend, to prevent the race condition
       in which the slave misses a wakeup due to the arrival of a message
       before its set to blocked */
    slvpool_mblock(job, srank);
  }
  pending_unlocksend(job, msg);

  if (hdr == NULL) {
    /* it has not been posted yet */
    response = 0;
    if (bdmq_send(job->c2sMQs[srank], &response, sizeof(int)) == -1)
      bdprintf("Failed to send a message to %d: %s\n", srank, strerror(errno));

#if 1
    /* Notify remote master that a receive request has been issued. */
    if (!babel_is_local(comm, msg->source)) {
      /* Only send the receive notification the first time that
       * BDMPI_TYPE_RECV is received. */
      if (1 == msg->new_request) {
        source_node = babel_get_node(comm, msg->source);

        mmsg.mcomm   = commid;
        mmsg.source  = msg->source;
        mmsg.dest    = msg->myrank;
        mmsg.msgtype = BDMPI_MSGTYPE_RECV;

        /* Send the message header using the global node number of wcomm */
        BDASSERT(MPI_SUCCESS == MPI_Send(&mmsg, sizeof(bdmsg_t), MPI_BYTE,\
          comm->wnranks[source_node], BDMPI_HDR_TAG, job->mpi_wcomm));
      }
    }
#endif
  }
  else {
    /* a matching send has been posted */
    response = 1;
    if (bdmq_send(job->c2sMQs[srank], &response, sizeof(int)) == -1)
      bdprintf("Failed to send a message to %d: %s\n", srank, strerror(errno));

    /* send msg info to the slave */
    xfer_out_scb(job->scbs[srank], &(hdr->msg), sizeof(bdmsg_t), BDMPI_BYTE);

    if (hdr->msg.fnum == -1) /* in memory message */
      xfer_out_scb(job->scbs[srank], hdr->buf, hdr->msg.count, hdr->msg.datatype);

    pending_freeheader(job, &hdr);

#if 1
    /* Notify remote master that a receive request has been completed. */
    if (!babel_is_local(comm, msg->source)) {
      /* Only send the receive notification if this is not the first time that
       * BDMPI_TYPE_RECV is received, meaning that a message of type
       * BDMPI_MSGTYPE_RECV was sent for this message. */
      if (0 == msg->new_request) {
        source_node = babel_get_node(comm, msg->source);

        mmsg.mcomm   = commid;
        mmsg.source  = msg->source;
        mmsg.dest    = msg->myrank;
        mmsg.msgtype = BDMPI_MSGTYPE_RECVD;

        /* Send the message header using the global node number of wcomm */
        BDASSERT(MPI_SUCCESS == MPI_Send(&mmsg, sizeof(bdmsg_t), MPI_BYTE,\
          comm->wnranks[source_node], BDMPI_HDR_TAG, job->mpi_wcomm));
      }
    }
#endif
  }

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_recv: source: %d, dest: %d [exiting]\n",
        job->mynode, msg->myrank, msg->source, msg->dest));

  BD_LET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

  gk_free((void **)&arg, LTERM);

  return NULL;
}


/*************************************************************************/
/*! Response to a BDMPI_Irecv.
    Checks if a previous send has already been posted, if yes it proceeds
    to service the request, otherwise it denies it.
*/
/*************************************************************************/
void *mstr_irecv(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  int response, srank, commid, source_node;
  bdmsg_t mmsg;
  header_t *hdr;
  bdmcomm_t *comm;

  BD_GET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_irecv: source: %d, dest: %d [entering]\n",
        job->mynode, msg->myrank, msg->source, msg->dest));

  /* hook to the key info */
  comm   = job->comms[msg->mcomm];
  commid = comm->mpi_commid;
  srank  = babel_get_srank(comm, msg->myrank);

  /* see if the send has been posted */
  if ((hdr = pending_getsend(job, msg, 1)) == NULL) { /* it has not been posted yet */
    M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_recv: mblocking srank: %d\n",
          job->mynode, msg->myrank, srank));

    response = 0;
    if (bdmq_send(job->c2sMQs[srank], &response, sizeof(int)) == -1)
      bdprintf("Failed to send a message to %d: %s\n", srank, strerror(errno));
  }
  else { /* a matching send has been posted */
    response = 1;
    if (bdmq_send(job->c2sMQs[srank], &response, sizeof(int)) == -1)
      bdprintf("Failed to send a message to %d: %s\n", srank, strerror(errno));

    /* send msg info to the slave */
    xfer_out_scb(job->scbs[srank], &(hdr->msg), sizeof(bdmsg_t), BDMPI_BYTE);

    if (hdr->msg.fnum == -1) /* in memory message */
      xfer_out_scb(job->scbs[srank], hdr->buf, hdr->msg.count, hdr->msg.datatype);

    pending_freeheader(job, &hdr);
  }

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_irecv: source: %d, dest: %d [exiting]\n",
        job->mynode, msg->myrank, msg->source, msg->dest));

  BD_LET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

  gk_free((void **)&arg, LTERM);

  return NULL;
}


/*************************************************************************/
/*! Response to a BDMPI_MSGTYPE_RECV.
    Increments the count of pending recvs for the appropriate slave.
*/
/*************************************************************************/
void *mstr_recv_remote(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  int slv, mcomm;
  bdmcomm_t *comm;

  mcomm = babel_get_my_mcomm(job, msg->mcomm);
  comm  = job->comms[mcomm];
  slv   = babel_get_srank(comm, msg->source);

  assert(!babel_is_local(comm, msg->dest));

  BD_GET_LOCK(job->schedule_lock);
  job->npending[slv]++;
  BD_LET_LOCK(job->schedule_lock);

  gk_free((void **)&arg, LTERM);

  return NULL;
}


/*************************************************************************/
/*! Response to a BDMPI_MSGTYPE_RECVD.
    Decrements the count of pending recvs for the appropriate slave.
*/
/*************************************************************************/
void *mstr_recvd_remote(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  int slv, mcomm;
  bdmcomm_t *comm;

  mcomm = babel_get_my_mcomm(job, msg->mcomm);
  comm  = job->comms[mcomm];
  slv   = babel_get_srank(comm, msg->source);

  assert(!babel_is_local(comm, msg->dest));

  BD_GET_LOCK(job->schedule_lock);
  BDASSERT(job->npending[slv] > 0);
  job->npending[slv]--;
  BD_LET_LOCK(job->schedule_lock);

  gk_free((void **)&arg, LTERM);

  return NULL;
}
