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
  int response, srank;
  header_t *hdr;
  bdmcomm_t *comm;

  BD_GET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_recv: source: %d, dest: %d [entering]\n", 
        job->mynode, msg->myrank, msg->source, msg->dest));

  /* hook to the key info */
  comm  = job->comms[msg->mcomm];
  srank = babel_get_srank(comm, msg->myrank);

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
  int response, srank;
  header_t *hdr;
  bdmcomm_t *comm;

  BD_GET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_irecv: source: %d, dest: %d [entering]\n", 
        job->mynode, msg->myrank, msg->source, msg->dest));

  /* hook to the key info */
  comm  = job->comms[msg->mcomm];
  srank = babel_get_srank(comm, msg->myrank);

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

