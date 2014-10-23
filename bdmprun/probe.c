/*!
\file
\brief Various functions for probe operations.
\date Started 6/9/2013
\author George
*/


#include "bdmprun.h"


/*************************************************************************/
/*! Response to a BDMPI_Probe.
    Checks if a previous send has already been posted, if yes it proceeds
    to extract the request information.
*/
/*************************************************************************/
void *mstr_probe(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  int flag, srank;
  header_t *hdr;
  bdmcomm_t *comm;

  BD_GET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_probe: source: %d, dest: %d [entering]\n",
        job->mynode, msg->myrank, msg->source, msg->dest));

  /* hook to the key info */
  comm  = job->comms[msg->mcomm];
  srank = babel_get_srank(comm, msg->myrank);

  /* see if the send has been posted */
  pending_locksend(job, msg);
  if ((hdr = pending_getsend(job, msg, 0)) == NULL)
    slvpool_mblock(job, srank); /* the locksend guard is for the same reason as recv */
  pending_unlocksend(job, msg);

  if (hdr == NULL) { /* it has not been posted yet */
    flag = 0;
    if (bdmq_send(job->c2sMQs[srank], &flag, sizeof(int)) == -1)
      bdprintf("Failed to send a message to %d: %s\n", srank, strerror(errno));
  }
  else { /* a matching send has been posted */
    flag = 1;
    if (bdmq_send(job->c2sMQs[srank], &flag, sizeof(int)) == -1)
      bdprintf("Failed to send a message to %d: %s\n", srank, strerror(errno));

    /* send msg info to the slave */
    xfer_out_scb(job->scbs[srank], &(hdr->msg), sizeof(bdmsg_t), BDMPI_BYTE);
  }

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_probe: source: %d, dest: %d, flag: %d [exiting]\n",
        job->mynode, msg->myrank, msg->source, msg->dest, flag));

  BD_LET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

  gk_free((void **)&arg, LTERM);

  return NULL;
}


/*************************************************************************/
/*! Response to a BDMPI_Iprobe.
    Checks if a previous send has already been posted, if yes it proceeds
    to extract the request information.
*/
/*************************************************************************/
void *mstr_iprobe(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  int flag, srank;
  header_t *hdr;
  bdmcomm_t *comm;

  BD_GET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_iprobe: source: %d, dest: %d [entering]\n",
        job->mynode, msg->myrank, msg->source, msg->dest));

  /* hook to the key info */
  comm  = job->comms[msg->mcomm];
  srank = babel_get_srank(comm, msg->myrank);

  if ((hdr = pending_getsend(job, msg, 0)) == NULL) { /* it has not been posted yet */
    flag = 0;
    if (bdmq_send(job->c2sMQs[srank], &flag, sizeof(int)) == -1)
      bdprintf("Failed to send a message to %d: %s\n", srank, strerror(errno));
  }
  else { /* a matching send has been posted */
    flag = 1;
    if (bdmq_send(job->c2sMQs[srank], &flag, sizeof(int)) == -1)
      bdprintf("Failed to send a message to %d: %s\n", srank, strerror(errno));

    /* send msg info to the slave */
    xfer_out_scb(job->scbs[srank], &(hdr->msg), sizeof(bdmsg_t), BDMPI_BYTE);
  }

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_iprobe: source: %d, dest: %d, flag: %d [exiting]\n",
        job->mynode, msg->myrank, msg->source, msg->dest, flag));

  BD_LET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

  gk_free((void **)&arg, LTERM);

  return NULL;
}
