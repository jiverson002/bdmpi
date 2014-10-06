/*!
\file
\brief Various functions for gather operations.
\date Started 5/3/2013
\author George
*/


#include "bdmprun.h"


/*************************************************************************/
/*! Response to a BDMPI_Gatherv. - Send part
    Protocol:
       Blocks the process.
       Decreases the counter associated with the counter.
       Copies the data that were sent by the slave (if any)
       If counter becomes 0, then moves all processes to runnable state
       and sets the counter back to the size of that communicator.

*/
/*************************************************************************/
void *mstr_gather_send(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  int i, srank, sleeping;
  char *buf;
  bdmcomm_t *comm;

  BD_GET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  /* hook to the key info */
  comm  = job->comms[msg->mcomm];
  srank = babel_get_srank(comm, msg->myrank);

  BDASSERT(comm->nnodes == 1);

  /* deal with copid */
  if (comm->counter == comm->lsize)
    comm->copid++;
  if (bdmq_send(job->c2sMQs[srank], &(comm->copid), sizeof(int)) == -1)
    bdprintf("Failed to send a message to %d: %s\n", srank, strerror(errno));
  msg->copid = comm->copid;


  /* if not the root, see if slave is sending you data */
  if (msg->myrank != msg->dest) {
    /* receive the specific info about this message */
    xfer_in_scb(job->scbs[srank], msg, sizeof(bdmsg_t), BDMPI_BYTE);

    if (msg->fnum == -1) {
      /* allocate memory and receive it */
      buf = gk_cmalloc(msg->count*bdmp_sizeof(msg->datatype), "gather: buf");
      xfer_in_scb(job->scbs[srank], buf, msg->count, msg->datatype);
    }
    else {
      buf = NULL;
    }

    pending_addgather(job, msg, buf, (buf == NULL ? 0 : msg->count*bdmp_sizeof(msg->datatype)));
  }

  /* block the slave, once you get the ok from it */
  xfer_in_scb(job->scbs[srank], &sleeping, sizeof(int), BDMPI_BYTE);
  slvpool_cblock(job, srank);

  /* take action of all slaves have called */
  if (--comm->counter == 0) {
    comm->counter = comm->lsize;
    for (i=0; i<comm->lsize; i++)
      slvpool_cunblock(job, comm->sranks[i]);
  }

  BD_LET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

  gk_free((void **)&arg, LTERM);

  return NULL;
}


/*************************************************************************/
/*! Response to a BDMPI_Gather - Recv step.
    Protocol:
      - Finds the headers of the individual messages (if they exist)
      - Reads and copies the data to the slave.
*/
/*************************************************************************/
void *mstr_gather_recv(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  size_t i, len, size, dsize;
  int srank;
  header_t *hdr;
  char *data;
  bdmcomm_t *comm;

  BD_GET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  /* hook to the key info */
  comm  = job->comms[msg->mcomm];
  srank = babel_get_srank(comm, msg->myrank);

  BDASSERT(comm->nnodes == 1);

  /* go and copy the data received from all the slaves to the root */
  for (i=0; i<comm->lsize-1; i++) {
    if ((hdr = pending_getgather(job, msg)) == NULL)
      slvpool_abort(1, "mstr_gather_recv: could not locate a pending_getgather.\n");

    /* send the specific info about the data been sent */
    xfer_out_scb(job->scbs[srank], &(hdr->msg), sizeof(bdmsg_t), BDMPI_BYTE);

    /* send the actual data */
    if (hdr->msg.fnum == -1)
      xfer_out_scb(job->scbs[srank], hdr->buf, hdr->msg.count, hdr->msg.datatype);

    pending_freeheader(job, &hdr);
  }

  BD_LET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

  gk_free((void **)&arg, LTERM);

  return NULL;
}
