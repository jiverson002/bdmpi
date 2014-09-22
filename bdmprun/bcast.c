/*!
\file
\brief Various functions for broadcast operations.
\date Started 4/6/2013
\author George
*/


#include "bdmprun.h"



/*************************************************************************/
/*! Response to a BDMPI_Bcast. - Init part
    Protocol:
       Blocks the process.
       Decreases the counter associated with the counter.
       Copies the data, if myrank==root
       If counter becomes 0, then moves all processes to runnable state
       and sets the counter back to the size of that communicator.

    Meaning of fields of msg:
       msg->source is the root of the broadcast (ie., all slaves in the
                 collective call should have the same root)
       msg->dest has the destination of the broadcast (ie., the slave that
                 sent a message to the master)

*/
/*************************************************************************/
void *mstr_bcast_init(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  int i, srank, sleeping;
  bdmcomm_t *comm;
  char *buf=NULL;
  size_t msize;
  header_t *hdr;

  BD_GET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  /* hook to the key info */
  comm  = job->comms[msg->mcomm];
  srank = babel_get_srank(comm, msg->myrank);

  /* message size in bytes */
  msize = bdmp_msize(msg->count, msg->datatype);

  /* deal with copid */
  if (comm->counter == comm->lsize) 
    comm->copid++;
  if (bdmq_send(job->c2sMQs[srank], &(comm->copid), sizeof(int)) == -1)
    bdprintf("Failed to send a message to %d: %s\n", srank, strerror(errno));
  msg->copid = comm->copid;


  /* see if the root send the broadcast */
  if (msg->myrank == msg->source) {
    /* allocate memory and get the data */
    buf = gk_cmalloc(msize, "bcast: buf");
    xfer_in_scb(job->scbs[srank], buf, msg->count, msg->datatype); 
    pending_addbcast(job, msg, buf, msize, comm->lsize);
  }

  /* block the slave, once you get the ok from it */
  xfer_in_scb(job->scbs[srank], &sleeping, sizeof(int), BDMPI_BYTE);
  slvpool_cblock(job, srank);

  /* take action if all slaves have called */
  if (--comm->counter == 0) {
    if (babel_is_local(comm, msg->source)) 
      hdr = pending_getbcast(job, msg, 1, NULL);

    if (comm->nnodes > 1) {  /* do something if more than one node is involved */
      BD_LET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

      if (babel_is_local(comm, msg->source)) {
        BDASSERT(MPI_Bcast(hdr->buf, msg->count, mpi_dt(msg->datatype), 
                     comm->mynode, comm->mpi_comm)
            == MPI_SUCCESS);

        if (hdr->counter == 0)
          pending_freeheader(job, &hdr);
      }
      else {
        buf = gk_cmalloc(msize, "bcast: buf");
        BDASSERT(MPI_Bcast(buf, msg->count, mpi_dt(msg->datatype), 
                     babel_get_node(comm, msg->source), comm->mpi_comm)
            == MPI_SUCCESS);

        pending_addbcast(job, msg, buf, msize, comm->lsize);
      }

      BD_GET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */
    }

    /* all processes have called the bcast, unblock them in order to get to 
       the second step of bcast */
    comm->counter = comm->lsize;
    for (i=0; i<comm->lsize; i++)
      slvpool_cunblock(job, comm->sranks[i]);
  }

  BD_LET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

  gk_free((void **)&arg, LTERM);

  return NULL;
}


/*************************************************************************/
/*! Response to a BDMPI_Bcast - Recv part.
    Protocol:
      - Finds the header of the bcast
      - Reads and copies the data to the slave.
    Note:
      - By construction the header of the bcast should be there! If not,
        then somerhing went very bad.
*/
/*************************************************************************/
void *mstr_bcast_recv(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  int srank, counter, response;
  header_t *hdr;
  bdmcomm_t *comm;

  BD_GET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] mstr_bacst_recv: entering: msg->myrank: %d.\n", 
        job->mynode, msg->myrank));

  /* hook to the key info */
  comm  = job->comms[msg->mcomm];
  srank = babel_get_srank(comm, msg->myrank);

  /* get the header */
  if ((hdr = pending_getbcast(job, msg, 1, &counter)) == NULL)
    slvpool_abort(1, "mstr_bcast_recv: failed to find a header for a bcast operation! [%d %d %d]\n", 
        msg->source, msg->dest, (int)msg->mcomm);

  /* send the data */
  xfer_out_scb(job->scbs[srank], hdr->buf, hdr->msg.count, hdr->msg.datatype); 

  BD_LET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

  /* wait for the ACK from the slave */
  bdmq_recv(job->c2mMQs[srank], &response, sizeof(int));

  if (counter == 0) { 
    BD_GET_WRLOCK(job->comms[msg->mcomm]->rwlock);
    pending_freeheader(job, &hdr);
    BD_LET_WRLOCK(job->comms[msg->mcomm]->rwlock);
  }

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] mstr_bacst_recv: exiting: msg->myrank: %d.\n", 
        job->mynode, msg->myrank));

  gk_free((void **)&arg, LTERM);

  return NULL;
}

