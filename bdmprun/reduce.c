/*!
\file
\brief Various functions for performing reduction operations.
\date Started 4/6/2013
\author George
*/


#include "bdmprun.h"



/*************************************************************************/
/*! Response to a BDMPI_Reduce. - Send-to-master part
    Protocol:
       Blocks the process.
       Decreases the counter associated with the counter.
       Reduces the data, if myrank!=root
       If counter becomes 0, then moves all processes to runnable state
       and sets the counter back to the size of that communicator.

    Meaning of fields of msg:
       msg->dest is the root of the reduction (ie., the slave that sent 
                 a message to the master)
*/
/*************************************************************************/
void *mstr_reduce_send(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  size_t i, dtsize, len, chunk;
  int srank, sleeping;
  header_t *hdr;
  char *buf=NULL, *rbuf=NULL;
  bdmcomm_t *comm;

  BD_GET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  /* hook to the key info */
  comm  = job->comms[msg->mcomm];
  srank = babel_get_srank(comm, msg->myrank);

  /* compute basic stats of the data size */
  dtsize = bdmp_sizeof(msg->datatype);
  chunk  = job->jdesc->smsize/dtsize;

  /* deal with copid */
  if (comm->counter == comm->lsize)
    comm->copid++;
  if (msg->myrank == msg->dest) {
    if (bdmq_send(job->c2sMQs[srank], &(comm->copid), sizeof(int)) == -1)
      bdprintf("Failed to send a message to %d: %s\n", srank, strerror(errno));
  }
  msg->copid = comm->copid;


  /* check to see if this is the first slave calling the reduction */
  if (comm->counter == comm->lsize) {
    /* allocate memory for the data on the master */
    buf = gk_cmalloc(msg->count*dtsize, "reduce: buf");

    /* get into a loop copying and storing the data */
    for (i=0; i<msg->count; i+=chunk) {
      len = (i+chunk < msg->count ? chunk : msg->count - i);
      bdscb_wait_full(job->scbs[srank]);
      memcpy(buf+i*dtsize, job->scbs[srank]->buf, len*dtsize);
      bdscb_post_empty(job->scbs[srank]);
    }

    pending_addreduce(job, msg, buf, msg->count*dtsize,
        (babel_is_local(comm, msg->dest) ? comm->lsize+1 : comm->lsize));
  }
  else { /* this is a set of data that needs to be reduced with the current data */
    hdr = pending_getreduce(job, msg);

    if (hdr == NULL) {
      /* this should not have happened */
      slvpool_abort(1, "Failed to find a header for a reduce operation! [%d %d %d]\n", 
          msg->source, msg->dest, msg->mcomm);
    }
    else {
      /* get into a loop copying, reducing, and updating data */
      for (i=0; i<msg->count; i+=chunk) {
        len = (i+chunk < msg->count ? chunk : msg->count - i);
        bdscb_wait_full(job->scbs[srank]);
        gk_startwctimer(job->aux3Tmr);
        reduce_op(hdr->buf+i*dtsize, job->scbs[srank]->buf, len, msg->datatype, msg->op);
        gk_stopwctimer(job->aux3Tmr);
        bdscb_post_empty(job->scbs[srank]);
      }
    }
  }

  /* block the slave, once you get the ok from it */
  xfer_in_scb(job->scbs[srank], &sleeping, sizeof(int), BDMPI_BYTE);
  slvpool_cblock(job, srank);

  /* do some work when all slaves have called */
  if (--comm->counter == 0) {
    hdr = pending_getreduce(job, msg);
    if (comm->nnodes > 1) { /* do something if more than one node is involved */
      BD_LET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

      if (babel_is_local(comm, msg->dest))
        rbuf = gk_cmalloc(msg->count*dtsize, "reduce: rbuf");

      MPI_Reduce(hdr->buf, rbuf, msg->count, mpi_dt(msg->datatype), mpi_op(msg->op),
          babel_get_node(comm, msg->dest), comm->mpi_comm);

      if (babel_is_local(comm, msg->dest)) {
        memcpy(hdr->buf, rbuf, msg->count*dtsize);
        gk_free((void **)&rbuf, LTERM);
      }

      BD_GET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */
    }
    if (hdr->counter == 0)
      pending_freeheader(job, &hdr);

    /* all processes have called reduce, unblock them in order to get to the second 
       step of reduce */
    comm->counter = comm->lsize;
    for (i=0; i<comm->lsize; i++)
      slvpool_cunblock(job, comm->sranks[i]);
  }

  BD_LET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

  gk_free((void **)&arg, LTERM);

  return NULL;
}


/*************************************************************************/
/*! Response to a BDMPI_Reduce. - Recv-from-root part
    Protocol:
       Sends the data to the root.

    Meaning of fields of msg:
       msg->dest is the root of the reduction (ie., the slave that sent 
                 a message to the master)

    Note: Only the root will ever call this.                   
*/
/*************************************************************************/
void *mstr_reduce_recv(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  int srank;
  header_t *hdr;
  bdmcomm_t *comm;

  BD_GET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  if (msg->myrank != msg->dest) 
    slvpool_abort(1, "The reduce_recv is not called from the root: root:%d myrank:%d\n",
        msg->dest, msg->myrank);

  /* hook to the key info */
  comm  = job->comms[msg->mcomm];
  srank = babel_get_srank(comm, msg->myrank);

  /* get the header */
  hdr = pending_getreduce(job, msg);

  /* send the data to the slave */
  xfer_out_scb(job->scbs[srank], hdr->buf, msg->count, msg->datatype); 

  pending_freeheader(job, &hdr);

  BD_LET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

  gk_free((void **)&arg, LTERM);

  return NULL;
}


/*************************************************************************/
/*! Response to a BDMPI_Allreduce. - Send-to-master part
    Protocol:
       Blocks the process.
       Decreases the counter associated with the counter.
       Reduces the data 
       If counter becomes 0, then moves all processes to runnable state
       and sets the counter back to the size of that communicator.

    Meaning of fields of msg:
       msg->dest is 0 which acts as the root of the virtual reduction.
*/
/*************************************************************************/
void *mstr_allreduce_send(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  size_t i, dtsize, len, chunk;
  int srank, sleeping;
  header_t *hdr;
  char *buf=NULL, *rbuf=NULL;
  bdmcomm_t *comm;

  BD_GET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  /* hook to the key info */
  comm  = job->comms[msg->mcomm];
  srank = babel_get_srank(comm, msg->myrank);

  if (comm->counter == comm->lsize)
    comm->copid++;
  if (bdmq_send(job->c2sMQs[srank], &(comm->copid), sizeof(int)) == -1)
    bdprintf("Failed to send a message to %d: %s\n", srank, strerror(errno));
  msg->copid = comm->copid;

  /* compute basic stats of the data size */
  dtsize = bdmp_sizeof(msg->datatype);
  chunk  = job->jdesc->smsize/dtsize;

  /* check to see if this is the first slave calling the reduction */
  if (comm->counter == comm->lsize) {
    /* allocate memory on the master */
    buf = gk_cmalloc(msg->count*dtsize, "allreduce: buf");

    /* get into a loop copying and storing the data */
    for (i=0; i<msg->count; i+=chunk) {
      len = (i+chunk < msg->count ? chunk : msg->count - i);
      bdscb_wait_full(job->scbs[srank]);
      memcpy(buf+i*dtsize, job->scbs[srank]->buf, len*dtsize);
      bdscb_post_empty(job->scbs[srank]);
    }

    pending_addallreduce(job, msg, buf, msg->count*dtsize);
  }
  else { /* this is a set of data that needs to be reduced with the current data */
    hdr = pending_getallreduce(job, msg, NULL);

    if (hdr == NULL) {
      /* this should not have happened */
      slvpool_abort(1, "Failed to find a header for an allreduce operation! [%d %d]\n", 
          msg->myrank, msg->mcomm);
    }
    else {
      /* get into a loop copying, reducing, and updating data */
      for (i=0; i<msg->count; i+=chunk) {
        len = (i+chunk < msg->count ? chunk : msg->count - i);
        bdscb_wait_full(job->scbs[srank]);
        reduce_op(hdr->buf+i*dtsize, job->scbs[srank]->buf, len, msg->datatype, msg->op);
        bdscb_post_empty(job->scbs[srank]);
      }
    }
  }

  /* block the slave, once you get the ok from it */
  xfer_in_scb(job->scbs[srank], &sleeping, sizeof(int), BDMPI_BYTE);
  slvpool_cblock(job, srank);

  /* take action when all slaves have called */
  if (--comm->counter == 0) {
    hdr = pending_getallreduce(job, msg, NULL);
    if (comm->nnodes > 1) { /* do something if more than one node is involved */
      BD_LET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

      rbuf = gk_cmalloc(msg->count*dtsize, "reduce: rbuf");

      MPI_Allreduce(hdr->buf, rbuf, msg->count, mpi_dt(msg->datatype), mpi_op(msg->op),
          comm->mpi_comm);

      memcpy(hdr->buf, rbuf, msg->count*dtsize);
      gk_free((void **)&rbuf, LTERM);

      BD_GET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */
    }

    /* all processes have called reduce, unblock them in order to get to the second 
       step of reduce */
    comm->counter = comm->lsize;
    for (i=0; i<comm->lsize; i++)
      slvpool_cunblock(job, comm->sranks[i]);
  }

  BD_LET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

  gk_free((void **)&arg, LTERM);

  return NULL;
}


/*************************************************************************/
/*! Response to a BDMPI_Allreduce - Send-to-all-slaves part.
    Protocol:
      - Finds the header of the allreduce 
      - Reads and copies the data to the slave.
    Note:
      - By construction the header of the allreduce should be there! If not,
        then somerhing went very bad.
*/
/*************************************************************************/
void *mstr_allreduce_recv(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  int srank, counter, response;
  header_t *hdr;
  bdmcomm_t *comm;

  BD_GET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  /* hook to the key info */
  comm  = job->comms[msg->mcomm];
  srank = babel_get_srank(comm, msg->myrank);

  /* see if the allreduce has been posted */
  if ((hdr = pending_getallreduce(job, msg, &counter)) == NULL)
    slvpool_abort(1, "Failed to find a header for a allreduce operation! [%d %d]\n", 
        msg->myrank, (int)msg->mcomm);

  /* send the data to the slave */
  xfer_out_scb(job->scbs[srank], hdr->buf, msg->count, msg->datatype); 

  BD_LET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

  /* wait for the ACK from the slave */
  bdmq_recv(job->c2mMQs[srank], &response, sizeof(int));

  if (counter == 0) { 
    BD_GET_WRLOCK(job->comms[msg->mcomm]->rwlock);
    pending_freeheader(job, &hdr);
    BD_LET_WRLOCK(job->comms[msg->mcomm]->rwlock);
  }

  gk_free((void **)&arg, LTERM);

  return NULL;
}

