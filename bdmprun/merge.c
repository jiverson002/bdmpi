/*!
\file
\brief Various functions for performing sparse reduction operations.
\date Started 4/6/2013
\author George
*/


#include "bdmprun.h"


/*************************************************************************/
/*! Response to a BDMPI_Merge. - Send-to-master part
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
void *mstr_merge_send(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  size_t i, dtsize, len, chunk;
  int srank, ndest, nleft, sleeping;
  header_t *hdr;
  char *buf=NULL, *rbuf=NULL;
  mergeinfo_t *minfo=NULL, *cminfo;
  bdmcomm_t *comm;
  MPI_Status status;

  BD_GET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_merge_send: counter: %d [entering]\n",
                        job->mynode, msg->myrank, job->comms[msg->mcomm]->counter));

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


  /* get the data */
  minfo = (mergeinfo_t *)gk_malloc(sizeof(mergeinfo_t), "minfo");
  memset(minfo, 0, sizeof(mergeinfo_t));
  minfo->len = msg->count;
  minfo->buf = gk_cmalloc(msg->count*dtsize, "merge: buf");
  minfo->ids = gk_imalloc(minfo->len, "merge: ids");

  gk_startwctimer(job->aux2Tmr);
  xfer_in_scb(job->scbs[srank], minfo->buf, msg->count, msg->datatype);
  xfer_in_scb(job->scbs[srank], minfo->ids, msg->count, BDMPI_INT);
  gk_stopwctimer(job->aux2Tmr);

  /* block the slave, once you get the ok from it */
  xfer_in_scb(job->scbs[srank], &sleeping, sizeof(int), BDMPI_BYTE);
  slvpool_cblock(job, srank);

  /* check to see if this is the first slave calling the reduction */
  if (comm->counter == comm->lsize) {
    pending_addmerge(job, msg, minfo, 
        (babel_is_local(comm, msg->dest) ? comm->lsize+1 : comm->lsize));
  }
  else { /* this is a set of data that needs to be reduced with the current data */
    if ((hdr = pending_getmerge(job, msg)) == NULL) {
      slvpool_abort(1, "Failed to find a header for a merge operation! [%d %d %d]\n", 
          msg->source, msg->dest, msg->mcomm);
    }

    cminfo = (mergeinfo_t *)hdr->buf;
    gk_startwctimer(job->aux1Tmr);
    merge_infos(cminfo, minfo, msg->datatype, msg->op);
    gk_stopwctimer(job->aux1Tmr);
    gk_free((void **)&minfo->buf, &minfo->ids, &minfo, LTERM);
  }


  /* do some work when all slaves have called */
  if (--comm->counter == 0) {
    if ((hdr = pending_getmerge(job, msg)) == NULL) 
      slvpool_abort(1, "Failed to find a header for a merge operation! [%d %d %d]\n", 
          msg->source, msg->dest, msg->mcomm);
    cminfo = (mergeinfo_t *)hdr->buf;

    if (comm->nnodes > 1) { /* do something if more than one node is involved */
      BD_LET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

      ndest = babel_get_node(comm, msg->dest);

      if (comm->mynode == ndest) {
        minfo = (mergeinfo_t *)gk_malloc(sizeof(mergeinfo_t), "minfo");
        memset(minfo, 0, sizeof(mergeinfo_t));

        nleft = comm->nnodes-1;
        while (nleft > 0) {
          MPI_Recv(&(minfo->len), sizeof(size_t), MPI_BYTE, MPI_ANY_SOURCE, 20000, 
              comm->mpi_comm, &status);
          i = status.MPI_SOURCE;

          minfo->buf = gk_cmalloc(minfo->len*dtsize, "merge: buf");
          minfo->ids = gk_imalloc(minfo->len, "merge: ids");

          MPI_Recv(minfo->buf, minfo->len*dtsize, MPI_BYTE, i, 20001, comm->mpi_comm, 
              &status);
          MPI_Recv(minfo->ids, minfo->len, MPI_INT, i, 20002, comm->mpi_comm, 
              &status);

          gk_startwctimer(job->aux3Tmr);
          merge_infos(cminfo, minfo, msg->datatype, msg->op);
          gk_stopwctimer(job->aux3Tmr);
          gk_free((void **)&minfo->buf, &minfo->ids, LTERM);
          nleft--;
        }

        gk_free((void **)&minfo, LTERM);
      }
      else {
        MPI_Send(&(cminfo->len), sizeof(size_t), MPI_BYTE, ndest, 20000, comm->mpi_comm);
        MPI_Send(cminfo->buf, cminfo->len*dtsize, MPI_BYTE, ndest, 20001, comm->mpi_comm);
        MPI_Send(cminfo->ids, cminfo->len, MPI_INT, ndest, 20002, comm->mpi_comm);
      }
      MPI_Barrier(comm->mpi_comm);

      BD_GET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */
    }

    if (hdr->counter == 0) {
      gk_free((void **)&cminfo->ids, &cminfo->buf, LTERM);
      pending_freeheader(job, &hdr);
    }

    /* all processes have called reduce, unblock them in order to get to the second 
       step of reduce */
    comm->counter = comm->lsize;
    for (i=0; i<comm->lsize; i++)
      slvpool_cunblock(job, comm->sranks[i]);
  }

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_merge_send: counter: %d [exiting]\n",
                        job->mynode, msg->myrank, job->comms[msg->mcomm]->counter));

  BD_LET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

  gk_free((void **)&arg, LTERM);

  return NULL;
}


/*************************************************************************/
/*! Response to a BDMPI_Merge. - Recv-from-root part
    Protocol:
       Sends the data to the root.

    Meaning of fields of msg:
       msg->dest is the root of the reduction (ie., the slave that sent 
                 a message to the master)

    Note: Only the root will ever call this.                   
*/
/*************************************************************************/
void *mstr_merge_recv(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  int srank;
  header_t *hdr;
  bdmcomm_t *comm;
  mergeinfo_t *cminfo;

  BD_GET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_merge_recv: counter: %d [entering]\n",
                        job->mynode, msg->myrank, job->comms[msg->mcomm]->counter));

  if (msg->myrank != msg->dest) 
    slvpool_abort(1, "The merge_recv is not called from the root: root:%d myrank:%d\n",
        msg->dest, msg->myrank);

  /* hook to the key info */
  comm  = job->comms[msg->mcomm];
  srank = babel_get_srank(comm, msg->myrank);

  /* get the header */
  if ((hdr = pending_getmerge(job, msg)) == NULL) 
    slvpool_abort(1, "Failed to find a header for a merge operation! [%d %d %d]\n", 
        msg->source, msg->dest, msg->mcomm);

  cminfo = (mergeinfo_t *)hdr->buf;

  if (msg->count < cminfo->len) 
    slvpool_abort(1, "The length of the merged info is greater than the provided recv buffer [%zu %zu].\n",
        msg->count, cminfo->len);

  /* send the data to the slave */
  xfer_out_scb(job->scbs[srank], &(cminfo->len), 1, BDMPI_SIZE_T); 
  xfer_out_scb(job->scbs[srank], cminfo->buf, cminfo->len, msg->datatype); 
  xfer_out_scb(job->scbs[srank], cminfo->ids, cminfo->len, BDMPI_INT); 

  gk_free((void **)&cminfo->ids, &cminfo->buf, LTERM);
  pending_freeheader(job, &hdr);

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_merge_recv: counter: %d [exiting]\n",
                        job->mynode, msg->myrank, job->comms[msg->mcomm]->counter));

  BD_LET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

  gk_free((void **)&arg, LTERM);

  return NULL;
}



/*************************************************************************/
/*! Performs a merge between a and b and returns the result in a */
/*************************************************************************/
void merge_infos(mergeinfo_t *a, mergeinfo_t *b, BDMPI_Datatype dtype, BDMPI_Op op)
{
  size_t i, j, k, maxlen, alen, blen;
  int *ids, *aids, *bids;
  double *buf, *abuf, *bbuf;

  if (dtype != BDMPI_DOUBLE && op != BDMPI_SUM)
    slvpool_abort(1, "Only the BDMPI_DOUBLE datatype and BDMPI_SUM operation has been implemented.\n");

  alen = a->len;
  aids = a->ids;
  abuf = (double *)a->buf;

  blen = b->len;
  bids = b->ids;
  bbuf = (double *)b->buf;

  maxlen = a->len+b->len;
  ids    = gk_imalloc(maxlen, "ids");
  buf    = gk_dmalloc(maxlen, "buf");

  for (i=0, j=0, k=0; i<alen && j<blen; k++) {
    if (aids[i] == bids[j]) {
      ids[k] = aids[i];
      buf[k] = abuf[i] + bbuf[j];
      i++; j++;
    }
    else if (aids[i] < bids[j]) {
      ids[k] = aids[i];
      buf[k] = abuf[i];
      i++;
    }
    else {
      ids[k] = bids[j];
      buf[k] = bbuf[j];
      j++;
    }
  }

  for (; i<alen; i++, k++) {
    ids[k] = aids[i];
    buf[k] = abuf[i];
  }

  for (; j<blen; j++, k++) {
    ids[k] = bids[j];
    buf[k] = bbuf[j];
  }

  gk_free((void **)&a->ids, &a->buf, LTERM);
  a->len = k;
  a->ids = ids;
  a->buf = buf;

  return;
}
