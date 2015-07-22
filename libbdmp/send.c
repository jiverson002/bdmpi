/*
 * send.c
 *
 * Implements the various send operations
 *
 * Started 4/4/2013
 * George
 *
 */

#include "bdmplib.h"


/*************************************************************************/
/* Performs a blocking send operation */
/*************************************************************************/
int bdmp_Send(sjob_t *job, void *buf, size_t count, BDMPI_Datatype datatype,
          int dest, int tag, BDMPI_Comm comm)
{
  bdmsg_t msg;
  int mype, response;

  S_IFSET(BDMPI_DBG_IPCS,
      bdprintf("BDMPI_Send: sending to %d [%zu] [goMQlen: %d]\n", dest, count,
        bdmq_length(job->goMQ)));

  /* some error checking */
  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }
  if (!datatype_isvalid(datatype)) {
    fprintf(stderr, "The datatype is invalid.\n");
    return BDMPI_ERR_TYPE;
  }
  if (dest == BDMPI_PROC_NULL) {
    return BDMPI_SUCCESS;
  } else if (dest < 0 || dest >= comm->size) {
    fprintf(stderr, "The dest rank is invalid.\n");
    return BDMPI_ERR_RANK;
  }

  mype = comm->rank;

  memset(&msg, 0, sizeof(bdmsg_t));
  msg.msgtype  = BDMPI_MSGTYPE_SEND;
  msg.mcomm    = comm->mcomm;
  msg.myrank   = mype;
  msg.fnum     = (bdmp_msize(count, datatype) > job->smallmsg ? xfer_getfnum() : -1);
  msg.tag      = tag;
  msg.source   = mype;
  msg.dest     = dest;
  msg.count    = count;
  msg.datatype = datatype;

  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Send: Copying %zu counts.\n", count));

  if (msg.fnum == -1) { /* notify & copy */
    bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));
    xfer_out_scb(job->scb, buf, count, datatype);
  }
  else { /* store & notify (remove race condition) */
    xfer_out_disk(msg.fnum, buf, count, datatype);
    bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));
  }

  /* wait for the master to tell us that he is done */
  bdmq_recv(job->c2sMQ, &response, sizeof(int));

  return BDMPI_SUCCESS;
}


/*************************************************************************/
/* Performs a non-blocking send operation */
/*************************************************************************/
int bdmp_Isend(sjob_t *job, void *buf, size_t count, BDMPI_Datatype datatype,
          int dest, int tag, BDMPI_Comm comm, BDMPI_Request *r_request)
{
  bdrequest_t *request;

  if (dest == BDMPI_PROC_NULL) {
    *r_request = BDMPI_REQUEST_NULL;
    return BDMPI_SUCCESS;
  }

  request = *r_request = (bdrequest_t *)bd_malloc(sizeof(bdrequest_t), "request");
  memset(request, 0, sizeof(bdrequest_t));

  request->type  = BDMPI_REQUEST_ISEND;
  request->state = bdmp_Send(job, buf, count, datatype, dest, tag, comm);

  request->status.BDMPI_SOURCE = comm->rank;
  request->status.BDMPI_TAG    = tag;
  request->status.BDMPI_ERROR  = request->state;
  request->status.MPI_SOURCE   = request->status.BDMPI_SOURCE;
  request->status.MPI_TAG      = request->status.BDMPI_TAG;
  request->status.MPI_ERROR    = request->status.BDMPI_ERROR;
  request->status.comm         = comm;
  request->status.count        = count;
  request->status.datatype     = datatype;

  return BDMPI_SUCCESS;
}
