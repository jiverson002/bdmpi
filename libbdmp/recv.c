/*!
\file
\brief Implements the various recv operations.
\date Started 4/6/2013
\author George
*/


#include "bdmplib.h"


/*************************************************************************/
/* Performs a blocking recv operation */
/*************************************************************************/
int bdmp_Recv(sjob_t *job, void *buf, size_t count, BDMPI_Datatype datatype,
          int source, int tag, BDMPI_Comm comm, BDMPI_Status *status)
{
  int mype, response, flag=0;
  bdmsg_t msg, rmsg, gomsg;

  S_IFSET(BDMPI_DBG_IPCS,
      bdprintf("BDMPI_Recv: Receiving from %d [goMQlen: %d]\n", source, bdmq_length(job->goMQ)));

  /* some error checking */
  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }
  if (!datatype_isvalid(datatype)) {
    fprintf(stderr, "The datatype is invalid.\n");
    return BDMPI_ERR_TYPE;
  }
  if (source == BDMPI_PROC_NULL) {
    return BDMPI_SUCCESS;
  } else if (source != BDMPI_ANY_SOURCE && (source < 0 || source >= comm->size)) {
    fprintf(stderr, "The source rank is invalid.\n");
    return BDMPI_ERR_RANK;
  }

  mype = comm->rank;

  /* save the state in case you go to sleep */
  bdmp_Iprobe(job, source, tag, comm, &flag, BDMPI_STATUS_IGNORE);
  if (flag == 0) {
    S_SB_IFSET(BDMPI_SB_DISCARD) {
      sbma_mclear(buf, bdmp_msize(count, datatype));
    }
    S_SB_IFSET(BDMPI_SB_SAVEALL) {
      if (job->jdesc->nr < job->jdesc->ns)
        sbma_mevictall();
    }
  }

  memset(&msg, 0, sizeof(bdmsg_t));
  msg.msgtype  = BDMPI_MSGTYPE_RECV;
  msg.mcomm    = comm->mcomm;
  msg.myrank   = mype;
  msg.tag      = tag;
  msg.source   = source;
  msg.dest     = mype;
  msg.count    = count;
  msg.datatype = datatype;

  for (;;) {
    /* notify the master that you want to receive a message */
    bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

    /* get the master's response  */
    bdmq_recv(job->c2sMQ, &response, sizeof(int));

    S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Recv: Response from server: %d\n", response));

    if (response == 1)
      break;  /* we got the go-ahead */

    /* prepare to go to sleep */
    S_SB_IFSET(BDMPI_SB_SAVEALL) {
      if (job->jdesc->nr < job->jdesc->ns)
        sbma_mevictall();
    }

    /* go to sleep... */
    BDMPL_SLEEP(job, gomsg, 1);
  }

  /* get the missing message info from the master */
  xfer_in_scb(job->scb, &rmsg, sizeof(bdmsg_t), BDMPI_BYTE);

  /* check that we have enough buffer space */
  if (bdmp_msize(rmsg.count, rmsg.datatype) > bdmp_msize(count, datatype))
    errexit("BDMPI_Recv: Received message is greater than provided buffer: recv: %zu; buffer: %zu\n",
        bdmp_msize(rmsg.count, rmsg.datatype), bdmp_msize(count, datatype));

  /* copy the data */
  if (rmsg.fnum == -1)
    xfer_in_scb(job->scb, buf, rmsg.count, rmsg.datatype);
  else
    xfer_in_disk(rmsg.fnum, buf, rmsg.count, rmsg.datatype, 1);

  if (status != BDMPI_STATUS_IGNORE) {
    status->BDMPI_SOURCE = rmsg.source;
    status->BDMPI_TAG    = rmsg.tag;
    status->BDMPI_ERROR  = BDMPI_SUCCESS;
    status->MPI_SOURCE   = status->BDMPI_SOURCE;
    status->MPI_TAG      = status->BDMPI_TAG;
    status->MPI_ERROR    = status->BDMPI_ERROR;
    status->count        = rmsg.count;
  }

  return BDMPI_SUCCESS;
}


/*************************************************************************/
/* Performs a non-blocking recv operation */
/*************************************************************************/
int bdmp_Irecv(sjob_t *job, void *buf, size_t count, BDMPI_Datatype datatype,
          int source, int tag, BDMPI_Comm comm, BDMPI_Request *r_request)
{
  int mype, response;
  bdmsg_t msg, rmsg;
  bdrequest_t *request;

  S_IFSET(BDMPI_DBG_IPCS,
      bdprintf("BDMPI_Irecv: Receiving from %d [goMQlen: %d]\n", source, bdmq_length(job->goMQ)));

  /* some error checking */
  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }
  if (!datatype_isvalid(datatype)) {
    fprintf(stderr, "The datatype is invalid.\n");
    return BDMPI_ERR_TYPE;
  }
  if (source == BDMPI_PROC_NULL) {
    *r_request = BDMPI_REQUEST_NULL;
    return BDMPI_SUCCESS;
  } else if (source != BDMPI_ANY_SOURCE && (source < 0 || source >= comm->size)) {
    fprintf(stderr, "The source rank is invalid.\n");
    return BDMPI_ERR_RANK;
  }

  request = *r_request = (bdrequest_t *)bd_malloc(sizeof(bdrequest_t), "request");
  memset(request, 0, sizeof(bdrequest_t));
  request->type = BDMPI_REQUEST_IRECV;

  mype = comm->rank;

  memset(&msg, 0, sizeof(bdmsg_t));
  msg.msgtype  = BDMPI_MSGTYPE_IRECV;
  msg.mcomm    = comm->mcomm;
  msg.myrank   = mype;
  msg.tag      = tag;
  msg.source   = source;
  msg.dest     = mype;
  msg.count    = count;
  msg.datatype = datatype;


  /* notify the master that you want to receive a message */
  bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

  /* get the master's response  */
  bdmq_recv(job->c2sMQ, &response, sizeof(int));

  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Recv: Response from server: %d\n", response));

  if (response == 1) { /* we got the go-ahead */
    /* get the missing message info from the master */
    xfer_in_scb(job->scb, &rmsg, sizeof(bdmsg_t), BDMPI_BYTE);

    /* check that we have enough buffer space */
    if (bdmp_msize(rmsg.count, rmsg.datatype) > bdmp_msize(count, datatype))
      errexit("BDMPI_Irecv: Received message is greater than provided buffer: recv: %zu; buffer: %zu\n",
          bdmp_msize(rmsg.count, rmsg.datatype), bdmp_msize(count, datatype));

    /* copy the data */
    if (rmsg.fnum == -1)
      xfer_in_scb(job->scb, buf, rmsg.count, rmsg.datatype);
    else
      xfer_in_disk(rmsg.fnum, buf, rmsg.count, rmsg.datatype, 1);

    request->state               = BDMPI_SUCCESS;
    request->status.BDMPI_SOURCE = rmsg.source;
    request->status.BDMPI_TAG    = rmsg.tag;
    request->status.BDMPI_ERROR  = BDMPI_SUCCESS;
    request->status.MPI_SOURCE   = request->status.BDMPI_SOURCE;
    request->status.MPI_TAG      = request->status.BDMPI_TAG;
    request->status.MPI_ERROR    = request->status.BDMPI_ERROR;
    request->status.count        = rmsg.count;
    request->status.datatype     = rmsg.datatype;
  }
  else {
    request->state               = BDMPI_INPROGRESS;
    request->buf                 = buf;
    request->status.comm         = comm;
    request->status.BDMPI_SOURCE = source;
    request->status.BDMPI_TAG    = tag;
    request->status.MPI_SOURCE   = request->status.BDMPI_SOURCE;
    request->status.MPI_TAG      = request->status.BDMPI_TAG;
    request->status.count        = count;
    request->status.datatype     = datatype;
  }

  return BDMPI_SUCCESS;
}
