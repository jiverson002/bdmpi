/*!
\file
\brief Implements the bcast operation.
\date Started 4/12/2013
\author George
*/


#include "bdmplib.h"


/*************************************************************************/
/* Performs BDMPI_Allgatherv() */
/*************************************************************************/
int bdmp_Allgatherv(sjob_t *job,
          void *sendbuf, size_t sendcount, BDMPI_Datatype sendtype,
          void *recvbuf, size_t *recvcounts,  size_t *displs,
          BDMPI_Datatype recvtype, BDMPI_Comm comm)
{
  int p, response, sleeping=1;
  bdmsg_t msg, rmsg, gomsg;
  size_t rdtsize;

  S_IFSET(BDMPI_DBG_IPCS,
      bdprintf("BDMPI_Allgather: entering: comm: %p [goMQlen: %d]\n",
          comm, bdmq_length(job->goMQ)));

  /* some error checking */
  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }
  if (!datatype_isvalid(sendtype)) {
    fprintf(stderr, "The sendtype is invalid.\n");
    return BDMPI_ERR_TYPE;
  }
  if (!datatype_isvalid(recvtype)) {
    fprintf(stderr, "The recvtype is invalid.\n");
    return BDMPI_ERR_TYPE;
  }

  rdtsize = bdmp_sizeof(recvtype);

  /* notify the master that you entering an allgather */
  msg.msgtype  = BDMPI_MSGTYPE_ALLGATHERI;
  msg.mcomm    = comm->mcomm;
  msg.myrank   = comm->rank;
  msg.count    = sendcount;
  msg.datatype = sendtype;
  msg.source   = comm->rank;
  msg.dest     = comm->rank;

  bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

  /* get the copid from the master */
  bdmq_recv(job->c2sMQ, &msg.copid, sizeof(int));

  /* everybody sends their data to the master */
  if (bdmp_msize(sendcount, sendtype) <= job->smallmsg) { /* to memory */
    msg.fnum = -1;

    /* tell the master that we are sending it some data */
    xfer_out_scb(job->scb, &msg, sizeof(bdmsg_t), BDMPI_BYTE);

    /* send the data */
    xfer_out_scb(job->scb, sendbuf, sendcount, sendtype);
  }
  else { /* to disk */
    msg.fnum = xfer_getfnum();

    /* save the data to the disk */
    xfer_out_disk(msg.fnum, sendbuf, sendcount, sendtype);

    /* tell the master that we are storing data on disk */
    xfer_out_scb(job->scb, &msg, sizeof(bdmsg_t), BDMPI_BYTE);
  }

  /* prepare to go to sleep */
  S_SB_IFSET(BDMPI_SB_SAVEALL) {
    if (job->jdesc->nr < job->jdesc->ns)
      sb_saveall();
  }
  xfer_out_scb(job->scb, &sleeping, sizeof(int), BDMPI_BYTE);

  /* go to sleep until everybody has called the allgather */
  BDMPL_SLEEP(job, gomsg);

  /* notify the master that you want to receive the allgathered data */
  msg.msgtype  = BDMPI_MSGTYPE_ALLGATHERF;
  msg.datatype = recvtype;
  bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

  for (p=0; p<comm->size; p++) {
    /* get size info of the received message */
    xfer_in_scb(job->scb, &rmsg, sizeof(bdmsg_t), BDMPI_BYTE);

    if (bdmp_msize(rmsg.count, rmsg.datatype) != bdmp_msize(recvcounts[p], recvtype))
      errexit("BDMPI_Allgatherv: Send/Recv [%d => %d] buffers are of different sizes: %zu %zu\n",
          p, comm->rank, bdmp_msize(sendcount, sendtype), bdmp_msize(recvcounts[p], recvtype));

    /* get the data */
    if (rmsg.fnum == -1)
      xfer_in_scb(job->scb, (char *)recvbuf+displs[p]*rdtsize, rmsg.count, rmsg.datatype);
    else
      xfer_in_disk(rmsg.fnum, (char *)recvbuf+displs[p]*rdtsize, rmsg.count, rmsg.datatype, 0);

    /* tell the master that you got the data */
    bdmq_send(job->c2mMQ, &response, sizeof(int));
  }

  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Allgather: exiting: comm: %p\n", comm));

  return BDMPI_SUCCESS;
}
