/*!
\file
\brief Implements the scatter operation.
\date Started 5/3/2013
\author George
*/


#include "bdmplib.h"


/*************************************************************************/
/* Performs BDMPI_Scatterv() when the communicator is within a single node */
/*************************************************************************/
int bdmp_Scatterv_node(sjob_t *job,
          void *sendbuf, size_t *sendcounts, size_t *sdispls, BDMPI_Datatype sendtype,
          void *recvbuf, size_t recvcount, BDMPI_Datatype recvtype, int root,
          BDMPI_Comm comm)
{
  size_t size, sdtsize, rdtsize;
  int npes, mype, p, sleeping=1;
  bdmsg_t msg, rmsg, gomsg;


  S_IFSET(BDMPI_DBG_IPCS,
      bdprintf("BDMPI_Scatterv_node: entering: comm: %p [goMQlen: %d]\n",
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
  if (root < 0 || root >= comm->size) {
    fprintf(stderr, "The root rank is invalid.\n");
    return BDMPI_ERR_ROOT;
  }

  npes = comm->size;
  mype = comm->rank;

  sdtsize = bdmp_sizeof(sendtype);
  rdtsize = bdmp_sizeof(recvtype);


  /* notify the master that you entering a scatter */
  msg.msgtype  = BDMPI_MSGTYPE_SCATTERI;
  msg.mcomm    = comm->mcomm;
  msg.myrank   = mype;
  msg.source   = root;
  msg.datatype = sendtype;

  bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

  /* get the copid from the master */
  bdmq_recv(job->c2sMQ, &msg.copid, sizeof(int));


  /* start sending data to the master or writing them to the disk */
  if (mype == root) {
    for (p=0; p<npes; p++) {
      if (p == mype) {
        if (sendcounts[p]*sdtsize > recvcount*rdtsize)
          errexit("BDMPI_Scatterv: Amount sent is greater than size of recv buffer.\n");
        memcpy((char *)recvbuf, (char *)sendbuf+sdispls[p]*sdtsize, sendcounts[p]*sdtsize);
      }
      else {
        msg.dest  = p;
        msg.count = sendcounts[p];

        if (sendcounts[p]*sdtsize <= job->smallmsg) { /* to memory */
          /* tell the master that we are sending it some data */
          msg.fnum = -1;
          xfer_out_scb(job->scb, &msg, sizeof(bdmsg_t), BDMPI_BYTE);

          /* copy the data */
          xfer_out_scb(job->scb, (char *)sendbuf+sdispls[p]*sdtsize, sendcounts[p], sendtype);
        }
        else { /* to disk */
          msg.fnum = xfer_getfnum();

          /* save the data to disk */
          xfer_out_disk(msg.fnum, (char *)sendbuf+sdispls[p]*sdtsize, sendcounts[p], sendtype);

          /* tell the master that we are storing the data on disk */
          xfer_out_scb(job->scb, &msg, sizeof(bdmsg_t), BDMPI_BYTE);
        }
      }
    }
  }

  /* prepare to go to sleep */
  S_SB_IFSET(BDMPI_SB_SAVEALL) {
    if (job->jdesc->nr < job->jdesc->ns)
      sb_saveall();
  }
  xfer_out_scb(job->scb, &sleeping, sizeof(int), BDMPI_BYTE);

  /* go to sleep until everybody has called the collective */
  BDMPL_SLEEP(job, gomsg);

  /*=====================================================================*/
  /* after waking up... */
  /*=====================================================================*/

  if (mype != root) {
    /* all but the root will send a SCATTERF request and get the data */
    msg.msgtype = BDMPI_MSGTYPE_SCATTERF;
    bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

    /* get information about from-memory data */
    xfer_in_scb(job->scb, &rmsg, sizeof(bdmsg_t), BDMPI_BYTE);

    if (bdmp_msize(rmsg.count, rmsg.datatype) > bdmp_msize(recvcount, recvtype))
      errexit("[%d]BDMPI_Scatterv: Amount of data to be received from %d is more than specified: %zu %zu\n",
          mype, root, bdmp_msize(rmsg.count, rmsg.datatype), bdmp_msize(recvcount, recvtype));

    /* get the data */
    if (rmsg.fnum == -1)
      xfer_in_scb(job->scb, recvbuf, rmsg.count, rmsg.datatype);
    else
      xfer_in_disk(rmsg.fnum, recvbuf, rmsg.count, rmsg.datatype, 1);
  }

  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Scatterv: exiting: comm: %p\n", comm));

  return BDMPI_SUCCESS;
}


/*************************************************************************/
/* Performs BDMPI_Scatterv() via a sequence of send/recv operations */
/*************************************************************************/
int bdmp_Scatterv_p2p(sjob_t *job,
          void *sendbuf, size_t *sendcounts, size_t *sdispls, BDMPI_Datatype sendtype,
          void *recvbuf, size_t recvcount, BDMPI_Datatype recvtype, int root,
          BDMPI_Comm comm)
{
  size_t k, sdtsize, rdtsize, size;
  int mype, npes, tag;

  S_IFSET(BDMPI_DBG_IPCS,
      bdprintf("BDMPI_Scatterv_p2p: entering: comm: %p [goMQlen: %d]\n",
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
  if (root < 0 || root >= comm->size) {
    fprintf(stderr, "The root rank is invalid.\n");
    return BDMPI_ERR_ROOT;
  }


  npes = comm->size;
  mype = comm->rank;

  tag = (++comm->copid)*BDMPL_COPID_MULT + BDMPL_SCATTER_TAG;

  sdtsize = bdmp_sizeof(sendtype);
  rdtsize = bdmp_sizeof(recvtype);

  /* the root sends data to everybody */
  if (mype == root) {
    for (k=0; k<npes; k++) {
      if (k != mype)
        bdmp_Send(job, (char *)sendbuf+sdispls[k]*sdtsize, sendcounts[k], sendtype,
            k, tag, comm);
    }

    BDASSERT(rdtsize*recvcount >= sdtsize*sendcounts[mype]);
    size = gk_min(sdtsize*sendcounts[mype], rdtsize*recvcount);
    memcpy((char *)recvbuf, (char *)sendbuf+sdispls[mype]*sdtsize, size);
  }

  /* sync to ensure collective semantics */
  bdmp_Barrier(job, comm);

  /* everybody receives data from the root */
  if (mype != root)
    bdmp_Recv(job, (char *)recvbuf, recvcount, recvtype, root, tag, comm,
        BDMPI_STATUS_IGNORE);


  return BDMPI_SUCCESS;
}
