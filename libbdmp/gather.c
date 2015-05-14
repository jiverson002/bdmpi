/*!
\file
\brief Implements the gather operation.
\date Started 5/3/2013
\author George
*/


#include "bdmplib.h"


/*************************************************************************/
/*! Performs BDMPI_Gatherv() when the communicator involves a single node */
/*************************************************************************/
int bdmp_Gatherv_node(sjob_t *job,
          void *sendbuf, size_t sendcount, BDMPI_Datatype sendtype, void *recvbuf,
          size_t *recvcounts, size_t *rdispls, BDMPI_Datatype recvtype, int root,
          BDMPI_Comm comm)
{
  size_t sdtsize, rdtsize;
  int npes, mype, i, p, sleeping=1;
  bdmsg_t msg, rmsg, gomsg;
  ssize_t myfnum=-1;

  S_IFSET(BDMPI_DBG_IPCS,
      bdprintf("BDMPI_Gatherv_node: entering: comm: %p [goMQlen: %d]\n",
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


  npes = comm->size;
  mype = comm->rank;

  sdtsize = bdmp_sizeof(sendtype);
  rdtsize = bdmp_sizeof(recvtype);

  /* notify the master that you entering an allgather */
  msg.msgtype  = BDMPI_MSGTYPE_GATHERI;
  msg.mcomm    = comm->mcomm;
  msg.myrank   = mype;
  msg.count    = sendcount;
  msg.dest     = root;
  msg.datatype = sendtype;

  bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

  /* get the copid from the master */
  bdmq_recv(job->c2sMQ, &msg.copid, sizeof(int));

  /* send the data to the master or write them to the disk */
  if (mype == root) {
    if (sendcount*sdtsize > recvcounts[mype]*rdtsize)
      errexit("BDMPI_Gatherv: Receive buffer for %d is not sufficient.\n", mype);

    if (sendcount*sdtsize > job->smallmsg) { /* store the data on disk */
      myfnum = xfer_getfnum();
      xfer_out_disk(myfnum, sendbuf, sendcount, sendtype);
    }
  }
  else {
    if (sendcount*sdtsize <= job->smallmsg) { /* to memory */
      /* tell the master that we are sending it some data */
      msg.fnum = -1;
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
  if (mype == root) {
    /* the root will send a GATHERF request and get the data */
    msg.msgtype = BDMPI_MSGTYPE_GATHERF;
    bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

    /* get the data */
    for (i=0; i<npes-1; i++) {
      xfer_in_scb(job->scb, &rmsg, sizeof(bdmsg_t), BDMPI_BYTE);
      p = rmsg.myrank;

      if (bdmp_msize(rmsg.count, rmsg.datatype) > bdmp_msize(recvcounts[p], recvtype))
        errexit("[%d]BDMPI_Gatherv: Amount of data to be received from %d is more than specified: %zu %zu\n",
            mype, p, bdmp_msize(rmsg.count, rmsg.datatype), bdmp_msize(recvcounts[p], recvtype));

      if (rmsg.fnum == -1)
        xfer_in_scb(job->scb, (char *)recvbuf+rdispls[p]*rdtsize, rmsg.count, rmsg.datatype);
      else
        xfer_in_disk(rmsg.fnum, (char *)recvbuf+rdispls[p]*rdtsize, rmsg.count, rmsg.datatype, 1);
    }

    if (myfnum == -1)
      memcpy((char *)recvbuf+rdispls[mype]*rdtsize, (char *)sendbuf, sendcount*sdtsize);
    else
      xfer_in_disk(myfnum, (char *)recvbuf+rdispls[mype]*rdtsize, sendcount, sendtype, 1);
  }


  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Gatherv: exiting: comm: %p\n", comm));

  return BDMPI_SUCCESS;
}


/*************************************************************************/
/* Performs BDMPI_Gatherv() via a sequence of send/recv operations */
/*************************************************************************/
int bdmp_Gatherv_p2p0(sjob_t *job,
          void *sendbuf, size_t sendcount, BDMPI_Datatype sendtype, void *recvbuf,
          size_t *recvcounts, size_t *rdispls, BDMPI_Datatype recvtype, int root,
          BDMPI_Comm comm)
{
  size_t k, sdtsize, rdtsize, size;
  int mype, npes;

  S_IFSET(BDMPI_DBG_IPCS,
      bdprintf("BDMPI_Gatherv_p2p: entering: comm: %p [goMQlen: %d]\n",
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

  npes = comm->size;
  mype = comm->rank;

  sdtsize = bdmp_sizeof(sendtype);
  rdtsize = bdmp_sizeof(recvtype);

  /* the non-root sends data to the root */
  if (mype != root)
    bdmp_Send(job, (char *)sendbuf, sendcount, sendtype, root, BDMPL_GATHER_TAG, comm);

  /* the root receives data from everybody */
  if (mype == root) {
    for (k=0; k<npes; k++) {
      if (k != mype)
        bdmp_Recv(job, (char *)recvbuf+rdispls[k]*rdtsize, recvcounts[k],
            recvtype, k, BDMPL_GATHER_TAG, comm, BDMPI_STATUS_IGNORE);
    }

    BDASSERT(rdtsize*recvcounts[mype] >= sdtsize*sendcount);
    size = gk_min(sdtsize*sendcount, rdtsize*recvcounts[mype]);
    memcpy((char *)recvbuf+rdispls[mype]*rdtsize, (char *)sendbuf, size);
  }

  /* sync to ensure collective semantics */
  bdmp_Barrier(job, comm);

  return BDMPI_SUCCESS;
}


/*************************************************************************/
/* Performs BDMPI_Gatherv() via a sequence of send/recv operations */
/*************************************************************************/
int bdmp_Gatherv_p2p(sjob_t *job,
          void *sendbuf, size_t sendcount, BDMPI_Datatype sendtype, void *recvbuf,
          size_t *recvcounts, size_t *rdispls, BDMPI_Datatype recvtype, int root,
          BDMPI_Comm comm)
{
  size_t i, p, sdtsize, rdtsize, size;
  int mype, npes, tag, response;
  bdmsg_t msg, rmsg, gomsg, *rmsgs;

  S_IFSET(BDMPI_DBG_IPCS,
      bdprintf("BDMPI_Gatherv_p2p: entering: comm: %p [goMQlen: %d]\n",
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

  npes = comm->size;
  mype = comm->rank;

  tag = (++comm->copid)*BDMPL_COPID_MULT + BDMPL_GATHER_TAG;

  sdtsize = bdmp_sizeof(sendtype);
  rdtsize = bdmp_sizeof(recvtype);

  /* the non-root sends data to the root */
  if (mype != root)
    bdmp_Send(job, (char *)sendbuf, sendcount, sendtype, root, tag, comm);

  /* save the data in case you go to sleep */
  S_SB_IFSET(BDMPI_SB_SAVEALL) {
    if (job->jdesc->nr < job->jdesc->ns)
      sb_saveall();
  }

  /* sync to ensure collective semantics */
  bdmp_Barrier(job, comm);

  /* the root receives data from everybody */
  if (mype == root) {
    rmsgs = (bdmsg_t *)bd_malloc(sizeof(bdmsg_t)*npes, "BDMPI_Gatherv_p2p: rmsgs");
    memset(rmsgs, 0, sizeof(bdmsg_t)*npes);

    /* deal with root's data */
    rmsgs[mype].fnum     = -1;
    rmsgs[mype].source   = mype;
    rmsgs[mype].count    = sendcount;
    rmsgs[mype].datatype = sendtype;
    if (sendcount*sdtsize > job->smallmsg) {
      rmsgs[mype].fnum = xfer_getfnum();
      xfer_out_disk(rmsgs[mype].fnum, (char *)sendbuf, sendcount, sendtype);
    }

    /* sbdiscard the incoming buffers */
    S_SB_IFSET(BDMPI_SB_DISCARD) {
      for (p=0; p<npes; p++)
        sb_discard((char *)recvbuf+rdispls[p]*rdtsize,
          bdmp_msize(recvcounts[p], recvtype));
    }

    /* receive data from everybody else */
    msg.msgtype  = BDMPI_MSGTYPE_RECV;
    msg.mcomm    = comm->mcomm;
    msg.myrank   = mype;
    msg.dest     = mype;
    msg.tag      = tag;
    msg.source   = BDMPI_ANY_SOURCE;
    msg.datatype = recvtype;

    msg.count = recvcounts[0];
    for (p=1; p<npes; p++)
      msg.count = (msg.count < recvcounts[p] ? recvcounts[p] : msg.count);

    for (i=0; i<npes-1; i++) {
      for (;;) {
        /* notify the master that you want to receive a message */
        bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

        /* get the master's response  */
        bdmq_recv(job->c2sMQ, &response, sizeof(int));

        if (response == 1)
          break;

        /* go to sleep... */
        S_SB_IFSET(BDMPI_SB_SAVEALL) {
          if (job->jdesc->nr < job->jdesc->ns)
            sb_saveall();
        }
        BDMPL_SLEEP(job, gomsg);
      }

      /* get the missing message info from the master */
      xfer_in_scb(job->scb, &rmsg, sizeof(bdmsg_t), BDMPI_BYTE);
      p = rmsg.source;
      rmsgs[p] = rmsg;

      /* check that we have enough buffer space */
      if (bdmp_msize(rmsgs[p].count, rmsgs[p].datatype) > recvcounts[p]*rdtsize)
        errexit("BDMPI_Gatherv_p2p: Received message is greater than provided buffer: recv: %zu; buffer: %zu\n",
            bdmp_msize(rmsgs[p].count, rmsgs[p].datatype), recvcounts[p]*rdtsize);

      /* copy the in-memory data for now and defer disk for later */
      if (rmsgs[p].fnum == -1)
        xfer_in_scb(job->scb, (char *)recvbuf+rdispls[p]*rdtsize, rmsgs[p].count,
            rmsgs[p].datatype);
    }


    /* once everything has been received, copy the data from the disk in memory */
    for (p=0; p<npes; p++) {
      if (rmsgs[p].fnum != -1)
        xfer_in_disk(rmsgs[p].fnum, (char *)recvbuf+rdispls[p]*rdtsize,
            rmsgs[p].count, rmsgs[p].datatype, 1);
    }

    /* deal with root's potential in-memory copy */
    if (rmsgs[mype].fnum == -1) {
      BDASSERT(rdtsize*recvcounts[mype] >= sdtsize*sendcount);
      size = gk_min(sdtsize*sendcount, rdtsize*recvcounts[mype]);
      memcpy((char *)recvbuf+rdispls[mype]*rdtsize, (char *)sendbuf, size);
    }

    bd_free((void **)&rmsgs, LTERM);
  }


  return BDMPI_SUCCESS;
}
