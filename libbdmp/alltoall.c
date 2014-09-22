/*!
\file
\brief Implements the alltoall operation.
\date Started 5/3/2013
\author George
*/

#include "bdmplib.h"


/*************************************************************************/
/* Performs BDMPI_Alltoallv() directly */
/*************************************************************************/
int bdmp_Alltoallv_node(sjob_t *job, 
          void *sendbuf, size_t *sendcounts, size_t *sdispls, BDMPI_Datatype sendtype, 
          void *recvbuf, size_t *recvcounts, size_t *rdispls, BDMPI_Datatype recvtype, 
          BDMPI_Comm comm)
{
  int npes, mype, sleeping=1;
  int i, p, response;
  size_t sdtsize, rdtsize;
  bdmsg_t msg, rmsg;
  ssize_t myfnum=-1;

  S_IFSET(BDMPI_DBG_IPCS, 
      bdprintf("BDMPI_Alltoallv_node: entering: comm: %p [goMQlen: %d]\n", 
          comm, bdmq_length(job->goMQ)));

  /* error checking */
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
  msg.msgtype  = BDMPI_MSGTYPE_ALLTOALLI;
  msg.mcomm    = comm->mcomm;
  msg.myrank   = mype;
  msg.datatype = sendtype;

  bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

  /* get the copid from the master */
  bdmq_recv(job->c2sMQ, &msg.copid, sizeof(int));

  /* start sending data to the master or writing them to the disk */
  for (p=0; p<npes; p++) {
    if (p == mype) {
      if (sendcounts[p]*sdtsize > recvcounts[p]*rdtsize)
        errexit("BDMPI_Alltoallv: Amount of data sent is greater than available buffer space.\n");

      if (sendcounts[p]*sdtsize > job->smallmsg) { /* store the data on disk */
        myfnum = xfer_getfnum();
        xfer_out_disk(myfnum, (char *)sendbuf+sdispls[p]*sdtsize, sendcounts[p], sendtype);
      }
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

        /* store the data */
        xfer_out_disk(msg.fnum, (char *)sendbuf+sdispls[p]*sdtsize, sendcounts[p], sendtype);

        /* tell the master about it */
        xfer_out_scb(job->scb, &msg, sizeof(bdmsg_t), BDMPI_BYTE);
      }
    }
  }

  /* prepare to go to sleep */
  if (job->jdesc->nr < job->jdesc->ns)
    sb_saveall();
  xfer_out_scb(job->scb, &sleeping, sizeof(int), BDMPI_BYTE);


  /* go to sleep until everybody has called the collective */
  bdmq_recv(job->goMQ, &response, sizeof(int));


  /*=====================================================================*/
  /* after waking up... */
  /*=====================================================================*/

  /* all will send a ALLGATHERF request and get the data */
  msg.msgtype = BDMPI_MSGTYPE_ALLTOALLF;
  bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));
  
  /* get the data */
  for (i=0; i<npes-1; i++) {
    //printf("%d slv-recv1 %2d.%2d\n", (int)time(0), mype, (int)i); 

    /* get msg info from master */
    xfer_in_scb(job->scb, &rmsg, sizeof(bdmsg_t), BDMPI_BYTE);
    p = rmsg.myrank;

    if (bdmp_msize(rmsg.count, rmsg.datatype) > bdmp_msize(recvcounts[p], recvtype))
      errexit("[%d]BDMPI_Alltoallv: Amount of data to be received from %d is more than specified: %zu %zu\n", 
          mype, p, bdmp_msize(rmsg.count, rmsg.datatype), bdmp_msize(recvcounts[p], recvtype));

    //printf("%d slv-recv2 %2d.%2d\n", (int)time(0), mype, (int)i); 
    if (rmsg.fnum == -1)
      xfer_in_scb(job->scb, (char *)recvbuf+rdispls[p]*rdtsize, rmsg.count, rmsg.datatype);
    else
      xfer_in_disk(rmsg.fnum, (char *)recvbuf+rdispls[p]*rdtsize, rmsg.count, rmsg.datatype, 1);
    //printf("%d slv-recv3 %2d.%2d\n", (int)time(0), mype, (int)i); 
  }

  /* copy the local data */
  if (myfnum == -1) 
    memcpy((char *)recvbuf+rdispls[mype]*rdtsize, (char *)sendbuf+sdispls[mype]*sdtsize, 
        sendcounts[mype]*sdtsize);
  else
    xfer_in_disk(myfnum, (char *)recvbuf+rdispls[mype]*rdtsize, sendcounts[mype], sendtype, 1);


  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Alltoallv: exiting: comm: %p\n", comm));

  return BDMPI_SUCCESS;
}


/*************************************************************************/
/* Performs BDMPI_Alltoallv() via a sequence of send/recv operations */
/*************************************************************************/
int bdmp_Alltoallv_p2p0(sjob_t *job, 
          void *sendbuf, size_t *sendcounts, size_t *sdispls, BDMPI_Datatype sendtype, 
          void *recvbuf, size_t *recvcounts, size_t *rdispls, BDMPI_Datatype recvtype, 
          BDMPI_Comm comm)
{
  size_t sdtsize, rdtsize, size;
  int i, p, mype, npes;
  BDMPI_Status status;

  S_IFSET(BDMPI_DBG_IPCS, 
      bdprintf("BDMPI_Alltoallv_p2p: entering: comm: %p [goMQlen: %d]\n", 
          comm, bdmq_length(job->goMQ)));

  /* error checking */
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

  /* everybody sends data to everybody */
  sdtsize = bdmp_sizeof(sendtype);
  for (i=1; i<npes; i++) {
    p = (mype+i)%npes;
    bdmp_Send(job, (char *)sendbuf+sdispls[p]*sdtsize, sendcounts[p], 
          sendtype, p, BDMPL_ALLTOALL_TAG, comm);
  }

  //bdmp_Barrier(job, comm);

#ifdef XXX
  /* everybody receives data from everybody */
  rdtsize = bdmp_sizeof(recvtype);
  for (i=1; i<npes; i++) {
    p = (mype+i)%npes;
    bdmp_Recv(job, (char *)recvbuf+rdispls[p]*rdtsize, recvcounts[p], 
          recvtype, p, BDMPL_ALLTOALL_TAG, comm, BDMPI_STATUS_IGNORE);
  }
#endif

  /* everybody receives data from everybody */
  rdtsize = bdmp_sizeof(recvtype);
  for (i=1; i<npes; i++) {
    bdmp_Probe(job, BDMPI_ANY_SOURCE, BDMPL_ALLTOALL_TAG, comm, &status);

    p = status.BDMPI_SOURCE;
    bdmp_Recv(job, (char *)recvbuf+rdispls[p]*rdtsize, recvcounts[p], 
          recvtype, p, BDMPL_ALLTOALL_TAG, comm, BDMPI_STATUS_IGNORE);
  }

  /* sync to ensure collective semantics */
  /* this is not required, as the success of the above receives ensures that */
  /* TODO: it is here due to some weird race condition */
  //bdmp_Barrier(job, comm);

  BDASSERT(rdtsize*recvcounts[mype] >= sdtsize*sendcounts[mype]);
  size = gk_min(sdtsize*sendcounts[mype], rdtsize*recvcounts[mype]);
  memcpy((char *)recvbuf+rdispls[mype]*rdtsize, (char *)sendbuf+sdispls[mype]*sdtsize, size);

  return BDMPI_SUCCESS;
}


/*************************************************************************/
/* Performs BDMPI_Alltoallv() via a sequence of send/recv operations */
/*************************************************************************/
int bdmp_Alltoallv_p2p(sjob_t *job, 
          void *sendbuf, size_t *sendcounts, size_t *sdispls, BDMPI_Datatype sendtype, 
          void *recvbuf, size_t *recvcounts, size_t *rdispls, BDMPI_Datatype recvtype, 
          BDMPI_Comm comm)
{
  size_t sdtsize, rdtsize, size;
  int i, p, mype, npes, tag, response;
  bdmsg_t msg, rmsg, *rmsgs;
  BDMPI_Status status;

  S_IFSET(BDMPI_DBG_IPCS, 
      bdprintf("BDMPI_Alltoallv_p2p: entering: comm: %p [goMQlen: %d]\n", 
          comm, bdmq_length(job->goMQ)));

  /* error checking */
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

  tag = (++comm->copid)*BDMPL_COPID_MULT + BDMPL_ALLTOALL_TAG;

  rmsgs = (bdmsg_t *)gk_malloc(sizeof(bdmsg_t)*npes, "BDMPI_Alltoallv_p2p: rmsgs");
  memset(rmsgs, 0, sizeof(bdmsg_t)*npes);

  sdtsize = bdmp_sizeof(sendtype);
  rdtsize = bdmp_sizeof(recvtype);


  /* send your data to everybody else */
  for (i=1; i<npes; i++) {
    p = (mype+i)%npes;
    bdmp_Send(job, (char *)sendbuf+sdispls[p]*sdtsize, sendcounts[p], sendtype, p, 
        tag, comm);
  }

  /* deal with your data */
  rmsgs[mype].fnum     = -1;
  rmsgs[mype].source   = mype;
  rmsgs[mype].count    = sendcounts[mype];
  rmsgs[mype].datatype = sendtype;
  if (sendcounts[mype]*sdtsize > job->smallmsg) { 
    rmsgs[mype].fnum = xfer_getfnum();
    xfer_out_disk(rmsgs[mype].fnum, (char *)sendbuf+sdispls[mype]*sdtsize, 
        sendcounts[mype], sendtype);
  }


  /* sbdiscard the incoming buffers */
  for (p=0; p<npes; p++) 
    sb_discard((char *)recvbuf+rdispls[p]*rdtsize, recvcounts[p]*rdtsize);

  /* save your address space before blocking */
  if (job->jdesc->nr < job->jdesc->ns)
    sb_saveall();

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
      if (job->jdesc->nr < job->jdesc->ns)
        sb_saveall();
      bdmq_recv(job->goMQ, &response, sizeof(int));
    }

    /* get the missing message info from the master */
    xfer_in_scb(job->scb, &rmsg, sizeof(bdmsg_t), BDMPI_BYTE);
    p = rmsg.source;
    rmsgs[p] = rmsg;

    /* check that we have enough buffer space */
    if (bdmp_msize(rmsgs[p].count, rmsgs[p].datatype) > recvcounts[p]*rdtsize)
      errexit("BDMPI_Alltoallv_p2p: Received message is greater than provided buffer: recv: %zu; buffer: %zu\n",
          bdmp_msize(rmsgs[p].count, rmsgs[p].datatype), recvcounts[p]*rdtsize);

    /* copy the in-memory data for now and defer disk for later */
    if (rmsgs[p].fnum == -1)
      xfer_in_scb(job->scb, (char *)recvbuf+rdispls[p]*rdtsize, rmsgs[p].count, rmsgs[p].datatype);
  }

  /* once everything has been received, copy the data from the disk in memory */
  for (p=0; p<npes; p++) {
    if (rmsgs[p].fnum != -1)
      xfer_in_disk(rmsgs[p].fnum, (char *)recvbuf+rdispls[p]*rdtsize, rmsgs[p].count, 
          rmsgs[p].datatype, 1);
  }

  /* deal with mype's potential in-memory copy */
  if (rmsgs[mype].fnum == -1) {
    BDASSERT(rdtsize*recvcounts[mype] >= sdtsize*sendcounts[mype]);
    size = gk_min(sdtsize*sendcounts[mype], rdtsize*recvcounts[mype]);
    memcpy((char *)recvbuf+rdispls[mype]*rdtsize, (char *)sendbuf+sdispls[mype]*sdtsize, size);
  }

  gk_free((void **)&rmsgs, LTERM);

  return BDMPI_SUCCESS;
}

