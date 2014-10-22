/*!
\file
\brief Implements the reduce operation.
\date Started 4/15/2013
\author George
*/


#include "bdmplib.h"


/*************************************************************************/
/*! Performs BDMPI_Reduce() */
/*************************************************************************/
int bdmp_Reduce(sjob_t *job, void *sendbuf, void *recvbuf, size_t count,
          BDMPI_Datatype datatype, BDMPI_Op op, int root, BDMPI_Comm comm)
{
  int mype, response, sleeping=1;
  bdmsg_t msg, gomsg;

  S_IFSET(BDMPI_DBG_IPCS,
      bdprintf("BDMPI_Reduce: entering: root: %d, comm: %p [goMQlen: %d]\n",
        root, comm, bdmq_length(job->goMQ)));

  /* some error checking */
  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }
  if (!datatype_isvalid(datatype)) {
    fprintf(stderr, "The datatype is invalid.\n");
    return BDMPI_ERR_TYPE;
  }
  if (!op_isvalid(op)) {
    fprintf(stderr, "The op is invalid.\n");
    return BDMPI_ERR_OP;
  }
  if (root < 0 || root >= comm->size) {
    fprintf(stderr, "The root rank is invalid.\n");
    return BDMPI_ERR_ROOT;
  }


  mype = comm->rank;

  /* notify the master that you entering a reduce */
  msg.msgtype  = BDMPI_MSGTYPE_REDUCEI;
  msg.mcomm    = comm->mcomm;
  msg.myrank   = mype;
  msg.fnum     = -1;
  msg.count    = count;
  msg.datatype = datatype;
  msg.dest     = root;
  msg.source   = mype;
  msg.op       = op;

  bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

  /* the root gets the copid message */
  if (mype == root)
    bdmq_recv(job->c2sMQ, &msg.copid, sizeof(int));


  /* everybody sends the data to the master */
  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Reduce: Copying %zu elements.\n", count));
  xfer_out_scb(job->scb, sendbuf, count, datatype);

  /* prepare to go to sleep */
  if (mype == root)
    sb_discard(recvbuf, bdmp_msize(count, datatype));
  /*if (job->jdesc->nr < job->jdesc->ns)
    sb_saveall();*/
  xfer_out_scb(job->scb, &sleeping, sizeof(int), BDMPI_BYTE);


  /* go to sleep until everybody has called the reduce */
  BDMPI_SLEEP(job, gomsg);

  /* the root sends a REDUCE_RECV request and get the data */
  if (mype == root) {
    /* notify the master that you want to receive the reduced data */
    msg.msgtype = BDMPI_MSGTYPE_REDUCEF;
    bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

    /* copy the data from the master */
    xfer_in_scb(job->scb, recvbuf, count, datatype);
  }

  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Reduce: exiting: comm: %p\n", comm));

  return BDMPI_SUCCESS;
}


/*************************************************************************/
/*! Performs BDMPI_Reduce() */
/*************************************************************************/
int bdmp_Reduce_init(sjob_t *job, void *sendbuf, void *recvbuf, size_t count,
          BDMPI_Datatype datatype, BDMPI_Op op, int root, BDMPI_Comm comm,
          BDMPI_Request *r_request)
{
  int mype, response;
  bdmsg_t msg;
  bdrequest_t *request;

  S_IFSET(BDMPI_DBG_IPCS,
      bdprintf("BDMPI_Reduce_init: entering: root: %d, comm: %p [goMQlen: %d]\n",
        root, comm, bdmq_length(job->goMQ)));

  *r_request = BDMPI_REQUEST_NULL;

  /* some error checking */
  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }
  if (!datatype_isvalid(datatype)) {
    fprintf(stderr, "The datatype is invalid.\n");
    return BDMPI_ERR_TYPE;
  }
  if (!op_isvalid(op)) {
    fprintf(stderr, "The op is invalid.\n");
    return BDMPI_ERR_OP;
  }
  if (root < 0 || root >= comm->size) {
    fprintf(stderr, "The root rank is invalid.\n");
    return BDMPI_ERR_ROOT;
  }

  mype = comm->rank;

  /* notify the master that you entering a reduce */
  msg.msgtype  = BDMPI_MSGTYPE_REDUCEI;
  msg.mcomm    = comm->mcomm;
  msg.myrank   = mype;
  msg.fnum     = -1;
  msg.count    = count;
  msg.datatype = datatype;
  msg.dest     = root;
  msg.source   = mype;
  msg.op       = op;

  bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

  /* the root gets the copid message */
  request = *r_request = (bdrequest_t *)gk_malloc(sizeof(bdrequest_t), "request");
  memset(request, 0, sizeof(bdrequest_t));
  request->type = BDMPI_REQUEST_REDUCEI;
  request->copid = -1;
  if (mype == root) {
    bdmq_recv(job->c2sMQ, &msg.copid, sizeof(int));
    request->copid = msg.copid;
  }

  /* everybody sends the data to the master */
  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Reduce: Copying %zu elements.\n", count));
  xfer_out_scb(job->scb, sendbuf, count, datatype);

  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Reduce_init: exiting: comm: %p\n", comm));

  return BDMPI_SUCCESS;
}


/*************************************************************************/
/*! Performs BDMPI_Reduce() */
/*************************************************************************/
int bdmp_Reduce_fine(sjob_t *job, void *recvbuf, size_t count,
          BDMPI_Datatype datatype, BDMPI_Op op, int root, BDMPI_Comm comm,
          BDMPI_Request *r_request)
{
  int mype, response;
  bdmsg_t msg, gomsg;
  bdrequest_t *request = *r_request;

  S_IFSET(BDMPI_DBG_IPCS,
      bdprintf("BDMPI_Reduce_fine: entering: root: %d, comm: %p [goMQlen: %d]\n",
        root, comm, bdmq_length(job->goMQ)));

  /* some error checking */
  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }
  if (!datatype_isvalid(datatype)) {
    fprintf(stderr, "The datatype is invalid.\n");
    return BDMPI_ERR_TYPE;
  }
  if (!op_isvalid(op)) {
    fprintf(stderr, "The op is invalid.\n");
    return BDMPI_ERR_OP;
  }
  if (root < 0 || root >= comm->size) {
    fprintf(stderr, "The root rank is invalid.\n");
    return BDMPI_ERR_ROOT;
  }
  if (request->type != BDMPI_REQUEST_REDUCEI) {
    fprintf(stderr, "Incorrect request type.\n");
    return BDMPI_ERR_ARG;
  }


  mype = comm->rank;

  /* go to sleep until everybody has called the reduce */
  /*if (job->jdesc->nr < job->jdesc->ns)
    sb_saveall();*/
  BDMPI_SLEEP(job, gomsg);

  /* the root sends a REDUCE_RECV request and get the data */
  if (mype == root) {
    msg.msgtype  = BDMPI_MSGTYPE_REDUCEF;
    msg.mcomm    = comm->mcomm;
    msg.myrank   = mype;
    msg.fnum     = -1;
    msg.count    = count;
    msg.datatype = datatype;
    msg.dest     = root;
    msg.source   = mype;
    msg.op       = op;
    msg.copid    = request->copid;

    /* notify the master that you want to receive the reduced data */
    bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

    /* copy the data from the master */
    xfer_in_scb(job->scb, recvbuf, count, datatype);
  }

  gk_free((void **)r_request, LTERM);
  *r_request = BDMPI_REQUEST_NULL;

  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Reduce_fine: exiting: comm: %p\n", comm));

  return BDMPI_SUCCESS;
}


/*************************************************************************/
/*! Performs BDMPI_Allreduce() */
/*************************************************************************/
int bdmp_Allreduce(sjob_t *job, void *sendbuf, void *recvbuf, size_t count,
          BDMPI_Datatype datatype, BDMPI_Op op, BDMPI_Comm comm)
{
  int mype, response, sleeping=1;
  bdmsg_t msg, gomsg;

  S_IFSET(BDMPI_DBG_IPCS,
      bdprintf("BDMPI_Allreduce: entering: comm: %p [goMQlen: %d]\n",
        comm, bdmq_length(job->goMQ)));

  /* some error checking */
  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }
  if (!datatype_isvalid(datatype)) {
    fprintf(stderr, "The datatype is invalid.\n");
    return BDMPI_ERR_TYPE;
  }
  if (!op_isvalid(op)) {
    fprintf(stderr, "The op is invalid.\n");
    return BDMPI_ERR_OP;
  }


  mype = comm->rank;

  /* notify the master that you entering a reduce */
  msg.msgtype  = BDMPI_MSGTYPE_ALLREDUCEI;
  msg.mcomm    = comm->mcomm;
  msg.myrank   = mype;
  msg.count    = count;
  msg.datatype = datatype;
  msg.op       = op;

  bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

  /* get the copid from the master */
  bdmq_recv(job->c2sMQ, &msg.copid, sizeof(int));

  /* everybody sends the data to the master */
  xfer_out_scb(job->scb, sendbuf, count, datatype);

  /* prepare to go to sleep */
  sb_discard(recvbuf, count*bdmp_sizeof(datatype));
  /*if (job->jdesc->nr < job->jdesc->ns)
    sb_saveall();*/
  xfer_out_scb(job->scb, &sleeping, sizeof(int), BDMPI_BYTE);

  /* go to sleep until everybody has called the reduce */
  BDMPI_SLEEP(job, gomsg);

  /* everybody sends a ALLREDUCE_RECV request and get the data */
  msg.msgtype = BDMPI_MSGTYPE_ALLREDUCEF;

  /* notify the master that you want to receive the reduced data */
  bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

  /* copy the data from the master */
  xfer_in_scb(job->scb, recvbuf, count, datatype);

  /* tell the master that we got the data */
  bdmq_send(job->c2mMQ, &response, sizeof(int));

  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Allreduce: exiting: comm: %p\n", comm));

  return BDMPI_SUCCESS;
}
