/*!
\file
\brief Implements the reduce operation.
\date Started 4/15/2013
\author George
*/


#include "bdmplib.h"


#define MERGE_APPLY_OP1(i, a, b, count, op) \
  do {\
    switch (op) {\
      case BDMPI_MIN:\
        for (i=0; i<count; i++) \
          (a)[i] = ((a)[i]<(b)[i] ? (a)[i] : (b)[i]);\
        break;\
      case BDMPI_MAX:\
        for (i=0; i<count; i++) \
          (a)[i] = ((a)[i]>(b)[i] ? (a)[i] : (b)[i]);\
        break;\
      case BDMPI_SUM:\
        for (i=0; i<count; i++) \
          (a)[i] += (b)[i];\
        break;\
      case BDMPI_PROD:\
        for (i=0; i<count; i++) \
          (a)[i] *= (b)[i];\
        break;\
      case BDMPI_LAND:\
        for (i=0; i<count; i++) \
          (a)[i] = (a)[i] && (b)[i];\
        break;\
      case BDMPI_BAND:\
        for (i=0; i<count; i++) \
          (a)[i] &= (b)[i];\
        break;\
      case BDMPI_LOR:\
        for (i=0; i<count; i++) \
          (a)[i] = (a)[i] || (b)[i];\
        break;\
      case BDMPI_BOR:\
        for (i=0; i<count; i++) \
          (a)[i] |= (b)[i];\
        break;\
      case BDMPI_BXOR:\
        for (i=0; i<count; i++) \
          (a)[i] ^= (b)[i];\
        break;\
      default:\
        errexit("BDMPI_Op: %d not implemented\n", op);\
    }\
  } while (0)

#define MERGE_APPLY_OP2(i, a, b, count, op) \
  do {\
    switch (op) {\
      case BDMPI_MIN:\
        for (i=0; i<count; i++) \
          (a)[i] = ((a)[i]<(b)[i] ? (a)[i] : (b)[i]);\
        break;\
      case BDMPI_MAX:\
        for (i=0; i<count; i++) \
          (a)[i] = ((a)[i]>(b)[i] ? (a)[i] : (b)[i]);\
        break;\
      case BDMPI_SUM:\
        for (i=0; i<count; i++) \
          (a)[i] += (b)[i];\
        break;\
      case BDMPI_PROD:\
        for (i=0; i<count; i++) \
          (a)[i] *= (b)[i];\
        break;\
      default:\
        errexit("BDMPI_Op: %d not implemented\n", op);\
    }\
  } while (0)


/*************************************************************************/
/*! Performs BDMPI_Merge() */
/*************************************************************************/
int bdmp_Merge(sjob_t *job, void *sendbuf, int *sendids, int sendcount,
          void *recvbuf, int *recvids, size_t *r_recvcount,
          BDMPI_Datatype datatype, BDMPI_Op op, int root, BDMPI_Comm comm)
{
  size_t i, len, chunk, size, dtsize;
  int mype, sleeping=1;
  bdmsg_t msg, gomsg;


  S_IFSET(BDMPI_DBG_IPCS,
      bdprintf("BDMPI_Merge: entering: root: %d, comm: %o [goMQlen: %d]\n",
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


  mype = comm->rank;

  /* notify the master that you entering a reduce */
  msg.msgtype  = BDMPI_MSGTYPE_MERGEI;
  msg.mcomm    = comm->mcomm;
  msg.myrank   = mype;
  msg.fnum     = -1;
  msg.count    = sendcount;
  msg.datatype = datatype;
  msg.dest     = root;
  msg.source   = mype;
  msg.op       = op;

  bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

  /* the root gets the copid message */
  if (mype == root)
    bdmq_recv(job->c2sMQ, &msg.copid, sizeof(int));


  /* everybody sends the data to the master */
  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Merge: Copying %zu elements.\n", sendcount));
  xfer_out_scb(job->scb, sendbuf, sendcount, datatype);
  xfer_out_scb(job->scb, sendids, sendcount, BDMPI_INT);

  /* prepare to go to sleep */
  S_SB_IFSET(BDMPI_SB_SAVEALL) {
    if (job->jdesc->nr < job->jdesc->ns)
      sbma_mevictall();
  }
  xfer_out_scb(job->scb, &sleeping, sizeof(int), BDMPI_BYTE);

  /* go to sleep until everybody has called the reduce */
  BDMPL_SLEEP(job, gomsg);

  /* the root sends a REDUCE_RECV request and get the data */
  if (mype == root) {
    /* notify the master that you want to receive the reduced data */
    msg.msgtype = BDMPI_MSGTYPE_MERGEF;
    msg.count   = *r_recvcount;
    bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

    /* copy the data from the master */
    xfer_in_scb(job->scb, r_recvcount, 1, BDMPI_SIZE_T);
    xfer_in_scb(job->scb, recvbuf, *r_recvcount, datatype);
    xfer_in_scb(job->scb, recvids, *r_recvcount, BDMPI_INT);
  }

  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Merge: exiting: comm: %p\n", comm));

  return BDMPI_SUCCESS;
}
