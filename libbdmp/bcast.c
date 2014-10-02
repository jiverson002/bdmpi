/*!
\file
\brief Implements the bcast operation.
\date Started 4/12/2013
\author George
*/

#include "bdmplib.h"


/*************************************************************************/
/* Performs BDMPI_Bcast() */
/*************************************************************************/
int bdmp_Bcast(sjob_t *job, void *buf, size_t count, BDMPI_Datatype datatype, 
          int root, BDMPI_Comm comm)
{
  int mype, response, sleeping=1;
  bdmsg_t msg, gomsg;

  if (bdmq_length(job->goMQ) != 0)
    bdprintf("BDMPI_Bcast: goMQ length is != 0.");

  S_IFSET(BDMPI_DBG_IPCS, 
      bdprintf("BDMPI_Bcast: entering: root: %d, comm: %p [goMQlen: %d]\n", 
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
  if (root < 0 || root >= comm->size) {
    fprintf(stderr, "The root rank is invalid.\n");
    return BDMPI_ERR_ROOT;
  }


  mype = comm->rank;

  /* notify the master that you entering a bcast */
  msg.msgtype  = BDMPI_MSGTYPE_BCASTI;
  msg.mcomm    = comm->mcomm;
  msg.myrank   = mype;
  msg.fnum     = -1;
  msg.count    = count;
  msg.datatype = datatype;
  msg.source   = root;
  msg.dest     = mype;

  bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

  /* get the copid from the master */
  bdmq_recv(job->c2sMQ, &msg.copid, sizeof(int));

  /* the root sends the data to the master */
  if (mype == root) {
    S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Bcast: Copying %zu elements.\n", count));
    xfer_out_scb(job->scb, buf, count, datatype);
  }


  /* prepare to go to sleep */
  /*if (mype != root)
    sb_discard(buf, bdmp_msize(count, datatype));*/
  /*if (job->jdesc->nr < job->jdesc->ns)
    sb_saveall();*/
  xfer_out_scb(job->scb, &sleeping, sizeof(int), BDMPI_BYTE);

  /* go to sleep until everybody has called the bcast */
  for (;;) {
    if (-1 == bdmq_recv(job->goMQ, &gomsg, sizeof(bdmsg_t)))
      bdprintf("Failed on trying to recv a go message: %s.\n", strerror(errno));
    if (BDMPI_MSGTYPE_PROCEED == gomsg.msgtype)
      break;
    slv_route(job, &gomsg);
  }

  /* all but the root will send a BCASTF request and get the data */
  if (mype != root) {
    /* notify the master that you want to receive the bcasted data */
    msg.msgtype = BDMPI_MSGTYPE_BCASTF;
    bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));
  
    /* copy the data from the master */
    xfer_in_scb(job->scb, buf, count, datatype);

    /* tell the master that you got the data */
    bdmq_send(job->c2mMQ, &response, sizeof(int));
  }

  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Bcast: exiting: comm: %p\n", comm));

  return BDMPI_SUCCESS;
}


