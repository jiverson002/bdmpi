/*!
\file
\brief Implements the various recv operations.
\date Started 4/6/2013
\author George
*/


#include "bdmplib.h"


/*************************************************************************/
/* Performs a blocking probe operation */
/*************************************************************************/
int bdmp_Probe(sjob_t *job, int source, int tag, BDMPI_Comm comm,
          BDMPI_Status *status)
{
  int response, flag=0;
  bdmsg_t msg, rmsg, gomsg;

  S_IFSET(BDMPI_DBG_IPCS,
      bdprintf("BDMPI_Probe: Receiving from %d [goMQlen: %d]\n", source, bdmq_length(job->goMQ)));

  /* some error checking */
  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }

  /* save the state in case you go to sleep */
  S_SB_IFSET(BDMPI_SB_SAVEALL) {
    bdmp_Iprobe(job, source, tag, comm, &flag, BDMPI_STATUS_IGNORE);
    if (flag == 0 && job->jdesc->nr < job->jdesc->ns)
      sb_saveall();
  }

  msg.msgtype  = BDMPI_MSGTYPE_PROBE;
  msg.mcomm    = comm->mcomm;
  msg.myrank   = comm->rank;
  msg.tag      = tag;
  msg.source   = source;
  msg.dest     = comm->rank;

  for (;;) {
    /* notify the master that you want to receive a message */
    bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

    /* get the master's response  */
    bdmq_recv(job->c2sMQ, &response, sizeof(int));

    if (response == 1)
      break;  /* we got the go-ahead */

    /* prepare to go to sleep */
    S_SB_IFSET(BDMPI_SB_SAVEALL) {
      if (job->jdesc->nr < job->jdesc->ns)
        sb_saveall();
    }

    /* go to sleep... */
    BDMPL_SLEEP(job, gomsg);
  }

  /* get the missing message info from the master */
  xfer_in_scb(job->scb, &rmsg, sizeof(bdmsg_t), BDMPI_BYTE);

  if (status != BDMPI_STATUS_IGNORE) {
    status->BDMPI_SOURCE = rmsg.source;
    status->BDMPI_TAG    = rmsg.tag;
    status->BDMPI_ERROR  = BDMPI_SUCCESS;
    status->MPI_SOURCE   = status->BDMPI_SOURCE;
    status->MPI_TAG      = status->BDMPI_TAG;
    status->MPI_ERROR    = status->BDMPI_ERROR;
    status->count        = rmsg.count;
    status->datatype     = rmsg.datatype;
  }

  return BDMPI_SUCCESS;
}


/*************************************************************************/
/* Performs a non-blocking probe operation */
/*************************************************************************/
int bdmp_Iprobe(sjob_t *job, int source, int tag, BDMPI_Comm comm,
          int *flag, BDMPI_Status *status)
{
  bdmsg_t msg, rmsg;

  S_IFSET(BDMPI_DBG_IPCS,
      bdprintf("BDMPI_Iprobe: Probing from %d [goMQlen: %d]\n", source, bdmq_length(job->goMQ)));

  /* some error checking */
  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }


  memset(&msg, 0, sizeof(bdmsg_t));
  msg.msgtype  = BDMPI_MSGTYPE_IPROBE;
  msg.mcomm    = comm->mcomm;
  msg.myrank   = comm->rank;
  msg.tag      = tag;
  msg.source   = source;
  msg.dest     = comm->rank;


  /* notify the master that you want to receive a message */
  bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

  /* get the master's response  */
  bdmq_recv(job->c2sMQ, flag, sizeof(int));

  if (*flag) { /* get the status info */
    xfer_in_scb(job->scb, &rmsg, sizeof(bdmsg_t), BDMPI_BYTE);

    if (status != BDMPI_STATUS_IGNORE) {
      status->BDMPI_SOURCE = rmsg.source;
      status->BDMPI_TAG    = rmsg.tag;
      status->BDMPI_ERROR  = BDMPI_SUCCESS;
      status->MPI_SOURCE   = status->BDMPI_SOURCE;
      status->MPI_TAG      = status->BDMPI_TAG;
      status->MPI_ERROR    = status->BDMPI_ERROR;
      status->count        = rmsg.count;
      status->datatype     = rmsg.datatype;
    }
  }

  return BDMPI_SUCCESS;
}
