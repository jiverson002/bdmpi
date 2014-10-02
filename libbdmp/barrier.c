/*!
\file
\brief Implements the barrier operation.
\date Started 4/12/2013
\author George
*/

#include "bdmplib.h"


/*************************************************************************/
/* Performs BDMPI_Barrier() */
/*************************************************************************/
int bdmp_Barrier(sjob_t *job, BDMPI_Comm comm)
{
  int response;
  bdmsg_t msg, gomsg;

  S_IFSET(BDMPI_DBG_IPCS, 
      bdprintf("BDMPI_Barrier: entering: comm: %p [goMQlen: %d]\n", comm, 
          bdmq_length(job->goMQ)));

  /* some error checking */
  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }

  msg.msgtype = BDMPI_MSGTYPE_BARRIER;
  msg.mcomm   = comm->mcomm;
  msg.myrank  = comm->rank;

  /* save the address space before blocking */
  /*if (job->jdesc->nr < job->jdesc->ns)
    sb_saveall();*/

  /* notify the master that you entering a barrier */
  bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

  /* go to sleep... */
  for (;;) {
    if (-1 == bdmq_recv(job->goMQ, &gomsg, sizeof(bdmsg_t)))
      bdprintf("Failed on trying to recv a go message: %s.\n", strerror(errno));
    if (BDMPI_MSGTYPE_PROCEED == gomsg.msgtype)
      break;
  }

  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Barrier: exiting: comm: %p\n", comm));

  return BDMPI_SUCCESS;
}

