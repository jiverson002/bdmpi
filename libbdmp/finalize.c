/*!
\file
\brief Implements the BDMPI_Finalize() function.
\date Started 4/5/2013
\author George
*/


#include "bdmplib.h"


/*************************************************************************/
/* Finalizes the BDMP library. */
/*************************************************************************/
int bdmp_Finalize(sjob_t *job)
{
  int i, response;
  bdmsg_t donemsg, gomsg;

  S_IFSET(BDMPI_DBG_IPCS, bdprintf("iBDMPI_Finalize: entering [goMQlen: %d]\n", bdmq_length(job->goMQ)));

  BDASSERT(BDMPI_Comm_free(&BDMPI_COMM_SELF) == BDMPI_SUCCESS);
  BDASSERT(BDMPI_Comm_free(&BDMPI_COMM_CWORLD) == BDMPI_SUCCESS);

  /* send a message to the slave telling it that you are leaving... */
  donemsg.msgtype = BDMPI_MSGTYPE_FINALIZE;
  donemsg.myrank  = job->rank;

  if (bdmq_send(job->reqMQ, &donemsg, sizeof(bdmsg_t)) == -1)
    bdprintf("Failed on sending a donemsg: %s.\n", strerror(errno));

  /* turn off sbmalloc */
  /* TODO: slave should subtract any memory it has from job->memrss */
  sb_finalize();

  /* wait for a go response from the master */
  for (;;) {
    if (-1 == bdmq_recv(job->goMQ, &gomsg, sizeof(bdmsg_t)))
      bdprintf("Failed on trying to recv a go message: %s.\n", strerror(errno));
    if (BDMPI_MSGTYPE_PROCEED == gomsg.msgtype)
      break;
    slv_route(job, &gomsg);
  }

  S_IFSET(BDMPI_DBG_IPCS, bdprintf("iBDMPI_Finalize: I got the following response: %d\n", gomsg.msgtype));

  /* close the global SMR */
  bdsm_close(job->globalSM);

  /* close the message queues */
  bdmq_close(job->reqMQ);
  bdmq_close(job->c2sMQ);
  bdmq_close(job->c2mMQ);
  bdmq_close(job->goMQ);

  /* close the mutexes */
  bdlock_close(job->mlockMX);
  bdlock_close(job->criticalMX);

  /* close the SCM for master-slave communication */
  bdscb_close(job->scb);

  /* free communicators */
  gk_free((void **)&BDMPI_COMM_WORLD, &BDMPI_COMM_NODE, LTERM);

  return BDMPI_SUCCESS;
}
