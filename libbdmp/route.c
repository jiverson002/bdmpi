/*!
 * \file    route.c
 * \brief   Implements the route operation for the slave nodes
 * \date    Started 9/22/2014
 * \author  Jeremy
 */


#include "bdmplib.h"


/*************************************************************************/
/*! Routes the go queue messages from the master to the appropriate
    handler */
/*************************************************************************/
void slv_route(sjob_t * const job, bdmsg_t const * const gomsg)
{
  int done=1;

  S_IFSET(BDMPI_DBG_IPCS, bdprintf("[%04d] slv_route: response: %d\n",
        job->rank, gomsg->msgtype));

  switch (gomsg->msgtype) {
    case BDMPI_MSGTYPE_MEMFREE:
      if (job->jdesc->nr < job->jdesc->ns)
        sb_saveall();
      bdmq_send(job->c2mMQ, &done, sizeof(int));
      break;

    default:
      errexit("Got go queue response %d\n", gomsg->tag);
  }

  return;
}
