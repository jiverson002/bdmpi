/*!
\file
\brief Master-side BDMPI_Init implementation.
\date Started 4/6/2013
\author George
*/


#include "bdmprun.h"


/*************************************************************************/
/*! Response to a completion of BDMPI_Init().
    - Increases the number of slaves in the pool.
    - When all have joined, it unblocks them and makes them runnable.
    It assumes that prior to initialization, all slaves are runnable and
    unblocked. The call to _Init() creates the co-operating execution.
*/
/*************************************************************************/
void mstr_init(mjob_t *job, bdmsg_t *msg)
{
  size_t i;

  /* block it */
  slvpool_cblock(job, msg->myrank-job->jdesc->soffset);

  BD_GET_LOCK(job->schedule_lock);

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_init: njoined: %d [entering]\n",
              job->mynode, msg->myrank, job->njoined));


  /* increament the njoined counter and take action accordingly */
  if (++job->njoined == job->ns) {
    /* all processes have called the init, unblock them and turn on the co-operating
     * scheduling constraint. */
    job->nR = job->nr_input;
    for (i=0; i<job->ns; i++)
      slvpool_cunblock(job, i);
  }

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_init: njoined: %d [exiting]\n",
              job->mynode, msg->myrank, job->njoined));

  BD_LET_LOCK(job->schedule_lock);

  return;
}
