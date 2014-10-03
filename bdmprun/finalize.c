/*!
\file
\brief Master-side BDMPI_Finalize operation.
\date Started 4/6/2013
\author George
*/


#include "bdmprun.h"


/*************************************************************************/
/*! Response to a BDMPI_Finalize().
    Increases the number of slaves in the pool and if all have joined,
    makes all of them runnable.
*/
/*************************************************************************/
void mstr_finalize(mjob_t *job, bdmsg_t *msg)
{
  size_t i;

  /* block it */
  slvpool_cblock(job, msg->myrank-job->jdesc->soffset);

  BD_GET_LOCK(job->schedule_lock);

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_finalize: njoined: %d [entering]\n",
        job->mynode, msg->myrank, job->njoined));

  /* decreament the njoined counter and take action accordingly */
  if (--job->njoined == 0) {
    /* all processes have called the init, unblock them and remove the co-operative
     * scheduling constraint. */
    job->nR = job->ns;
    for (i=0; i<job->ns; i++)
      slvpool_cunblock(job, i);
  }

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_finalize: njoined: %d [exiting]\n",
        job->mynode, msg->myrank, job->njoined));

  BD_LET_LOCK(job->schedule_lock);

  return;
}
