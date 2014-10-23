/*!
\file
\brief The master-side BDMPI_Barrier operation.
\date Started 4/6/2013
\author George
*/


#include "bdmprun.h"


/*************************************************************************/
/*! Response to a BDMPI_Barrier.
    Protocol:
       Blocks the process.
       Decreases the counter associated with the counter.
       If counter becomes 0, then moves all processes to runnable state
       and sets the counter back to the size of that communicator.
*/
/*************************************************************************/
void *mstr_barrier(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  size_t i;

  BD_GET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_barrier: counter: %d [entering]\n",
                job->mynode, msg->myrank, job->comms[msg->mcomm]->counter));

  /* block the slave */
  slvpool_cblock(job, babel_get_srank(job->comms[msg->mcomm], msg->myrank));

  /* see if all slaves have called the barrier */
  if (--job->comms[msg->mcomm]->counter == 0) {
    if (job->comms[msg->mcomm]->nnodes > 1) {
      BD_LET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */
      MPI_Barrier(job->comms[msg->mcomm]->mpi_comm);
      M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_barrier: counter: %d [done with MPI_Barrier]\n",
                  job->mynode, msg->myrank, job->comms[msg->mcomm]->counter));
      BD_GET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */
    }

    /* unblock slaves */
    job->comms[msg->mcomm]->counter = job->comms[msg->mcomm]->lsize;
    for (i=0; i<job->comms[msg->mcomm]->lsize; i++)
      slvpool_cunblock(job, job->comms[msg->mcomm]->sranks[i]);
  }

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_barrier: counter: %d [exiting]\n",
                job->mynode, msg->myrank, job->comms[msg->mcomm]->counter));

  BD_LET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

  gk_free((void **)&arg, LTERM);

  return NULL;
}
