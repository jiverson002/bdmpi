/*!
 *  \file   memory.c
 *  \brief  Various functions controlling memory usage.
 *  \date   Started 9/22/2014
 *  \author Jeremy
 */


#include "bdmprun.h"


/*************************************************************************/
/*! Response to a BDMPI_MSGTYPE_MEMRQST
    Checks if a the memory being requested will over-subscribe the system
    memory, in which case it chooses a slave to wake and release its memory,
    otherwise it notifies the requesting slave that it can safely allocate the
    amount of memory requested.
*/
/*************************************************************************/
void * mstr_mem_rqst(void * const arg)
{
  mjob_t * const job = ((ptarg_t const *)arg)->job;
  bdmsg_t const * const msg = &(((ptarg_t const *)arg)->msg);
  int const source = msg->source;
  size_t const count = msg->count;
  bdmsg_t gomsg;

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_mem_rqst: count: "
    "%zu [entering]\n", job->mynode, msg->source, msg->count));

  BD_GET_LOCK(job->memory_lock);
  job->memrss += count;
  BD_LET_LOCK(job->memory_lock);

  memory_wakeup_some(job);

  gomsg.msgtype = BDMPI_MSGTYPE_PROCEED;
  if (-1 == bdmq_send(job->goMQs[source], &gomsg, sizeof(bdmsg_t)))
    bdprintf("Failed to send a go message to %d: %s\n", source, strerror(errno));

  gk_free((void **)&arg, LTERM);

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_mem_rqst: count: "
    "%zu [exiting]\n", job->mynode, msg->source, msg->count));

  return NULL;
}


/*************************************************************************/
/*! Response to a BDMPI_MSGTYPE_MEMRLSD
    Checks if a the memory being requested will over-subscribe the system
    memory, in which case it chooses a slave to wake and release its memory,
    otherwise it notifies the requesting slave that it can safely allocate the
    amount of memory requested.
*/
/*************************************************************************/
void * mstr_mem_rlsd(void * const arg)
{
  mjob_t * const job = ((ptarg_t const *)arg)->job;
  bdmsg_t const * const msg = &(((ptarg_t const *)arg)->msg);
  int const source = msg->source;
  size_t const count = msg->count;
  bdmsg_t gomsg;

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_mem_rlsd: count: "
    "%zu [entering]\n", job->mynode, msg->source, msg->count));

  BD_GET_LOCK(job->memory_lock);
  job->memrss -= count;
  BD_LET_LOCK(job->memory_lock);

  gk_free((void **)&arg, LTERM);

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_mem_rlsd: count: "
    "%zu [exiting]\n", job->mynode, msg->source, msg->count));

  return NULL;
}


/*************************************************************************/
/*! If the number of running slaves is less that nr, then some runnable
    slaves are told to go ahead and execute. */
/*************************************************************************/
void memory_wakeup_some(mjob_t * const job)
{
  int i, itogo, togo;
  size_t memrss;
  bdmsg_t msg, gomsg;

  BD_GET_LOCK(job->memory_lock);
  memrss = job->memrss;
  BD_LET_LOCK(job->memory_lock);

  gomsg.msgtype = BDMPI_MSGTYPE_MEMFREE;

  BD_GET_LOCK(job->schedule_lock);

  while ((job->nrunnable+job->nmblocked+job->ncblocked) > 0 &&
         memrss > job->memmax)
  {
    /* select a blocked slave to wakeup for memory free'ing */
    itogo = slvpool_select_task_to_wakeup(job, BDMPRUN_WAKEUP_VRSS);
    /*itogo = slvpool_select_task_to_wakeup(job, BDMPRUN_WAKEUP_MINRSS);*/

    if (itogo < job->nrunnable)
      togo = job->runnablelist[itogo];
    else if (itogo < job->nmblocked)
      togo = job->mblockedlist[itogo-job->nrunnable];
    else
      togo = job->cblockedlist[itogo-job->nrunnable-job->nmblocked];

    /* send that slave a go message */
    M_IFSET(BDMPI_DBG_IPCM,
        bdprintf("[MSTR%04d] Telling slave [%d:%d] to free memory [msg:%d]\n",
                 job->mynode, togo, (int)job->spids[togo], gomsg.tag));
    if (-1 == bdmq_send(job->goMQs[togo], &gomsg, sizeof(bdmsg_t)))
      bdprintf("Failed to send a go message to %d: %s\n", togo, strerror(errno));

    BD_GET_LOCK(job->memory_lock);
    if (memrss == job->memrss) {
      BD_LET_LOCK(job->memory_lock);
      break;
    }
    memrss = job->memrss;
    BD_LET_LOCK(job->memory_lock);
  }

  BD_LET_LOCK(job->schedule_lock);
}
