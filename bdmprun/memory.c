/*!
 *  \file   memory.c
 *  \brief  Various functions controlling memory usage.
 *  \date   Started 9/22/2014
 *  \author Jeremy
 */


#include "bdmprun.h"


/*************************************************************************/
/*! Response to a BDMPI_MSGTYPE_MEMRQST or BDMPI_MSGTYPE_MEMLOAD
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

  BD_GET_LOCK(job->schedule_lock);
  BD_GET_LOCK(job->memory_lock);

  job->memrss += count;

  switch (msg->msgtype) {
  case BDMPI_MSGTYPE_MEMRQST:
    job->slvtot[source] += count;
  case BDMPI_MSGTYPE_MEMLOAD:
    job->slvrss[source] += count;
    break;
  default:
    bdprintf("[MSTR%04d.%04d] mstr_mem_rqst: invalid message %d\n",
      job->mynode, source, msg->msgtype);
    break;
  }

  memory_wakeup_some(job);

  BD_LET_LOCK(job->memory_lock);
  BD_LET_LOCK(job->schedule_lock);

  gomsg.msgtype = BDMPI_MSGTYPE_PROCEED;
  if (-1 == bdmq_send(job->goMQs[source], &gomsg, sizeof(bdmsg_t)))
    bdprintf("Failed to send a go message to %d: %s\n", source, strerror(errno));

  gk_free((void **)&arg, LTERM);

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_mem_rqst: count: "
    "%zu [exiting]\n", job->mynode, msg->source, msg->count));

  return NULL;
}


/*************************************************************************/
/*! Response to a BDMPI_MSGTYPE_MEMRLSD or BDMPI_MSGTYPE_MEMSAVE
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

  //BDASSERT(job->memrss >= count);
  if (job->memrss < count)
    bdprintf("mstr_mem_rlsd.0\n");
  job->memrss -= count;

  switch (msg->msgtype) {
  case BDMPI_MSGTYPE_MEMRLSD:
    //BDASSERT(job->slvtot[source] >= count);
    if (job->slvtot[source] < count)
      bdprintf("mstr_mem_rlsd.1\n");
    job->slvtot[source] -= count;
  case BDMPI_MSGTYPE_MEMSAVE:
    //BDASSERT(job->slvrss[source] >= count);
    if (job->slvrss[source] < count)
      bdprintf("mstr_mem_rlsd.2\n");
    job->slvrss[source] -= count;
    break;
  default:
    bdprintf("[MSTR%04d.%04d] mstr_mem_rlsd: invalid message %d\n",
      job->mynode, source, msg->msgtype);
    break;
  }

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
  size_t count;
  bdmsg_t msg;

  msg.msgtype = BDMPI_MSGTYPE_MEMFREE;

  BD_GET_LOCK(job->schedule_lock);
  BD_GET_LOCK(job->memory_lock);

  while ((job->nrunnable+job->nmblocked+job->ncblocked) > 0 &&
         job->memrss > job->memmax)
  {
    /* select a blocked slave to wakeup for memory free'ing */
    itogo = memory_select_task_to_wakeup(job, BDMPRUN_WAKEUP_VMEM);

    if (itogo < job->nrunnable)
      togo = job->runnablelist[itogo];
    else if (itogo < job->nmblocked)
      togo = job->mblockedlist[itogo-job->nrunnable];
    else
      togo = job->cblockedlist[itogo-job->nrunnable-job->nmblocked];

    /* send that slave a go message */
    if (-1 == bdmq_send(job->goMQs[togo], &msg, sizeof(bdmsg_t)))
      bdprintf("Failed to send a go message to %d: %s\n", togo, strerror(errno));
    if (-1 == bdmq_recv(job->c2mMQs[togo], &count, sizeof(size_t)))
      bdprintf("Failed to recv a done message to %d: %s\n", togo, strerror(errno));

    if (0 == count)
      break;

    if (job->memrss < count)
      bdprintf("memory_wakeup_some.0 %d %zu %zu\n", togo, job->memrss, count);
    job->memrss -= count;
    if (job->slvrss[togo] < count)
      bdprintf("memory_wakeup_some.1 %d %zu %zu\n", togo, job->slvrss[togo], count);
    job->slvrss[togo] -= count;
  }

  BD_LET_LOCK(job->memory_lock);
  BD_LET_LOCK(job->schedule_lock);
}


/*************************************************************************/
/*! Selects a runnable task to wake up.
    Returns the index within the a list of the selected task. */
/*************************************************************************/
int memory_select_task_to_wakeup(mjob_t *job, int type)
{
  int i, itogo;
  size_t size, resident;
  float ifres, cfres;
  char fname[1024];
  FILE *fp;

  switch (type) {
    case BDMPRUN_WAKEUP_VMEM:
      ifres = -1.0;
      itogo = 0;
      for (i=0; i<job->nrunnable; i++) {
        resident = job->slvrss[job->runnablelist[i]];
        size = job->slvtot[job->runnablelist[i]];
        cfres = 1.0*resident/size;
        if (cfres > ifres) {
          itogo = i;
          ifres = cfres;
        }
      }
      for (i=0; i<job->ncblocked; i++) {
        resident = job->slvrss[job->cblockedlist[i]];
        size = job->slvtot[job->cblockedlist[i]];
        cfres = 1.0*resident/size;
        if (cfres > ifres) {
          itogo = i + job->nrunnable;
          ifres = cfres;
        }
      }
      for (i=0; i<job->nmblocked; i++) {
        resident = job->slvrss[job->mblockedlist[i]];
        size = job->slvtot[job->mblockedlist[i]];
        cfres = 1.0*resident/size;
        if (cfres > ifres) {
          itogo = i + job->nrunnable + job->ncblocked;
          ifres = cfres;
        }
      }
      break;

    default:
      itogo = 0;
  }

  return itogo;
}
