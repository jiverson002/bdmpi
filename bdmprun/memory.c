/*!
 *  \file   memory.c
 *  \brief  Various functions controlling memory usage.
 *  \date   Started 9/22/2014
 *  \author Jeremy
 */


#include "bdmprun.h"


/*************************************************************************/
/*! Response to a BDMPI_MSGTYPE_MEMLOAD
    Checks if a the memory being requested will over-subscribe the system
    memory, in which case it chooses a slave to wake and release its memory,
    otherwise it notifies the requesting slave that it can safely load the
    amount of memory requested.
*/
/*************************************************************************/
void * mstr_mem_load(void * const arg)
{
  mjob_t * const job = ((ptarg_t const *)arg)->job;
  bdmsg_t const * const msg = &(((ptarg_t const *)arg)->msg);
  bdmsg_t gomsg;

  BD_GET_LOCK(job->schedule_lock);
  //bdprintf("[%04d] mstr_mem_load %10zu / %10zu / %10zu\n", msg->source,
  //  msg->count, job->memrss, job->memmax);

  job->memrss += msg->count;
  job->slvrss[msg->source] += msg->count;

  if (job->memrss > job->memmax)
    memory_wakeup_some(job, msg->count);

  gomsg.msgtype = BDMPI_MSGTYPE_PROCEED;
  if (-1 == bdmq_send(job->goMQs[msg->source], &gomsg, sizeof(bdmsg_t)))
    bdprintf("Failed to send a go message to %d: %s\n", msg->source, strerror(errno));
  BD_LET_LOCK(job->schedule_lock);

  gk_free((void **)&arg, LTERM);

  return NULL;
}


/*************************************************************************/
/*! Response to a BDMPI_MSGTYPE_MEMSAVE
    Decreases resident set size for job.
*/
/*************************************************************************/
void * mstr_mem_save(void * const arg)
{
  mjob_t * const job = ((ptarg_t const *)arg)->job;
  bdmsg_t const * const msg = &(((ptarg_t const *)arg)->msg);
  bdmsg_t gomsg;

  BD_GET_LOCK(job->schedule_lock);
  //bdprintf("[%04d] mstr_mem_save %10zu / %10zu / %10zu\n", msg->source,
  //  msg->count, job->memrss, job->memmax);

  job->memrss -= msg->count;
  job->slvrss[msg->source] -= msg->count;

  gomsg.msgtype = BDMPI_MSGTYPE_PROCEED;
  if (-1 == bdmq_send(job->goMQs[msg->source], &gomsg, sizeof(bdmsg_t)))
    bdprintf("Failed to send a go message to %d: %s\n", msg->source, strerror(errno));
  BD_LET_LOCK(job->schedule_lock);

  gk_free((void **)&arg, LTERM);

  return NULL;
}


/*************************************************************************/
/*! If the number of running slaves is less that nr, then some runnable
    slaves are told to go ahead and execute. */
/*************************************************************************/
void memory_wakeup_some(mjob_t * const job, size_t const size)
{
  int i, itogo, togo=0;
  size_t count=0;
  bdmsg_t msg;

  msg.msgtype = BDMPI_MSGTYPE_MEMFREE;

  BD_GET_LOCK(job->schedule_lock);
  if (0 != job->nrunnable) {
    itogo = memory_select_task_to_wakeup(job, BDMPRUN_WAKEUP_VRSS);

    if (-1 != itogo) {
      togo = job->runnablelist[itogo];
      job->runnablelist[itogo] = job->runnablelist[--job->nrunnable];
      job->runnablemap[togo] = -1;
      BD_LET_LOCK(job->schedule_lock);
#if 1

      /* send that slave a go free memory message */
      if (-1 == bdmq_send(job->goMQs[togo], &msg, sizeof(bdmsg_t)))
        bdprintf("Failed to send a go message to %d: %s\n", togo, strerror(errno));
      /* recv from slave a done message */
      if (-1 == bdmq_recv(job->c2mMQs[togo], &count, sizeof(size_t)))
        bdprintf("Failed to recv a done message from %d: %s\n", togo, strerror(errno));

#endif
      BD_GET_LOCK(job->schedule_lock);
      job->memrss -= count;
      job->slvrss[togo] -= count;

      job->runnablelist[job->nrunnable] = togo;
      job->runnablemap[togo] = job->nrunnable++;
    }
  }
  BD_LET_LOCK(job->schedule_lock);
}


/*************************************************************************/
/*! Selects a runnable task to wake up.
    Returns the index within the a list of the selected task. */
/*************************************************************************/
int memory_select_task_to_wakeup(mjob_t *job, int type)
{
  int i, itogo;
  size_t size, resident, isize;
  float ifres, cfres;
  char fname[1024];
  FILE *fp;

  BD_GET_LOCK(job->schedule_lock);
  switch (type) {
    case BDMPRUN_WAKEUP_VRSS:
      ifres = 0;
      isize = 0;
      itogo = -1;
      for (i=0; i<job->nrunnable; i++) {
        sprintf(fname, "/proc/%d/statm", job->spids[job->runnablelist[i]]);

        if (NULL == (fp = fopen(fname, "r")))
          slvpool_abort(1, "Failed to open %s.\n", fname);
        if (2 != fscanf(fp, "%zu %zu", &size, &resident))
          slvpool_abort(1, "Failed to read to values from %s.\n", fname);
        if (0 != fclose(fp))
          slvpool_abort(1, "Failed to close %s.\n", fname);
        resident = job->slvrss[job->runnablelist[i]];

        cfres = 1.0*resident/(size*sysconf(_SC_PAGESIZE));
        //bdprintf("%s %10zu %10zu %5.4f\n", fname, size, resident, cfres);
        if (cfres > ifres) {
          itogo = i;
          ifres = cfres;
          isize = size;
        }
      }

      if (-1 != itogo) {
        bdprintf("[%04d] memory_wakeup_some\n", job->runnablelist[itogo]);
        bdprintf("[%04d] %10zu / %10zu\n", job->runnablelist[itogo],
          job->slvrss[job->runnablelist[itogo]], isize*sysconf(_SC_PAGESIZE));
      }

      break;

    default:
      itogo = 0;
  }
  BD_LET_LOCK(job->schedule_lock);

  return itogo;
}
