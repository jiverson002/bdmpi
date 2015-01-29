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
    amount of memory requested. */
/*************************************************************************/
void * mstr_mem_load(void * const arg)
{
  mjob_t * const job = ((ptarg_t const *)arg)->job;
  bdmsg_t const * const msg = &(((ptarg_t const *)arg)->msg);
  bdmsg_t gomsg;

  BD_GET_LOCK(job->schedule_lock);
  job->memrss += msg->count;
  job->slvrss[msg->source] += msg->count;

  //printf("LOAD (%d) %10zu / %10zu / %10zu\n", msg->source, msg->count, job->memrss, job->memmax);
  //fflush(stdout);

  if (job->memrss > job->memmax) {
    //printf("  WAKE (%d)\n", msg->source);
    //fflush(stdout);
    memory_wakeup_some(job, msg->source, msg->count);
  }

  gomsg.msgtype = BDMPI_MSGTYPE_PROCEED;
  if (-1 == bdmq_send(job->goMQs[msg->source], &gomsg, sizeof(bdmsg_t)))
    bdprintf("Failed to send a go message to %d: %s\n", msg->source, strerror(errno));
  BD_LET_LOCK(job->schedule_lock);

  gk_free((void **)&arg, LTERM);

  return NULL;
}


/*************************************************************************/
/*! Response to a BDMPI_MSGTYPE_MEMSAVE
    Decreases resident set size for job. */
/*************************************************************************/
void * mstr_mem_save(void * const arg)
{
  mjob_t * const job = ((ptarg_t const *)arg)->job;
  bdmsg_t const * const msg = &(((ptarg_t const *)arg)->msg);
  bdmsg_t gomsg;

  BD_GET_LOCK(job->schedule_lock);
  job->memrss -= msg->count;
  job->slvrss[msg->source] -= msg->count;

  //printf("SAVE (%d) %10zu / %10zu / %10zu\n", msg->source, msg->count, job->memrss, job->memmax);
  //fflush(stdout);

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
void memory_wakeup_some(mjob_t * const job, int const source,
                        size_t const size)
{
  int i, itogo, iitogo, togo=0, type=0;
  size_t count=0, resident, ires;
  bdmsg_t msg;

  msg.msgtype = BDMPI_MSGTYPE_MEMFREE;

  while (job->memrss > job->memmax &&
         0 != job->nrunnable+job->nmblocked+job->ncblocked)
  {
#if 0
    itogo = memory_select_task_to_wakeup(job, BDMPRUN_WAKEUP_VRSS);
#else
    for (itogo=-1,ires=0,i=0; i<job->nrunnable; i++) {
      resident = job->slvrss[job->runnablelist[i]];

      if (source != job->runnablelist[i] && resident > ires) {
        itogo = i;
        ires = resident;
      }
    }
    for (i=0; i<job->nmblocked; i++) {
      resident = job->slvrss[job->mblockedlist[i]];

      if (source != job->mblockedlist[i] && resident > ires) {
        itogo = i+job->nrunnable;
        ires = resident;
      }
    }
    for (i=0; i<job->ncblocked; i++) {
      resident = job->slvrss[job->cblockedlist[i]];

      if (source != job->cblockedlist[i] && resident > ires) {
        itogo = i+job->nrunnable+job->nmblocked;
        ires = resident;
      }
    }
#endif

    if (-1 == itogo) {
      /*printf("  NOMEM [");
      for (i=0; i<job->nrunnable; ++i)
        printf("%d ", job->runnablelist[i]);
      printf("| ");
      for (i=0; i<job->nmblocked; ++i)
        printf("%d ", job->mblockedlist[i]);
      printf("| ");
      for (i=0; i<job->ncblocked; ++i)
        printf("%d ", job->cblockedlist[i]);
      printf("] (");
      for (i=0; i<4; ++i)
        printf("%zu ", job->slvrss[i]);
      printf(")\n");
      fflush(stdout);*/
      break;
    }

#if 0
    if (itogo < job->nrunnable) {
      iitogo = itogo;
      togo = job->runnablelist[iitogo];
      job->runnablelist[iitogo] = job->runnablelist[--job->nrunnable];
      job->runnablemap[togo] = -1;
      type = 0;
    }
    else if (itogo < job->nrunnable+job->nmblocked) {
      iitogo = itogo - job->nrunnable;
      togo = job->mblockedlist[iitogo];
      job->mblockedlist[iitogo] = job->mblockedlist[--job->nmblocked];
      job->mblockedmap[togo] = -1;
      type = 1;
    }
    else {
      iitogo = itogo - job->nrunnable - job->nmblocked;
      togo = job->cblockedlist[iitogo];
      job->cblockedlist[iitogo] = job->cblockedlist[--job->ncblocked];
      job->cblockedmap[togo] = -1;
      type = 2;
    }
    BD_LET_LOCK(job->schedule_lock);
#else
    //for (i=0; i<4; ++i)
    //  printf("[%zu] ", job->slvrss[i]);
    if (itogo < job->nrunnable) {
      iitogo = itogo;
      togo = job->runnablelist[iitogo];
    }
    else if (itogo < job->nrunnable+job->nmblocked) {
      iitogo = itogo-job->nrunnable;
      togo = job->mblockedlist[iitogo];
    }
    else {
      iitogo = itogo-job->nrunnable-job->nmblocked;
      togo = job->cblockedlist[iitogo];
    }
#endif

#if 1
    //printf("WAKE (%d) %d %d %d\n", togo, source, itogo, iitogo);
    //fflush(stdout);

    //printf("telling %d to free memory\n", togo);
    /* send that slave a go free memory message */
    if (-1 == bdmq_send(job->goMQs[togo], &msg, sizeof(bdmsg_t)))
      bdprintf("Failed to send a go message to %d: %s\n", togo, strerror(errno));
    /* recv from slave a done message */
    if (-1 == bdmq_recv(job->c2mMQs[togo], &count, sizeof(size_t)))
      bdprintf("Failed to recv a done message from %d: %s\n", togo, strerror(errno));
    //printf("%d free'd %zu bytes of memory\n", togo, count);
#endif

#if 0
    BD_GET_LOCK(job->schedule_lock);
    if (0 == type) {
      job->runnablelist[job->nrunnable] = togo;
      job->runnablemap[togo] = job->nrunnable++;
    }
    else if (1 == type) {
      job->mblockedlist[job->nmblocked] = togo;
      job->mblockedmap[togo] = job->nmblocked++;
    }
    else if (2 == type) {
      job->cblockedlist[job->ncblocked] = togo;
      job->cblockedmap[togo] = job->ncblocked++;
    }
#endif

#if 1
    if (0 == count)
      break;
#endif

    job->memrss -= count;
    job->slvrss[togo] -= count;

    //printf("SAVE (%d) %10zu / %10zu / %10zu\n", togo, count, job->memrss, job->memmax);
    //fflush(stdout);
  }
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
    case BDMPRUN_WAKEUP_VRSS:
      ifres = 0;
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
        if (0 != resident && cfres > ifres) {
          itogo = i;
          ifres = cfres;
        }
      }
      for (i=0; i<job->nmblocked; i++) {
        sprintf(fname, "/proc/%d/statm", job->spids[job->mblockedlist[i]]);

        if (NULL == (fp = fopen(fname, "r")))
          slvpool_abort(1, "Failed to open %s.\n", fname);
        if (2 != fscanf(fp, "%zu %zu", &size, &resident))
          slvpool_abort(1, "Failed to read to values from %s.\n", fname);
        if (0 != fclose(fp))
          slvpool_abort(1, "Failed to close %s.\n", fname);
        resident = job->slvrss[job->mblockedlist[i]];

        cfres = 1.0*resident/(size*sysconf(_SC_PAGESIZE));
        if (0 != resident && cfres > ifres) {
          itogo = i+job->nrunnable;
          ifres = cfres;
        }
      }
      for (i=0; i<job->ncblocked; i++) {
        sprintf(fname, "/proc/%d/statm", job->spids[job->cblockedlist[i]]);

        if (NULL == (fp = fopen(fname, "r")))
          slvpool_abort(1, "Failed to open %s.\n", fname);
        if (2 != fscanf(fp, "%zu %zu", &size, &resident))
          slvpool_abort(1, "Failed to read to values from %s.\n", fname);
        if (0 != fclose(fp))
          slvpool_abort(1, "Failed to close %s.\n", fname);
        resident = job->slvrss[job->cblockedlist[i]];

        cfres = 1.0*resident/(size*sysconf(_SC_PAGESIZE));
        if (0 != resident && cfres > ifres) {
          itogo = i+job->nrunnable+job->nmblocked;
          ifres = cfres;
        }
      }
      break;

    default:
      itogo = 0;
  }

  return itogo;
}
