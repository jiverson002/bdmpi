/*!
\file
\brief Various functions from listening and coordinating the pool of slaves.
\date Started 4/6/2013
\author George
*/


#include "bdmprun.h"



/*************************************************************************/
/*! It repeatedly checks/routes messages from the slaves and watches for
    slaves that have terminated. */
/*************************************************************************/
void slvpool_listen(mjob_t *job)
{
  int status;
  ssize_t len;
  pid_t cpid;
  bdmsg_t msg;

  /* get into the loop */
  for (;;) {
    //slvpool_check_deadlock(job);

    M_IFSET(BDMPI_DBG_SLVPOOL, slvpool_showstats(job));

    /* get into a loop to check if there are any exited slave processes */
    while ((cpid = waitpid(-1, &status, WNOHANG|WUNTRACED)) > 0) {
      if (WIFEXITED(status)) {
        M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] slvpool_listen: cpid: %d terminated with exit.\n",
            job->mynode, (int)cpid));
        slvpool_remove_pid(job, cpid);
        if (WEXITSTATUS(status) != EXIT_SUCCESS)
          slvpool_kill_all(job);
      }
      else if (WIFSIGNALED(status)) {
        M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] slvpool_listen: cpid: %d terminated with signal.\n",
            job->mynode, (int)cpid));
        slvpool_abort(1, "Process %d killed by signal %d.\n", (int)cpid, WTERMSIG(status));
      } else if (WIFSTOPPED(status)) {
        M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] slvpool_listen: cpid: %d terminated with stop.\n",
            job->mynode, (int)cpid));
        slvpool_abort(1, "Process %d stopped by signal %d.\n", (int)cpid, WSTOPSIG(status));
      }
      else {
        M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] slvpool_listen: cpid: %d exited with status: %d\n",
            job->mynode, (int)cpid, status));
      }
    }

    /* this is what gets us out of the loop */
    if (job->nalive == 0)
      return;

    /* see if there is room to schedule some more */
    slvpool_wakeup_some(job);

    /* service any requests to the master node */
    if (job->mynode == 0 && job->nnodes > 1)
      mnode_service_xcomms(job);

    /* service any inter-node requests */
    slvpool_service_xcomms(job);

    /* see if there is room to schedule some more */
    slvpool_wakeup_some(job);

    /* wait a bit for some messages */
    len = bdmq_timedrecv(job->reqMQ, &msg, sizeof(bdmsg_t), 100000);
    if (len == -1) {
      if (errno != ETIMEDOUT) {
        bdprintf("Failed to timedrecv(job->reqMQ). EBADF: %d, EINTR: %d, EINVAL: %d, EMSGSIZE:%d, %s\n",
            (errno == EBADF), (errno == EINTR), (errno == EINVAL), (errno == EMSGSIZE),
            strerror(errno));
        slvpool_kill_all(job);
      }
    }
    else {
      if (len != sizeof(bdmsg_t)) {
        bdprintf("Failed on timedrecv(job->reqMQ). Incorrect msgsize: %zd %zd\n", sizeof(bdmsg_t), len);
        slvpool_kill_all(job);
      }
      slvpool_route(job, &msg);
    }
  }

  return;
}


/*************************************************************************/
/*! Displays various statistics related to the slave pool
*/
/*************************************************************************/
void slvpool_showstats(mjob_t *job)
{
  int i, len;
  header_t *curr;

  BD_GET_LOCK(job->schedule_lock);

  bdprintf("[MSTR%04d] nalive: %d, nrunning: %d, nrunnable: %d, n[m/c]blocked: %d/%d\n",
      job->mynode, job->nalive, job->nrunning, job->nrunnable, job->nmblocked, job->ncblocked);

  bdprintf("[MSTR%04d] %6s %3s ARRBB %4s Headers\n", job->mynode, "SPID", "RNK", "SQL");
  for (i=0; i<job->ns; i++) {
    for (curr=job->psends[i], len=0; curr!=NULL; curr=curr->next, len++);

    bdprintf("[MSTR%04d] %6d %3d %d%d%d%d%d %4d",
        job->mynode,
        (int)job->spids[i], (int)i,
        (job->alivemap[i] == -1 ? 0 : 1),
        (job->runningmap[i] == -1 ? 0 : 1),
        (job->runnablemap[i] == -1 ? 0 : 1),
        (job->mblockedmap[i] == -1 ? 0 : 1),
        (job->cblockedmap[i] == -1 ? 0 : 1),
        len);

    for (curr=job->psends[i]; curr!=NULL; curr=curr->next)
      bdprintf(" [%3d %3d %3d]", curr->msg.source, curr->msg.tag, curr->msg.mcomm);
    bdprintf("\n");
  }

  BD_LET_LOCK(job->schedule_lock);

}


/*************************************************************************/
/*! Checks if a deadlock may have occured.
*/
/*************************************************************************/
void slvpool_check_deadlock(mjob_t *job)
{
  int i;

  if (job->nrunning == 0 && job->nrunnable == 0) {
    bdprintf("[MSTR%04d] nalive: %d, nrunning: %d, nrunnable: %d, n[m/c]blocked: %d/%d\n",
        job->mynode, job->nalive, job->nrunning, job->nrunnable, job->nmblocked, job->ncblocked);

    /*
    for (i=0; i<job->ns; i++) {
      if (job->mblockedmap[i] != -1 && job->psends[i] != NULL)
        bdprintf("[MSTR%04d] DEADLOCK %6d %3d %d%d%d%d%d\n",
            job->mynode,
            (int)job->spids[i], (int)i,
            (job->alivemap[i] == -1 ? 0 : 1),
            (job->runningmap[i] == -1 ? 0 : 1),
            (job->runnablemap[i] == -1 ? 0 : 1),
            (job->mblockedmap[i] == -1 ? 0 : 1),
            (job->cblockedmap[i] == -1 ? 0 : 1));
    }
    */
  }
}


/*************************************************************************/
/*! Routes the messages from the slave to the appropriate handler
*/
/*************************************************************************/
void slvpool_route(mjob_t *job, bdmsg_t *msg)
{
  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] slvpool_route: msgtype: %d from: %d\n",
        job->mynode, msg->msgtype, msg->myrank));

  gk_startwctimer(job->routeTmr);

  if (msg->msgtype == BDMPI_MSGTYPE_INIT) {
    mstr_init(job, msg);
  }
  else if (msg->msgtype == BDMPI_MSGTYPE_FINALIZE) {
    mstr_finalize(job, msg);
  }
  else {
    pthread_t thread;
    ptarg_t *arg;

    arg = (ptarg_t *)gk_malloc(sizeof(ptarg_t), "route: arg");
    arg->job = job;
    arg->msg = *msg;

    switch (msg->msgtype) {
      case BDMPI_MSGTYPE_SEND:
        BDASSERT(pthread_create(&thread, NULL, mstr_send, arg) == 0);
        break;

      case BDMPI_MSGTYPE_RECV:
        BDASSERT(pthread_create(&thread, NULL, mstr_recv, arg) == 0);
        break;

      case BDMPI_MSGTYPE_IRECV:
        BDASSERT(pthread_create(&thread, NULL, mstr_irecv, arg) == 0);
        break;

      case BDMPI_MSGTYPE_PROBE:
        BDASSERT(pthread_create(&thread, NULL, mstr_probe, arg) == 0);
        break;

      case BDMPI_MSGTYPE_IPROBE:
        BDASSERT(pthread_create(&thread, NULL, mstr_iprobe, arg) == 0);
        break;

      case BDMPI_MSGTYPE_BARRIER:
        BDASSERT(pthread_create(&thread, NULL, mstr_barrier, arg) == 0);
        break;

      case BDMPI_MSGTYPE_BCASTI:
        BDASSERT(pthread_create(&thread, NULL, mstr_bcast_init, arg) == 0);
        break;

      case BDMPI_MSGTYPE_BCASTF:
        BDASSERT(pthread_create(&thread, NULL, mstr_bcast_recv, arg) == 0);
        break;

      case BDMPI_MSGTYPE_ALLGATHERI:
        BDASSERT(pthread_create(&thread, NULL, mstr_allgather_send, arg) == 0);
        break;

      case BDMPI_MSGTYPE_ALLGATHERF:
        BDASSERT(pthread_create(&thread, NULL, mstr_allgather_recv, arg) == 0);
        break;

      case BDMPI_MSGTYPE_REDUCEI:
        BDASSERT(pthread_create(&thread, NULL, mstr_reduce_send, arg) == 0);
        break;

      case BDMPI_MSGTYPE_REDUCEF:
        BDASSERT(pthread_create(&thread, NULL, mstr_reduce_recv, arg) == 0);
        break;

      case BDMPI_MSGTYPE_MERGEI:
        BDASSERT(pthread_create(&thread, NULL, mstr_merge_send, arg) == 0);
        break;

      case BDMPI_MSGTYPE_MERGEF:
        BDASSERT(pthread_create(&thread, NULL, mstr_merge_recv, arg) == 0);
        break;

      case BDMPI_MSGTYPE_ALLREDUCEI:
        BDASSERT(pthread_create(&thread, NULL, mstr_allreduce_send, arg) == 0);
        break;

      case BDMPI_MSGTYPE_ALLREDUCEF:
        BDASSERT(pthread_create(&thread, NULL, mstr_allreduce_recv, arg) == 0);
        break;

      case BDMPI_MSGTYPE_GATHERI:
        BDASSERT(pthread_create(&thread, NULL, mstr_gather_send, arg) == 0);
        break;

      case BDMPI_MSGTYPE_GATHERF:
        BDASSERT(pthread_create(&thread, NULL, mstr_gather_recv, arg) == 0);
        break;

      case BDMPI_MSGTYPE_SCATTERI:
        BDASSERT(pthread_create(&thread, NULL, mstr_scatter_send, arg) == 0);
        break;

      case BDMPI_MSGTYPE_SCATTERF:
        BDASSERT(pthread_create(&thread, NULL, mstr_scatter_recv, arg) == 0);
        break;

      case BDMPI_MSGTYPE_ALLTOALLI:
        BDASSERT(pthread_create(&thread, NULL, mstr_alltoall_send, arg) == 0);
        break;

      case BDMPI_MSGTYPE_ALLTOALLF:
        BDASSERT(pthread_create(&thread, NULL, mstr_alltoall_recv, arg) == 0);
        break;

      case BDMPI_MSGTYPE_COMMDUP:
        BDASSERT(pthread_create(&thread, NULL, mstr_comm_dup, arg) == 0);
        break;

      case BDMPI_MSGTYPE_COMMFREE:
        BDASSERT(pthread_create(&thread, NULL, mstr_comm_free, arg) == 0);
        break;

      case BDMPI_MSGTYPE_COMMSPLIT:
        BDASSERT(pthread_create(&thread, NULL, mstr_comm_split, arg) == 0);
        break;

      case BDMPI_MSGTYPE_MEMRQST:
      case BDMPI_MSGTYPE_MEMLOAD:
        BDASSERT(pthread_create(&thread, NULL, mstr_mem_rqst, arg) == 0);
        break;

      case BDMPI_MSGTYPE_MEMRLSD:
      case BDMPI_MSGTYPE_MEMSAVE:
        BDASSERT(pthread_create(&thread, NULL, mstr_mem_rlsd, arg) == 0);
        break;

      default:
        slvpool_abort(1, "Got message: %d\n", msg->msgtype);
        return;
    }

    pthread_detach(thread);
  }

  gk_stopwctimer(job->routeTmr);

  return;
}


/*************************************************************************/
/*! Services any requests from other masters, and if yes, invokes the
    appropriate handler. */
/*************************************************************************/
void slvpool_service_xcomms(mjob_t *job)
{
  int i, flag;
  MPI_Message message;
  MPI_Status status;
  pthread_t thread;
  ptarg_t *arg;

  /* Check for abort requests */
  BDASSERT(MPI_Improbe(MPI_ANY_SOURCE, BDMPI_ABORT_TAG, job->mpi_wcomm,
               &flag, &message, &status)
      == MPI_SUCCESS);

  if (flag)
    slvpool_abort(0, "I was asked by another node to abort.\n");

  /* Check for HDR requests */
  for (i=0; i<job->nnodes; i++) {
    BDASSERT(MPI_Improbe(MPI_ANY_SOURCE, BDMPI_HDR_TAG, job->mpi_wcomm,
                 &flag, &message, &status)
        == MPI_SUCCESS);

    if (!flag)
      break;

    arg = (ptarg_t *)gk_malloc(sizeof(ptarg_t), "route: arg");
    arg->job = job;

    BDASSERT(MPI_Mrecv(&(arg->msg), sizeof(bdmsg_t), MPI_BYTE, &message, &status)
        == MPI_SUCCESS);

    switch (arg->msg.msgtype) {
      case BDMPI_MSGTYPE_SEND:
        BDASSERT(pthread_create(&thread, NULL, mstr_send_remote, arg) == 0);
        break;

      default:
        slvpool_abort(1, "Got message: %d\n", arg->msg.msgtype);
        return;
    }

    pthread_detach(thread);
    //pthread_join(thread, NULL);
  }

  return;
}


/*************************************************************************/
/*! If the number of running slaves is less that nr, then some runnable
    slaves are told to go ahead and execute. */
/*************************************************************************/
void slvpool_wakeup_some(mjob_t *job)
{
  int i, itogo, togo, msg=1;
  bdmsg_t gomsg;

  gomsg.msgtype = BDMPI_MSGTYPE_PROCEED;

  BD_GET_LOCK(job->schedule_lock);

  while (job->nrunnable > 0 && job->nrunning < job->nR) {
    /* select a runnable slave */
    //itogo = slvpool_select_task_to_wakeup(job, BDMPRUN_WAKEUP_LAST);
    //itogo = slvpool_select_task_to_wakeup(job, BDMPRUN_WAKEUP_FIRST);
    //itogo = slvpool_select_task_to_wakeup(job, BDMPRUN_WAKEUP_LIFO);
    //itogo = slvpool_select_task_to_wakeup(job, BDMPRUN_WAKEUP_VRSS);
    itogo = slvpool_select_task_to_wakeup(job, BDMPRUN_WAKEUP_VMEM);

    togo = job->runnablelist[itogo];
    job->runnablelist[itogo] = job->runnablelist[--job->nrunnable];
    job->runnablemap[togo] = -1;

    /* move it to runninglist */
    job->runninglist[job->nrunning] = togo;
    job->runningmap[togo] = job->nrunning++;

    /* send that slave a go message */
    M_IFSET(BDMPI_DBG_IPCM,
        bdprintf("[MSTR%04d] Telling slave [%d:%d] to go [msg:%d]\n",
          job->mynode, togo, (int)job->spids[togo], msg));
    if (-1 == bdmq_send(job->goMQs[togo], &gomsg, sizeof(bdmsg_t)))
      bdprintf("Failed to send a go message to %d: %s\n", togo, strerror(errno));
  }

  BD_LET_LOCK(job->schedule_lock);
}


/*************************************************************************/
/*! Selects a runnable task to wake up.
    Returns the index within the runnablelist of the selected task. */
/*************************************************************************/
int slvpool_select_task_to_wakeup(mjob_t *job, int type)
{
  int i, itogo;
  size_t size, resident;
  float ifres, cfres;
  char fname[1024];
  FILE *fp;

  BD_GET_LOCK(job->schedule_lock);

  switch (type) {
    case BDMPRUN_WAKEUP_FIRST:
      itogo = 0;
      break;

    case BDMPRUN_WAKEUP_LAST:
      itogo = job->nrunnable-1;
      break;

    case BDMPRUN_WAKEUP_LIFO:
      itogo = 0;
      for (i=1; i<job->nrunnable; i++) {
        if (job->blockedts[job->runnablelist[i]] > job->blockedts[job->runnablelist[itogo]])
          itogo = i;
      }
      break;

    case BDMPRUN_WAKEUP_VRSS:
      ifres = -1.0;
      itogo = 0;
      for (i=0; i<job->nrunnable; i++) {
        sprintf(fname, "/proc/%d/statm", job->spids[job->runnablelist[i]]);

        if ((fp = fopen(fname, "r")) == NULL)
          slvpool_abort(1, "Failed to open %s.\n", fname);
        if (fscanf(fp, "%zu %zu", &size, &resident) != 2)
          slvpool_abort(1, "Failed to read to values from %s.\n", fname);
        if (fclose(fp) != 0)
          slvpool_abort(1, "Failed to close%s.\n", fname);

        cfres = 1.0*resident/size;
        //bdprintf("%s %10zu %10zu %5.4f\n", fname, size, resident, cfres);
        if (cfres > ifres) {
          itogo = i;
          ifres = cfres;
        }
      }
      break;

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
      break;

    default:
      itogo = 0;
  }

  BD_LET_LOCK(job->schedule_lock);

  return itogo;

}


/*************************************************************************/
/*! Updates the slave lists/maps to reflect that a slave has exited.
*/
/*************************************************************************/
void slvpool_remove(mjob_t *job, int rank)
{
  int i;

  BD_GET_LOCK(job->schedule_lock);

  /* remove it from runninglist */
  if (job->runningmap[rank] != -1) {
    i = job->runningmap[rank];
    job->runningmap[rank] = -1;
    job->nrunning--;
    if (job->nrunning != i) {
      job->runninglist[i] = job->runninglist[job->nrunning];
      job->runningmap[job->runninglist[i]] = i;
    }
  }

  /* remove it from runnablelist */
  if (job->runnablemap[rank] != -1) {
    i = job->runnablemap[rank];
    job->runnablemap[rank] = -1;
    job->nrunnable--;
    if (job->nrunnable != i) {
      job->runnablelist[i] = job->runnablelist[job->nrunnable];
      job->runnablemap[job->runnablelist[i]] = i;
    }
  }

  /* remove it from mblockedlist */
  if (job->mblockedmap[rank] != -1) {
    i = job->mblockedmap[rank];
    job->mblockedmap[rank] = -1;
    job->nmblocked--;
    if (job->nmblocked != i) {
      job->mblockedlist[i] = job->mblockedlist[job->nmblocked];
      job->mblockedmap[job->mblockedlist[i]] = i;
    }
  }

  /* remove it from cblockedlist */
  if (job->cblockedmap[rank] != -1) {
    i = job->cblockedmap[rank];
    job->cblockedmap[rank] = -1;
    job->ncblocked--;
    if (job->ncblocked != i) {
      job->cblockedlist[i] = job->cblockedlist[job->ncblocked];
      job->cblockedmap[job->cblockedlist[i]] = i;
    }
  }


  BD_LET_LOCK(job->schedule_lock);

  return;
}


/*************************************************************************/
/*! Updates the slave lists/maps to reflect that a slave has exited.
*/
/*************************************************************************/
void slvpool_remove_pid(mjob_t *job, pid_t cpid)
{
  int i, srank=-1;

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] slvpool_remove_pid: cpid: %d.\n",
      job->mynode, (int)cpid));

  BD_GET_LOCK(job->schedule_lock);

  for (i=0; i<job->nalive; i++) {
    srank = job->alivelist[i];
    if (job->spids[srank] == cpid)
      break;
  }
  if (i == job->nalive)
    slvpool_abort(1, "I could not find cpid: %d for removal.\n", (int)cpid);

  /* remove it from alivelist */
  i = job->alivemap[srank];
  job->alivemap[srank] = -1;
  job->nalive--;
  if (job->nalive != i) {
    job->alivelist[i] = job->alivelist[job->nalive];
    job->alivemap[job->alivelist[i]] = i;
  }

  slvpool_remove(job, srank);

  BD_LET_LOCK(job->schedule_lock);

  return;
}


/*************************************************************************/
/*! Kills all running slaves
*/
/*************************************************************************/
void slvpool_kill_all(mjob_t *job)
{
  int i, srank;

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] slvpool_kill_all:.\n", job->mynode));

  BD_GET_LOCK(job->schedule_lock);

  for (i=0; i<job->nalive; i++) {
    srank = job->alivelist[i];
    kill(job->spids[srank], 9);
  }

  BD_LET_LOCK(job->schedule_lock);
}


/*************************************************************************/
/*! This function "gracefully" shutsdown the execution environment of
    a bdmpi job. */
/*************************************************************************/
void slvpool_abort(int all, char *f_str,...)
{
  int i, status=0;
  pid_t cpid;
  char errstr[16384];

  /* print the error message */
  va_list argp;
  va_start(argp, f_str);
  vsnprintf(errstr, 16383, f_str, argp);
  va_end(argp);

  fprintf(stderr, "[MSTR%04d]%s", job->mynode, errstr);
  fprintf(stderr, "[MSTR%04d]Aborting the bdmpi job.\n", job->mynode);
  fprintf(stderr, "[MSTR%04d]Killing local running slaves.\n", job->mynode);
  fflush(stderr);
  slvpool_kill_all(job);

  if (all) {
    fprintf(stderr, "[MSTR%04d]Notifying the other masters.\n", job->mynode);
    fflush(stderr);
    for (i=0; i<job->nnodes; i++) {
      if (job->mynode != i) {
        if (MPI_Send(&status, 1, MPI_INT, i, BDMPI_ABORT_TAG, job->mpi_wcomm) != MPI_SUCCESS) {
          fprintf(stderr, "[MSTR%04d]MPI_Send for BDMPI_ABORT_TAG to %d failed.\n",
              job->mynode, i);
          fflush(stderr);
        }
      }
    }
  }


  /* wait for everyone to exit */
  while ((cpid = wait(&status)) != -1) {
    /*
    bdprintf("Child: %d returned with status: %d\n", (int)cpid, status);
    if (WIFEXITED(status)) {
      bdprintf("exited, status=%d\n", WEXITSTATUS(status));
    } else if (WIFSIGNALED(status)) {
      bdprintf("killed by signal %d\n", WTERMSIG(status));
    } else if (WIFSTOPPED(status)) {
      bdprintf("stopped by signal %d\n", WSTOPSIG(status));
    } else if (WIFCONTINUED(status)) {
      bdprintf("continued\n");
    }
    */
  }

  /* cleanup before exiting */
  cleanup_master(job);

  MPI_Finalize();

  exit(EXIT_FAILURE);
}


/*************************************************************************/
/*! Removes a slave from the runnable/running list and adds it into the
    mblocked list. The slave should not be cblocked.
*/
/*************************************************************************/
void slvpool_mblock(mjob_t *job, int rank)
{
  int i;
  static int ts=0;

  BD_GET_LOCK(job->schedule_lock);

  BDASSERT(rank>=0 && rank<job->ns);
  M_IFSET(BDMPI_DBG_IPCM,
      bdprintf("[MSTR%04d] mblocking.0: %d[%+d:%+d]; nalive: %d, nrunning: %d, nrunnable: %d, n[c:m]blocked: %d:%d\n",
              job->mynode, rank, job->mblockedmap[rank], job->cblockedmap[rank],
              job->nalive, job->nrunning, job->nrunnable, job->ncblocked,
              job->nmblocked));

  BDASSERT(job->cblockedmap[rank] == -1);  /* it should not be cblocked */
  BDASSERT(job->mblockedmap[rank] == -1);  /* it should not be mblocked */
  BDASSERT(job->runnablemap[rank] != -1 || job->runningmap[rank] != -1);

  /* move it to blockedlist */
  if (job->mblockedmap[rank] == -1) {
    job->blockedts[rank] = ts++;
    job->mblockedlist[job->nmblocked] = rank;
    job->mblockedmap[rank] = job->nmblocked++;
  }

  /* TODO: Need to put some better invariant checking here. */

  /* remove it from runninglist */
  if (job->runningmap[rank] != -1) {
    i = job->runningmap[rank];
    job->runningmap[rank] = -1;
    job->nrunning--;
    if (job->nrunning != i) {
      job->runninglist[i] = job->runninglist[job->nrunning];
      job->runningmap[job->runninglist[i]] = i;
    }
  }

  /* remove it from runnablelist */
  if (job->runnablemap[rank] != -1) {
    i = job->runnablemap[rank];
    job->runnablemap[rank] = -1;
    job->nrunnable--;
    if (job->nrunnable != i) {
      job->runnablelist[i] = job->runnablelist[job->nrunnable];
      job->runnablemap[job->runnablelist[i]] = i;
    }
  }


  M_IFSET(BDMPI_DBG_IPCM,
      bdprintf("[MSTR%04d] mblocking.1: %d[%+d:%+d]; nalive: %d, nrunning: %d, nrunnable: %d, n[c:m]blocked: %d:%d\n",
              job->mynode, rank, job->mblockedmap[rank], job->cblockedmap[rank],
              job->nalive, job->nrunning, job->nrunnable, job->ncblocked,
              job->nmblocked));

  BD_LET_LOCK(job->schedule_lock);

  return;
}


/*************************************************************************/
/*! Removes a slave from the runnable/running list and adds it into the
    cblocked list. The slave should not be mblocked.
*/
/*************************************************************************/
void slvpool_cblock(mjob_t *job, int rank)
{
  int i;
  static int ts=0;

  BD_GET_LOCK(job->schedule_lock);

  BDASSERT(rank>=0 && rank<job->ns);
  M_IFSET(BDMPI_DBG_IPCM,
      bdprintf("[MSTR%04d] cblocking.0: %d[%+d:%+d]; nalive: %d, nrunning: %d, nrunnable: %d, n[c:b]blocked: %d:%d, \n",
              job->mynode, rank, job->cblockedmap[rank], job->mblockedmap[rank],
              job->nalive, job->nrunning, job->nrunnable, job->ncblocked,
              job->nmblocked));

  BDASSERT(job->cblockedmap[rank] == -1);  /* it should not be cblocked */
  BDASSERT(job->mblockedmap[rank] == -1);  /* it should not be mblocled */
  BDASSERT(job->runnablemap[rank] != -1 || job->runningmap[rank] != -1);

  /* move it to cblockedlist */
  if (job->cblockedmap[rank] == -1) {
    job->blockedts[rank] = ts++;
    job->cblockedlist[job->ncblocked] = rank;
    job->cblockedmap[rank] = job->ncblocked++;
  }

  /* TODO: Need to put some better invariant checking here. */

  /* remove it from runninglist */
  if (job->runningmap[rank] != -1) {
    i = job->runningmap[rank];
    job->runningmap[rank] = -1;
    job->nrunning--;
    if (job->nrunning != i) {
      job->runninglist[i] = job->runninglist[job->nrunning];
      job->runningmap[job->runninglist[i]] = i;
    }
  }

  /* remove it from runnablelist */
  if (job->runnablemap[rank] != -1) {
    i = job->runnablemap[rank];
    job->runnablemap[rank] = -1;
    job->nrunnable--;
    if (job->nrunnable != i) {
      job->runnablelist[i] = job->runnablelist[job->nrunnable];
      job->runnablemap[job->runnablelist[i]] = i;
    }
  }

  M_IFSET(BDMPI_DBG_IPCM,
      bdprintf("[MSTR%04d] cblocking.1: %d[%+d:%+d]; nalive: %d, nrunning: %d, nrunnable: %d, n[c:b]blocked: %d:%d\n",
              job->mynode, rank, job->cblockedmap[rank], job->mblockedmap[rank],
              job->nalive, job->nrunning, job->nrunnable, job->ncblocked,
              job->nmblocked));

  BD_LET_LOCK(job->schedule_lock);

  return;
}


/*************************************************************************/
/*! Removes a slave from the mblocked list and adds it into the runnable
    list as long as it is not in the cblocked list (only if it is blocked).
*/
/*************************************************************************/
void slvpool_munblock(mjob_t *job, int rank)
{
  int i;

  BD_GET_LOCK(job->schedule_lock);

  BDASSERT(rank>=0 && rank<job->ns);
  M_IFSET(BDMPI_DBG_IPCM,
      bdprintf("[MSTR%04d] munblocking.0: %d[%+d, %+d]; nalive: %d, nrunning: %d, nrunnable: %d, n[c:m]blocked: %d:%d\n",
              job->mynode, rank, job->mblockedmap[rank], job->cblockedmap[rank],
              job->nalive, job->nrunning, job->nrunnable, job->ncblocked,
              job->nmblocked));

  /* remove it from mblockedlist */
  if (job->mblockedmap[rank] != -1) {
    i = job->mblockedmap[rank];
    job->mblockedmap[rank] = -1;
    job->nmblocked--;
    if (job->nmblocked != i) {
      job->mblockedlist[i] = job->mblockedlist[job->nmblocked];
      job->mblockedmap[job->mblockedlist[i]] = i;
    }

    /* move it to runnablelist if it is not blocked due to collectives */
    BDASSERT(job->runnablemap[rank] == -1);
    if (job->cblockedmap[rank] == -1 && job->runnablemap[rank] == -1) {
      job->runnablelist[job->nrunnable] = rank;
      job->runnablemap[rank] = job->nrunnable++;
    }
  }

  M_IFSET(BDMPI_DBG_IPCM,
      bdprintf("[MSTR%04d] munblocking.1: %d[%+d, %+d]; nalive: %d, nrunning: %d, nrunnable: %d, n[c:m]blocked: %d:%d\n",
              job->mynode, rank, job->mblockedmap[rank], job->cblockedmap[rank],
              job->nalive, job->nrunning, job->nrunnable, job->ncblocked,
              job->nmblocked));

  BD_LET_LOCK(job->schedule_lock);
}


/*************************************************************************/
/*! Removes a slave from the cblocked list and adds it into the runnable
    list. The slave should be cblocked and should not be mblocked.
*/
/*************************************************************************/
void slvpool_cunblock(mjob_t *job, int rank)
{
  int i;

  BD_GET_LOCK(job->schedule_lock);

  BDASSERT(rank>=0 && rank<job->ns);
  M_IFSET(BDMPI_DBG_IPCM,
      bdprintf("[MSTR%04d] cunblocking.0: %d[%+d, %+d]; nalive: %d, nrunning: %d, nrunnable: %d, n[c:m]blocked: %d:%d\n",
              job->mynode, rank, job->cblockedmap[rank], job->mblockedmap[rank],
              job->nalive, job->nrunning, job->nrunnable, job->ncblocked,
              job->nmblocked));

  BDASSERT(job->cblockedmap[rank] != -1);  /* it should be cblocked */
  BDASSERT(job->mblockedmap[rank] == -1);  /* it should not be mblocked */

  /* remove it from cblockedlist */
  if (job->cblockedmap[rank] != -1) {
    i = job->cblockedmap[rank];
    job->cblockedmap[rank] = -1;
    job->ncblocked--;
    if (job->ncblocked != i) {
      job->cblockedlist[i] = job->cblockedlist[job->ncblocked];
      job->cblockedmap[job->cblockedlist[i]] = i;
    }

    /* move it to runnablelist if it is not blocked due to collectives */
    BDASSERT(job->runnablemap[rank] == -1); /* it should not be in runnable */
    if (job->runnablemap[rank] == -1) {
      job->runnablelist[job->nrunnable] = rank;
      job->runnablemap[rank] = job->nrunnable++;
    }
  }

  M_IFSET(BDMPI_DBG_IPCM,
      bdprintf("[MSTR%04d] cunblocking.1: %d[%+d, %+d]; nalive: %d, nrunning: %d, nrunnable: %d, n[c:m]blocked: %d:%d\n",
              job->mynode, rank, job->cblockedmap[rank], job->mblockedmap[rank],
              job->nalive, job->nrunning, job->nrunnable, job->ncblocked,
              job->nmblocked));

  BD_LET_LOCK(job->schedule_lock);
}
