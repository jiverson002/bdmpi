/*
 * Copyright 2013, Regents of the University of Minnesota
 *
 * bdmprun.c
 *
 * The entry point of the jobrun program.
 *
 * Started 3/31/2013
 * George
 *
 * $Id: gpmetis.c 13900 2013-03-24 15:27:07Z karypis $
 *
 */

#define BDMPRUN_MAIN
#include "bdmprun.h"


/*************************************************************************/
/*! Let the game begin! */
/*************************************************************************/
int main(int argc, char *argv[])
{
  int status=0, provided;
  pid_t cpid;

  mallopt(M_TRIM_THRESHOLD, 64*4096);
  mallopt(M_MMAP_THRESHOLD, 64*4096);

  setbuf(stdout, NULL);
  setbuf(stderr, NULL);

  BDASSERT(MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided)
           == MPI_SUCCESS);
  if (provided != MPI_THREAD_MULTIPLE)
    errexit("Failed to get the required thread support: provided: %d\n", provided);

  job = parse_cmdline(argc, argv);

  /* start the slaves */
  spawn_slaves(job);

  /* start listening to requests until all slaves have finished */
  slvpool_listen(job);

  /* wait for everyone to exit [this should not be needed at this point] */
  while ((cpid = wait(&status)) != -1) {
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
  }

  /* cleanup before exiting */
  cleanup_master(job);

  MPI_Finalize();

  exit(EXIT_SUCCESS);
}


/*************************************************************************/
/*! Forks the slaves processes to perform the exec() */
/*************************************************************************/
void spawn_slaves(mjob_t *job)
{
  int i;
  pid_t cpid;

  /* pre-forking setup steps */
  setup_master_prefork(job);

  /* lock the global shared memory, as you will need it later to set
     up the pid-to-rank mappings during BDMPI_Init() */
  bdsm_lock(job->globalSM);

  /* print the command-line arguments */
  if (job->dbglvl>0 && job->mynode == 0) {
    for (i=0; job->exeargv[i]!=NULL; i++)
      bdprintf("exeargv[%d]=%s\n", i, job->exeargv[i]);
  }


  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR] Ready to spawn\n"));

  /* spawn the slave processes */
  for (i=0; i<job->ns; i++) {
    switch (cpid = fork()) {
      case -1:
        errexit("Failed on fork(): %s\n", strerror(errno));
        break;

      case 0:
        execvp(job->exefile, job->exeargv);
        errexit("Failed on execvp(): %s\n", strerror(errno));
        break;

      default:
        job->spids[i] = cpid;
        M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR] Child process: %d\n", (int)cpid));
        break;
    }
  }

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR] Done spawning\n"));

  /* post-forking setup steps */
  setup_master_postfork(job);

  /* setup communicators */
  comm_setup(job);

  /* unlock the global shared memory */
  bdsm_unlock(job->globalSM);

  return;
}


/*************************************************************************/
/*! Sets up various master process info prior to fork() */
/*************************************************************************/
void setup_master_prefork(mjob_t *job)
{
  int p;
  char *wdir, *tmpstr;
  pthread_mutexattr_t mtx_attr;

  /* initialize timers */
  gk_clearwctimer(job->totalTmr);
  gk_clearwctimer(job->routeTmr);
  gk_clearwctimer(job->sendTmr);
  gk_clearwctimer(job->recvTmr);
  gk_clearwctimer(job->colTmr);
  gk_clearwctimer(job->barrierTmr);
  gk_clearwctimer(job->commTmr);
  gk_clearwctimer(job->aux1Tmr);
  gk_clearwctimer(job->aux2Tmr);
  gk_clearwctimer(job->aux3Tmr);

  /* start the overall timer */
  gk_startwctimer(job->totalTmr);

  /* populate the various fields for the job */
  job->mpid    = getpid();
  job->imsize *= sysconf(_SC_PAGESIZE);
  job->mmsize *= sysconf(_SC_PAGESIZE);
  job->smsize *= sysconf(_SC_PAGESIZE);

  /* setup MPI-related information */
  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR] Setting up MPI environment\n"));
  MPI_Comm_dup(MPI_COMM_WORLD, &(job->mpi_wcomm));
  MPI_Comm_size(job->mpi_wcomm, &job->nnodes);
  MPI_Comm_rank(job->mpi_wcomm, &job->mynode);

  job->slvdist = gk_imalloc(job->nnodes+1, "setup_master_prefork: slvdist");
  MPI_Allgather(&job->ns, 1, MPI_INT, job->slvdist, 1, MPI_INT, job->mpi_wcomm);
  MAKECSR(p, job->nnodes, job->slvdist);

  /* create the global shared memory */
  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR] Creating globalSM\n"));
  job->globalSM = bdsm_create(BDMPI_GLOBALSMSIZE, "global", (int)job->mpid);

  /* create the MQ for the slaves to send requests to the master */
  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR] Creating reqMQ for masterpid: %d\n", (int)job->mpid));
  job->reqMQ = bdmq_create("reqMQ", (int)job->mpid);

  /* create the mutexes */
  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR] Creating mutexes for masterpid: %d\n", (int)job->mpid));
  job->mlockMX    = bdlock_create("mlockMX", (int)job->mpid, 1);
  job->criticalMX = bdlock_create("criticalMX", (int)job->mpid, job->nc);

  /* allocate memory for jdesc and setup whatever info you can setup now. */
  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR] Allocating job->jdesc\n"));
  job->jdesc = (bdjdesc_t *)bdsm_malloc(job->globalSM, sizeof(bdjdesc_t), "job->jdesc");
  memset(job->jdesc, 0, sizeof(bdjdesc_t));

  job->jdesc->nnodes  = job->nnodes;
  job->jdesc->mynode  = job->mynode;
  job->jdesc->np      = job->slvdist[job->nnodes];
  job->jdesc->soffset = job->slvdist[job->mynode];
  job->jdesc->ns      = job->ns;
  job->jdesc->nr      = job->nr_input;
  job->jdesc->smsize  = job->smsize;
  job->jdesc->imsize  = job->imsize;
  job->jdesc->sbsize  = job->sbsize;
  job->jdesc->dbglvl  = job->dbglvl;
  job->jdesc->mpid    = job->mpid;

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR] nnodes: %d, mynode: %d, np: %d, ns: %d\n",
        job->nnodes, job->mynode, job->slvdist[job->nnodes], job->ns));

  /* allocate memory for spids in the global SM */
  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR] Allocating job->spids\n"));
  job->spids = (pid_t *)bdsm_malloc(job->globalSM, sizeof(pid_t)*job->ns, "job->spids");
  memset(job->spids, 0, sizeof(pid_t)*job->ns);


  /* setup the working directory */
  snprintf(job->jdesc->wdir, BDMPI_WDIR_LEN, "%s/%d", job->iwdir, (int)job->mpid);
  gk_free((void **)&job->iwdir, LTERM);

  if (gk_mkpath(job->jdesc->wdir) == -1)
    errexit("Failed to create the wdir: %s\n", job->jdesc->wdir);
  if (job->dbglvl>0)
    printf("Creating working directory: %s\n", job->jdesc->wdir);

  xfer_setwdir(job->jdesc->wdir);

#ifdef XXX
  /* setup the message storing directories */
  tmpstr = gk_cmalloc(strlen(job->jdesc->wdir)+100, "setup_master: tmpstr");
  sprintf(tmpstr, "%s/col", job->jdesc->wdir);
  if (gk_mkpath(tmpstr) == -1)
    errexit("Failed to create the directory: %s\n", tmpstr);

  sprintf(tmpstr, "%s/p2p", job->jdesc->wdir);
  if (gk_mkpath(tmpstr) == -1)
    errexit("Failed to create the directory: %s\n", tmpstr);

  /*  No per-pe directories at this time
  for (p=0; p<job->ns; p++) {
    sprintf(tmpstr, "%s/col/%d", job->jdesc->wdir, p);
    if (gk_mkpath(tmpstr) == -1)
      errexit("Failed to create the directory: %s\n", tmpstr);

    sprintf(tmpstr, "%s/p2p/%d", job->jdesc->wdir, p);
    if (gk_mkpath(tmpstr) == -1)
      errexit("Failed to create the directory: %s\n", tmpstr);
  }
  */
#endif

#ifdef XXX
  /* setup the msgdb */
  job->nextmsgid = 0;
  sprintf(tmpstr, "%s/msgdb", job->jdesc->wdir);
  job->msgdbfname = tmpstr;
  if ((job->msgdb = db_create(job->msgdbfname)) == NULL)
    errexit("Failed to create the msgdb: %s\n", job->msgdbfname);
#endif

  /* setup the pending requests */
  pending_setup(job);

  /* initialize pthreads */
  job->schedule_lock = (pthread_mutex_t *)gk_malloc(sizeof(pthread_mutex_t), "schedule_lock");
  job->comm_lock = (pthread_mutex_t *)gk_malloc(sizeof(pthread_mutex_t), "comm_lock");

  BDASSERT(pthread_mutexattr_init(&mtx_attr) == 0);
  BDASSERT(pthread_mutexattr_settype(&mtx_attr, PTHREAD_MUTEX_RECURSIVE) == 0);

  BDASSERT(pthread_mutex_init(job->schedule_lock, &mtx_attr) == 0);
  BDASSERT(pthread_mutex_init(job->comm_lock, &mtx_attr) == 0);
  BDASSERT(pthread_mutexattr_destroy(&mtx_attr) == 0);

  return;
}


/*************************************************************************/
/*! Sets up various master process info post to fork() */
/*************************************************************************/
void setup_master_postfork(mjob_t *job)
{
  int i;

  /* allocate the per-slave shared memory */
  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR] Allocating SCBs\n"));
  job->scbs = (bdscb_t **)gk_malloc(sizeof(bdscb_t *)*job->ns, "job->scbs");
  for (i=0; i<job->ns; i++)
    job->scbs[i] = bdscb_create(job->smsize, "slave", (int)job->spids[i]);

  /* create the per-slave message queues */
  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR] Creating per-slave message queues\n"));
  job->goMQs = (bdmq_t **)gk_malloc(sizeof(bdmq_t *)*job->ns, "job->goMQs");
  for (i=0; i<job->ns; i++)
    job->goMQs[i] = bdmq_create("goMQ", (int)job->spids[i]);

  /* create the per-slave message queues for communication to/from master */
  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR] Creating the 1-1 queues\n"));
  job->c2sMQs = (bdmq_t **)gk_malloc(sizeof(bdmq_t)*job->ns, "job->c2sMQs");
  job->c2mMQs = (bdmq_t **)gk_malloc(sizeof(bdmq_t)*job->ns, "job->c2mMQs");
  for (i=0; i<job->ns; i++) {
    job->c2sMQs[i] = bdmq_create("c2sMQ", (int)job->spids[i]);
    job->c2mMQs[i] = bdmq_create("c2mMQ", (int)job->spids[i]);
  }

  /* allocate the various lists/maps to store the state of the slaves */
  job->alivelist    = gk_imalloc(job->ns, "alivelist");
  job->runnablelist = gk_imalloc(job->ns, "runnablelist");
  job->runninglist  = gk_imalloc(job->ns, "runninglist");
  job->mblockedlist = gk_imalloc(job->ns, "mblockedlist");
  job->cblockedlist = gk_imalloc(job->ns, "cblockedlist");
  job->alivemap     = gk_ismalloc(job->ns, -1, "alivemap");
  job->runnablemap  = gk_ismalloc(job->ns, -1, "runnablemap");
  job->runningmap   = gk_ismalloc(job->ns, -1, "runningmap");
  job->mblockedmap  = gk_ismalloc(job->ns, -1, "mblockedmap");
  job->cblockedmap  = gk_ismalloc(job->ns, -1, "cblockedmap");
  job->blockedts    = gk_ismalloc(job->ns, -1, "blockedts");

  /* everybody is alive and running (i.e., no co-operating scheduling) and
     runnablelist and blockedlist are empty. */
  job->nalive    = job->ns;
  job->nrunning  = job->ns;
  job->nrunnable = 0;
  job->nmblocked = 0;
  job->ncblocked = 0;
  for (i=0; i<job->ns; i++) {
    job->alivelist[i] = i;
    job->alivemap[i]  = i;

    job->runninglist[i] = i;
    job->runningmap[i]  = i;
  }

  job->njoined = 0;  /* controls the init/finalize barrier */

  job->nR = job->ns; /* initially, no co-operating scheduling */

  /* populate the various memory statistics */
  if (job->rmsize > 63)
    errexit("Invalid rmsize: %d\n", job->rmsize);
  job->memrss = 0;
  job->memmax = 1LLU<<job->rmsize;
  job->slvrss = (size_t *)gk_malloc(job->ns*sizeof(size_t), "slvrss");
  job->slvtot = (size_t *)gk_malloc(job->ns*sizeof(size_t), "slvtot");
  memset(job->slvrss, 0, job->ns*sizeof(size_t));
  memset(job->slvtot, 0, job->ns*sizeof(size_t));

  return;
}


/*************************************************************************/
/*! Cleans up various items created by the job */
/*************************************************************************/
void cleanup_master(mjob_t *job)
{
  int i;
  double maxtmr;

  gk_stopwctimer(job->totalTmr);

  sleep(1);
  bdprintf("------------------------------------------------\n");
  bdprintf("Master %d is done.\n", job->mynode);
  bdprintf("Memory stats [%zu / %zu]\n", job->memrss, job->memmax);

  gk_rmpath(job->jdesc->wdir);

  /* clean up the various per-slave message queues and shared memory regions */
  for (i=0; i<job->ns; i++) {
    bdprintf("       [%3d] [%zu / %zu]\n", i, job->slvrss[i], job->slvtot[i]);
    bdscb_destroy(job->scbs[i]);
    bdmq_destroy(job->goMQs[i]);
    bdmq_destroy(job->c2sMQs[i]);
    bdmq_destroy(job->c2mMQs[i]);
  }

  bdmq_destroy(job->reqMQ);
  bdsm_destroy(job->globalSM);

  bdlock_destroy(job->mlockMX);
  bdlock_destroy(job->criticalMX);

  comm_cleanup(job);

  pending_cleanup(job);

  /* print timing info */
  if (job->mynode == 0) {
    bdprintf("------------------------------------------------\n");
    bdprintf("Master timings\n");
  }

  MPI_Reduce(&job->totalTmr, &maxtmr, 1, MPI_DOUBLE, MPI_MAX, 0, job->mpi_wcomm);
  if (job->mynode == 0)
    bdprintf("    totalTmr:   %8.3lfs\n", maxtmr);

  MPI_Reduce(&job->routeTmr, &maxtmr, 1, MPI_DOUBLE, MPI_MAX, 0, job->mpi_wcomm);
  if (job->mynode == 0)
    bdprintf("    routeTmr:   %8.3lfs\n", maxtmr);

  MPI_Reduce(&job->sendTmr, &maxtmr, 1, MPI_DOUBLE, MPI_MAX, 0, job->mpi_wcomm);
  if (job->mynode == 0 && maxtmr>0)
    bdprintf("     sendTmr:   %8.3lfs\n", maxtmr);

  MPI_Reduce(&job->recvTmr, &maxtmr, 1, MPI_DOUBLE, MPI_MAX, 0, job->mpi_wcomm);
  if (job->mynode == 0 && maxtmr>0)
    bdprintf("     recvTmr:   %8.3lfs\n", maxtmr);

  MPI_Reduce(&job->colTmr, &maxtmr, 1, MPI_DOUBLE, MPI_MAX, 0, job->mpi_wcomm);
  if (job->mynode == 0 && maxtmr>0)
    bdprintf("      colTmr:   %8.3lfs\n", maxtmr);

  MPI_Reduce(&job->barrierTmr, &maxtmr, 1, MPI_DOUBLE, MPI_MAX, 0, job->mpi_wcomm);
  if (job->mynode == 0 && maxtmr>0)
    bdprintf("  barrierTmr:   %8.3lfs\n", maxtmr);

  MPI_Reduce(&job->commTmr, &maxtmr, 1, MPI_DOUBLE, MPI_MAX, 0, job->mpi_wcomm);
  if (job->mynode == 0 && maxtmr>0)
    bdprintf("     commTmr:   %8.3lfs\n", maxtmr);

  MPI_Reduce(&job->aux1Tmr, &maxtmr, 1, MPI_DOUBLE, MPI_MAX, 0, job->mpi_wcomm);
  if (job->mynode == 0 && maxtmr>0)
    bdprintf("     aux1Tmr:   %8.3lfs\n", maxtmr);

  MPI_Reduce(&job->aux2Tmr, &maxtmr, 1, MPI_DOUBLE, MPI_MAX, 0, job->mpi_wcomm);
  if (job->mynode == 0 && maxtmr>0)
    bdprintf("     aux2Tmr:   %8.3lfs\n", maxtmr);

  MPI_Reduce(&job->aux3Tmr, &maxtmr, 1, MPI_DOUBLE, MPI_MAX, 0, job->mpi_wcomm);
  if (job->mynode == 0 && maxtmr>0)
    bdprintf("     aux3Tmr:   %8.3lfs\n", maxtmr);


  if (job->mynode == 0)
    bdprintf("------------------------------------------------\n");

  pthread_mutex_destroy(job->schedule_lock);
  pthread_mutex_destroy(job->comm_lock);

  gk_free((void **)&job->scbs, &job->goMQs, &job->c2sMQs, &job->c2mMQs,
      &job->alivelist, &job->runnablelist, &job->runninglist,
      &job->mblockedlist, &job->cblockedlist,
      &job->alivemap, &job->runnablemap, &job->runningmap,
      &job->mblockedmap, &job->cblockedmap,
      &job->blockedts,
      &job->slvdist,
      &job->slvrss, &job->slvtot,
      &job->schedule_lock, &job->comm_lock,
      &job, LTERM);

  return;
}
