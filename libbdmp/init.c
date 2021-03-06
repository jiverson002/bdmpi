/*
 * init.c
 *
 * Implements BDMPI_Init()
 *
 * Started 4/3/2013
 * George
 *
 */

#define BDMPLIB_INIT_C


#include "bdmplib.h"


/*************************************************************************/
/* Initializes the BDMP library. */
/*************************************************************************/
int bdmp_Init(sjob_t **r_job, int *argc, char **argv[])
{
  int i, opts;
  pid_t mpid;
  sjob_t *job;
  bdmsg_t donemsg, gomsg;

  mallopt(M_TRIM_THRESHOLD, 64*4096);
  mallopt(M_MMAP_THRESHOLD, 64*4096);

  job = *r_job = (sjob_t *)bd_malloc(sizeof(sjob_t), "bdmp_Init: job");
  memset(job, 0, sizeof(sjob_t));

  mpid       = getppid();
  job->mypid = getpid();

  /* open the global SMR */
  job->globalSM = bdsm_open(BDMPI_GLOBALSMSIZE, "global", mpid);

  /* lock and then right away unlock the global SMR in order to
     ensure that the master has finished setting it up. */
  /* need to use the libc variations of these, since sbma is not initialized
   * yet. */
  bdsm_lock(job->globalSM);
  bdsm_unlock(job->globalSM);

  /* hook job->jdesc into the global shared memory */
  job->jdesc = (bdjdesc_t *)bdsm_malloc(job->globalSM, sizeof(bdjdesc_t), "job->jdesc");

  /* hook job->spids into the global shared memory */
  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Init: Hooking into global spids\n"));

  job->spids = (pid_t *)bdsm_malloc(job->globalSM, sizeof(pid_t)*job->jdesc->ns, "job->spids");

  /* hook job->mallinfo into the global shared memory */
  job->mallinfo = (struct mallinfo*)bdsm_malloc(job->globalSM,
    sizeof(struct mallinfo)*job->jdesc->ns, "job->mallinfo");
  /* hook job->timeinfo into the global shared memory */
  job->timeinfo = (struct sbma_timeinfo*)bdsm_malloc(job->globalSM,
    sizeof(struct sbma_timeinfo)*job->jdesc->ns, "job->timeinfo");

  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Init: Info: ns: %d, mpids: %d:%d\n",
        job->jdesc->ns, (int)job->jdesc->mpid, (int)mpid));

  /* do some sanity checks */
  if (mpid != job->jdesc->mpid) {
    bdprintf("The master pid and the parent pid did not match! %d %d\n",
        (int)mpid, (int)job->jdesc->mpid);
    exit(EXIT_SUCCESS);
  }

  /* determine intra rank */
  for (i=0; i<job->jdesc->ns; i++) {
    S_IFSET(BDMPI_DBG_IPCS, bdprintf("spids[%d]=%d\n", i, (int)job->spids[i]));
    if (job->mypid == job->spids[i]) {
      job->lrank = i;
      job->rank  = job->jdesc->soffset + i;
      break;
    }
  }
  if (i == job->jdesc->ns) {
    bdprintf("Mypid is not in the job->spids! [%d]\n", (int)job->mypid);
    exit(EXIT_SUCCESS);
  }

  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Init: Info: np: %d, ns: %d, rank: %d, lrank: %d\n",
        job->jdesc->np, job->jdesc->ns, job->rank, job->lrank));

  /* open the message queues that the master has created */
  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Init: Opening MQs\n"));

  job->reqMQ = bdmq_open("reqMQ", (int)mpid);
  job->c2sMQ = bdmq_open("c2sMQ", (int)job->mypid);
  job->c2mMQ = bdmq_open("c2mMQ", (int)job->mypid);
  job->goMQ  = bdmq_open("goMQ", (int)job->mypid);

  /* open the mutexes that the master has created */
  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Init: Opening MXs\n"));

  job->mlockMX    = bdlock_open("mlockMX", (int)mpid);
  job->criticalMX = bdlock_open("criticalMX", (int)mpid);


  /* open the SCM for master-slave communication */
  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Init: Opening SCB\n"));

  job->scb = bdscb_open(job->jdesc->smsize, "slave", (int)job->mypid);

  /* setup communicators */
  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Init: Setting up communicators\n"));

  /* BDMPI_COMM_WORLD */
  BDMPI_COMM_WORLD = (bdscomm_t *)bd_malloc(sizeof(bdscomm_t), "BDMPI_COMM_WORLD");
  BDMPI_COMM_WORLD->mcomm = 0;  /* This is hard-coded and is OK */
  BDMPI_COMM_WORLD->size  = job->jdesc->np;
  BDMPI_COMM_WORLD->rank  = job->rank;
  BDMPI_COMM_WORLD->lsize = job->jdesc->ns;
  BDMPI_COMM_WORLD->lrank = job->lrank;
  BDMPI_COMM_WORLD->nsize = job->jdesc->nnodes;
  BDMPI_COMM_WORLD->nrank = job->jdesc->mynode;
  BDMPI_COMM_WORLD->rrank = job->jdesc->soffset;
  BDMPI_COMM_WORLD->copid = 1;

  /* BDMPI_COMM_NODE */
  BDMPI_COMM_NODE = (bdscomm_t *)bd_malloc(sizeof(bdscomm_t), "BDMPI_COMM_NODE");
  BDMPI_COMM_NODE->mcomm = 1;  /* This is hard-coded and is OK */
  BDMPI_COMM_NODE->size  = job->jdesc->ns;
  BDMPI_COMM_NODE->rank  = job->lrank;
  BDMPI_COMM_NODE->lsize = job->jdesc->ns;
  BDMPI_COMM_NODE->lrank = job->lrank;
  BDMPI_COMM_NODE->nsize = 1;
  BDMPI_COMM_NODE->nrank = 0;
  BDMPI_COMM_NODE->rrank = 0;
  BDMPI_COMM_NODE->copid = 1;


  /* other.... */
  job->smallmsg = gk_max(sizeof(bdmsg_t)*(job->jdesc->ns+3), job->jdesc->imsize);
  //job->smallmsg = 0;

  xfer_setwdir(job->jdesc->wdir);

  /* tell the master that you are done with init */
  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Init: Notify the master that you are done with init\n"));

  memset(&donemsg, 0, sizeof(bdmsg_t));
  donemsg.msgtype = BDMPI_MSGTYPE_INIT;
  donemsg.myrank  = job->rank;

  if (bdmq_send(job->reqMQ, &donemsg, sizeof(bdmsg_t)) == -1)
    bdprintf("Failed on sending a donemsg: %s.\n", strerror(errno));

  /* wait for a message from the master to go */
  S_IFSET(BDMPI_DBG_IPCS,
      bdprintf("BDMPI_Init: Waiting for a go message [goMQlen: %d]\n", bdmq_length(job->goMQ)));

  /* sleep... */
  BDMPL_SLEEP(job, gomsg, 0);

  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Init: I got the following gomsg: "
    "%d\n", gomsg.msgtype));

  /* ====================================================================== */
  /* everything above here must have been allocated via the libc interface. */
  /* ====================================================================== */

  if (-1 == SBMA_init(job->jdesc->wdir, job->jdesc->mpid,\
    job->jdesc->pgsize*sysconf(_SC_PAGESIZE), job->jdesc->ns,\
    job->jdesc->rmsize, SBMA_parse_optstr(job->jdesc->sboptstr)))
  {
    bdprintf("Failed to init sbma\n");
  }

  /* Need a barrier here to ensure that SBMA is properly initialized */
  BDMPI_Barrier(BDMPI_COMM_WORLD);

  /* create additional standard communicators -- must come after memory
   * management environments are created */
  BDASSERT(BDMPI_Comm_split(BDMPI_COMM_WORLD, BDMPI_COMM_WORLD->rank, 1,\
    &BDMPI_COMM_SELF) == BDMPI_SUCCESS);
  BDASSERT(BDMPI_Comm_split(BDMPI_COMM_WORLD, 1, BDMPI_COMM_NODE->rank,\
    &BDMPI_COMM_CWORLD) == BDMPI_SUCCESS);

  return BDMPI_SUCCESS;
}
