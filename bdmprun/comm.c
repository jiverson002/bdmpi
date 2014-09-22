/*
 * Copyright 2013, Regents of the University of Minnesota
 *
 * comm.c
 *
 * Various functions for dealing with communicators.
 *
 * Started 4/4/2013
 * George
 *
 */

#include "bdmprun.h"



/*************************************************************************/
/*! Setup the communicator information and create BDMPI_COMM_WORLD */
/*************************************************************************/
void comm_setup(mjob_t *job)
{
  int i, j;
  bdmcomm_t *comm;

  /* this will be used to determine the unique ID of a communicator across
     the nodes that are involved in the communicator */
  job->next_mpi_commid = 10;

  job->next_mpi_tag = 1000;  /* used to uniquely tag the send's out of a node */ 

  job->maxncomm   = BDMPI_INIT_MAXNCOMM;
  job->nfreecomm  = job->maxncomm - 2; /* to account for the WORLD/NODE communicators */

  /* create and populate the free list of communicators */
  job->commfreelist = gk_imalloc(job->maxncomm, "comm_init: commfreelist");
  for (i=0; i<job->nfreecomm; i++)
    job->commfreelist[i] = i+1;

  /* create the array of communicator head pointers */
  job->comms = (bdmcomm_t **)gk_malloc(sizeof(bdmcomm_t *)*job->maxncomm, "comm_init: comms");
  for (i=0; i<job->maxncomm; i++)
    job->comms[i] = NULL;


  /* --------------------------------------- */
  /* create the BDMPI_COMM_WORLD communicator */
  /* --------------------------------------- */
  comm = job->comms[BDMPI_COMM_WORLD] = (bdmcomm_t *)gk_malloc(sizeof(bdmcomm_t), "comms[WORLD]");
  memset(comm, 0, sizeof(bdmcomm_t));

  /* intra-node component */
  comm->mcomm   = BDMPI_COMM_WORLD;
  comm->lsize   = job->ns;
  comm->counter = job->ns;
  comm->copid   = 0;
  comm->sranks  = gk_imalloc(comm->lsize, "comm->sranks");
  for (i=0; i<comm->lsize; i++)
    comm->sranks[i] = i;

  /* inter-node component */
  comm->mpi_commid = BDMPI_COMM_WORLD;

  MPI_Comm_dup(job->mpi_wcomm, &(comm->mpi_comm));
  MPI_Comm_size(comm->mpi_comm, &(comm->nnodes));
  MPI_Comm_rank(comm->mpi_comm, &(comm->mynode));

  comm->wnranks = gk_imalloc(comm->nnodes, "wnranks");
  MPI_Allgather(&(job->mynode), 1, MPI_INT, comm->wnranks, 1, MPI_INT, comm->mpi_comm);

  /* inter<->intra coupling */
  comm->slvdist = gk_imalloc(comm->nnodes+1, "comm->slvdist");
  MPI_Allgather(&(comm->lsize), 1, MPI_INT, comm->slvdist, 1, MPI_INT, comm->mpi_comm);
  MAKECSR(i, comm->nnodes, comm->slvdist);

  comm->gsize = comm->slvdist[comm->nnodes];

  comm->g2nmap = gk_imalloc(comm->gsize, "comm->g2nmap");
  comm->g2lmap = gk_imalloc(comm->gsize, "comm->g2lmap");
  comm->l2gmap = gk_imalloc(comm->lsize, "comm->l2gmap");

  for (i=0; i<comm->nnodes; i++) {
    for (j=comm->slvdist[i]; j<comm->slvdist[i+1]; j++) {
      comm->g2nmap[j] = i;
      comm->g2lmap[j] = j-comm->slvdist[i];
    }
  }
  for (i=0; i<comm->lsize; i++)
    comm->l2gmap[i] = comm->slvdist[comm->mynode]+i;

  /* auxiliary */
  comm->rwlock = (pthread_rwlock_t *)gk_malloc(sizeof(pthread_rwlock_t), "comm->lock");
  BDASSERT(pthread_rwlock_init(comm->rwlock, NULL) == 0);


  /* -------------------------------------- */
  /* create the BDMPI_COMM_NODE communicator */
  /* -------------------------------------- */
  comm = job->comms[BDMPI_COMM_NODE] = (bdmcomm_t *)gk_malloc(sizeof(bdmcomm_t), "comms[WORLD]");
  memset(comm, 0, sizeof(bdmcomm_t));

  /* intra-node component */
  comm->mcomm   = BDMPI_COMM_NODE;
  comm->lsize   = job->ns;
  comm->counter = job->ns;
  comm->copid   = 0;
  comm->sranks  = gk_imalloc(comm->lsize, "comm->sranks");
  for (i=0; i<comm->lsize; i++)
    comm->sranks[i] = i;

  /* inter-node component */
  comm->mpi_commid = BDMPI_COMM_NODE;

  MPI_Comm_dup(MPI_COMM_SELF, &(comm->mpi_comm));
  MPI_Comm_size(comm->mpi_comm, &(comm->nnodes));
  MPI_Comm_rank(comm->mpi_comm, &(comm->mynode));

  comm->wnranks = gk_imalloc(comm->nnodes, "wnranks");
  MPI_Allgather(&(job->mynode), 1, MPI_INT, comm->wnranks, 1, MPI_INT, comm->mpi_comm);

  /* inter<->intra coupling */
  comm->slvdist = gk_imalloc(comm->nnodes+1, "comm->slvdist");
  MPI_Allgather(&(comm->lsize), 1, MPI_INT, comm->slvdist, 1, MPI_INT, comm->mpi_comm);
  MAKECSR(i, comm->nnodes, comm->slvdist);

  comm->gsize = comm->slvdist[comm->nnodes];

  comm->g2nmap = gk_imalloc(comm->gsize, "comm->g2nmap");
  comm->g2lmap = gk_imalloc(comm->gsize, "comm->g2lmap");
  comm->l2gmap = gk_imalloc(comm->lsize, "comm->l2gmap");

  for (i=0; i<comm->nnodes; i++) {
    for (j=comm->slvdist[i]; j<comm->slvdist[i+1]; j++) {
      comm->g2nmap[j] = i;
      comm->g2lmap[j] = j-comm->slvdist[i];
    }
  }
  for (i=0; i<comm->lsize; i++)
    comm->l2gmap[i] = comm->slvdist[comm->mynode]+i;

  /* auxiliary */
  comm->rwlock = (pthread_rwlock_t *)gk_malloc(sizeof(pthread_rwlock_t), "comm->lock");
  BDASSERT(pthread_rwlock_init(comm->rwlock, NULL) == 0);


  return;
}


/*************************************************************************/
/*! Cleanup the communicator information prior to exit */
/*************************************************************************/
void comm_cleanup(mjob_t *job)
{
  int i;

  BD_GET_LOCK(job->comm_lock);

  for (i=0; i<job->maxncomm; i++) {
    if (job->comms[i] != NULL)
      comm_free(job, i);
  }

  gk_free((void **)&job->commfreelist, &job->comms, LTERM);

  BD_LET_LOCK(job->comm_lock);

  return;
}


/*************************************************************************/
/*! Frees a communicator */
/*************************************************************************/
void comm_free(mjob_t *job, int comm)
{

  BD_GET_LOCK(job->comm_lock);

  if (job->maxncomm < comm)
    slvpool_abort(1, "Trying to free a non-existing communicator [%d %d]!\n", 
        comm, job->maxncomm);

  if (job->comms[comm] == NULL)
    slvpool_abort(1, "Trying to free a non-existing communicator!\n");

  /* we assume that this is a local operation. If not, it may deadlock! */
  MPI_Comm_free(&(job->comms[comm]->mpi_comm));

  pthread_rwlock_destroy(job->comms[comm]->rwlock);
  
  gk_free((void **)&(job->comms[comm]->sranks), 
          &(job->comms[comm]->slvdist),
          &(job->comms[comm]->g2nmap),
          &(job->comms[comm]->g2lmap),
          &(job->comms[comm]->l2gmap),
          &(job->comms[comm]->wnranks),
          &(job->comms[comm]->rwlock),
          &(job->comms[comm]), 
          LTERM);

  job->comms[comm] = NULL;  /* this is already done in gk_free, but just in case :) */
  job->commfreelist[job->nfreecomm++] = comm;

  BD_LET_LOCK(job->comm_lock);

  return;
}


/*************************************************************************/
/* Reallocates more memory for the communicators */
/*************************************************************************/
void comm_realloc_master(mjob_t *job)
{
  size_t i;

  BD_GET_LOCK(job->comm_lock);

  job->comms = (bdmcomm_t **)gk_realloc(job->comms, sizeof(bdmcomm_t *)*2*job->maxncomm, "comm_realloc: comms");
  for (i=job->maxncomm; i<2*job->maxncomm; i++)
    job->comms[i] = NULL;

  job->commfreelist = gk_irealloc(job->commfreelist, 2*job->maxncomm, "comm_realloc: commfreelist");
  for (i=0; i<job->maxncomm; i++)
    job->commfreelist[i] = job->maxncomm+i;

  job->nfreecomm = job->maxncomm;
  job->maxncomm *= 2;

  BD_LET_LOCK(job->comm_lock);

  return;
}


/*************************************************************************/
/*! Response to a BDMPI_Comm_dup.
    Protocol:
       Blocks the process.
       Decreases the counter associated with the counter.
       If counter becomes 0, then a new master communicator is created 
       Sets the counter back to the size of that communicator.
       Copies the per-slave communicator info in the per-slave scbs
*/
/*************************************************************************/
void *mstr_comm_dup(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  size_t i;
  int srank, oldcomm, newcomm, mpi_commid;
  bdscomm_t scomm;

  /* block the slave */
  slvpool_cblock(job, babel_get_srank(job->comms[msg->mcomm], msg->myrank));

  BD_GET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  M_IFSET(BDMPI_DBG_IPCM, 
      bdprintf("[MSTR%04d] comm_dup: mcomm: %d, counter: %d [entering]\n", 
        job->mynode, msg->mcomm, job->comms[msg->mcomm]->counter));

  /* see if all local slaves called this collective */
  if (--job->comms[msg->mcomm]->counter == 0) {
    BD_LET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

    oldcomm = msg->mcomm;

    BD_GET_LOCK(job->comm_lock);
    if (job->nfreecomm == 0)
      comm_realloc_master(job);

    newcomm = job->commfreelist[--job->nfreecomm];
    job->comms[newcomm] = (bdmcomm_t *)gk_malloc(sizeof(bdmcomm_t), "comms[newcomm]");
    memset(job->comms[newcomm], 0, sizeof(bdmcomm_t));
    BD_LET_LOCK(job->comm_lock);

    /* intra-node component */
    job->comms[newcomm]->mcomm   = newcomm;
    job->comms[newcomm]->lsize   = job->comms[oldcomm]->lsize;
    job->comms[newcomm]->counter = job->comms[newcomm]->lsize;
    job->comms[newcomm]->copid   = 0;
    job->comms[newcomm]->sranks  = gk_imalloc(job->comms[newcomm]->lsize, "comms[newcomm]->sranks");
    for (i=0; i<job->comms[newcomm]->lsize; i++)
      job->comms[newcomm]->sranks[i] = job->comms[oldcomm]->sranks[i];

    /* inter-node component */
    MPI_Comm_dup(job->comms[oldcomm]->mpi_comm, &(job->comms[newcomm]->mpi_comm));
    MPI_Comm_size(job->comms[newcomm]->mpi_comm, &(job->comms[newcomm]->nnodes));
    MPI_Comm_rank(job->comms[newcomm]->mpi_comm, &(job->comms[newcomm]->mynode));

    job->comms[newcomm]->wnranks = gk_imalloc(job->comms[newcomm]->nnodes, "wnranks");
    gk_icopy(job->comms[newcomm]->nnodes, job->comms[oldcomm]->wnranks, 
        job->comms[newcomm]->wnranks);

    /* inter<->intra component */
    job->comms[newcomm]->gsize   = job->comms[oldcomm]->gsize;

    job->comms[newcomm]->slvdist = gk_imalloc(job->comms[newcomm]->nnodes+1, "slvdist");
    gk_icopy(job->comms[newcomm]->nnodes+1, job->comms[oldcomm]->slvdist, 
        job->comms[newcomm]->slvdist);

    job->comms[newcomm]->g2nmap = gk_imalloc(job->comms[newcomm]->gsize, "g2nmap");
    gk_icopy(job->comms[newcomm]->gsize, job->comms[oldcomm]->g2nmap, 
        job->comms[newcomm]->g2nmap);

    job->comms[newcomm]->g2lmap = gk_imalloc(job->comms[newcomm]->gsize, "g2lmap");
    gk_icopy(job->comms[newcomm]->gsize, job->comms[oldcomm]->g2lmap, 
        job->comms[newcomm]->g2lmap);

    job->comms[newcomm]->l2gmap = gk_imalloc(job->comms[newcomm]->lsize, "l2gmap");
    gk_icopy(job->comms[newcomm]->lsize, job->comms[oldcomm]->l2gmap, 
        job->comms[newcomm]->l2gmap);


    /* get the commid from the master node */
    if (job->comms[newcomm]->mynode == 0) 
      job->comms[newcomm]->mpi_commid = mnode_get_next_commid(job);
    MPI_Bcast(&(job->comms[newcomm]->mpi_commid), 1, MPI_INT, 0, job->comms[newcomm]->mpi_comm);

    /* lock */
    job->comms[newcomm]->rwlock = (pthread_rwlock_t *)gk_malloc(sizeof(pthread_rwlock_t), "comm->lock");
    BDASSERT(pthread_rwlock_init(job->comms[newcomm]->rwlock, NULL) == 0);


    /* slave component */
    scomm.mcomm = job->comms[newcomm]->mcomm;
    scomm.size  = job->comms[newcomm]->gsize;
    scomm.lsize = job->comms[newcomm]->lsize;
    scomm.nsize = job->comms[newcomm]->nnodes;
    scomm.nrank = job->comms[newcomm]->mynode;
    scomm.rrank = job->comms[newcomm]->l2gmap[0];

    BD_GET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

    /* reset the counter */
    job->comms[oldcomm]->counter = job->comms[oldcomm]->lsize;

    /* all processes have called, unblock them */
    for (i=0; i<job->comms[newcomm]->lsize; i++) {
      scomm.rank  = job->comms[newcomm]->l2gmap[i];
      scomm.lrank = i;

      srank = job->comms[newcomm]->sranks[i];
      xfer_out_scb(job->scbs[srank], &scomm, sizeof(bdscomm_t), BDMPI_BYTE);
      slvpool_cunblock(job, srank);
    }
  }

  M_IFSET(BDMPI_DBG_IPCM, 
      bdprintf("[MSTR%04d] comm_dup: mcomm: %d, counter: %d [exiting]\n", 
        job->mynode, msg->mcomm, job->comms[msg->mcomm]->counter));

  BD_LET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

  gk_free((void **)&arg, LTERM);

  return NULL;
}


/*************************************************************************/
/*! Response to a BDMPI_Comm_free.
    Protocol:
       Blocks the process.
       Decreases the counter associated with the counter.
       If counter becomes 0, then the master communicator is deleted.
       Sets the counter back to the size of that communicator.
*/
/*************************************************************************/
void *mstr_comm_free(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  size_t i;
  int comm;

  /* block the slave */
  slvpool_cblock(job, babel_get_srank(job->comms[msg->mcomm], msg->myrank));

  BD_GET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  M_IFSET(BDMPI_DBG_IPCM, 
      bdprintf("[MSTR%04d.%04d] comm_free: mcomm: %d, counter: %d [entering]\n", 
        job->mynode, msg->myrank, msg->mcomm, job->comms[msg->mcomm]->counter));

  if (--job->comms[msg->mcomm]->counter == 0) {
    comm = msg->mcomm;

    /* all processes have called, unblock them */
    for (i=0; i<job->comms[comm]->lsize; i++) 
      slvpool_cunblock(job, job->comms[comm]->sranks[i]);

    M_IFSET(BDMPI_DBG_IPCM, 
        bdprintf("[MSTR%04d.%04d] comm_free: mcomm: %d [exiting]\n", 
            job->mynode, msg->myrank, msg->mcomm));

    BD_LET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

    comm_free(job, comm);
  }
  else {
    M_IFSET(BDMPI_DBG_IPCM, 
        bdprintf("[MSTR%04d.%04d] comm_free: mcomm: %d [exiting]\n", 
            job->mynode, msg->myrank, msg->mcomm));

    BD_LET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */
  }

  gk_free((void **)&arg, LTERM);

  return NULL;
}


/*************************************************************************/
/*! Response to a BDMPI_Comm_split.
    Protocol:
       Blocks the process.
       Decreases the counter associated with the counter.
       If counter becomes 0, then a new master communicator is created 
       Sets the counter back to the size of that communicator.
       Copies the per-slave communicator info in the per-slave scbs

    Note:
       msg->source = color
       msg->dest   = key
*/
/*************************************************************************/
void *mstr_comm_split(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  size_t i, j, k;
  int oldnnodes, oldsize, oldlsize, newsize, newlsize, srank, lrank, newcnum, itmp;
  bdmcomm_t *oldcomm, *newcomm;
  bdscomm_t scomm;
  splitdata_t *skeys, *allskeys;
  int *counts, *displs, *oldslvdist;
  MPI_Comm null_mpi_comm;

  /*TODO: There may be a deadlock in this call due to the persistent WRLOCK */

  /* block it */
  slvpool_cblock(job, babel_get_srank(job->comms[msg->mcomm], msg->myrank));

  BD_GET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  M_IFSET(BDMPI_DBG_IPCM, 
      bdprintf("[MSTR%04d.%04d] comm_split: mcomm: %d, counter: %d [entering]\n", 
        job->mynode, msg->myrank, msg->mcomm, job->comms[msg->mcomm]->counter));

  oldcomm = job->comms[msg->mcomm];

  /* allocate memory for colors/keys */
  if (oldcomm->counter == oldcomm->lsize) 
    oldcomm->skeys = (splitdata_t *)gk_malloc(oldcomm->lsize*sizeof(splitdata_t), "skeys");

  /* populate the colors & keys */
  lrank = babel_get_lrank(oldcomm, msg->myrank);
  oldcomm->skeys[lrank].color = msg->source;
  oldcomm->skeys[lrank].key   = msg->dest;
  oldcomm->skeys[lrank].rank  = msg->myrank;


  if (--oldcomm->counter == 0) { /* time to do some work! */
    oldcomm->counter = oldcomm->lsize;

    oldnnodes  = oldcomm->nnodes;
    oldslvdist = oldcomm->slvdist;
    skeys      = oldcomm->skeys;
    oldsize    = oldcomm->gsize;
    oldlsize   = oldcomm->lsize;

    /* get all information to all processors */
    allskeys = (splitdata_t *)gk_malloc(oldsize*sizeof(splitdata_t), "allskeys");
    counts   = gk_imalloc(oldnnodes+1, "counts");
    displs   = gk_imalloc(oldnnodes, "displs");
    for (i=0; i<oldnnodes; i++) {
      counts[i] = (oldslvdist[i+1]-oldslvdist[i])*sizeof(splitdata_t);
      displs[i] = oldslvdist[i]*sizeof(splitdata_t);
    }

    BDASSERT(MPI_Allgatherv(skeys, sizeof(splitdata_t)*oldlsize, MPI_BYTE,
                   allskeys, counts, displs, MPI_BYTE, oldcomm->mpi_comm)
             == MPI_SUCCESS);

    /*
    bdprintf("[%04d] oldlsize: %d, oldsize: %d\n", job->mynode, oldlsize, oldsize);
    for (i=0; i<oldlsize; i++)
      bdprintf("[%04d] skeys[%zu]: %d %d %d\n", job->mynode, i, skeys[i].color, skeys[i].key, skeys[i].rank);
    for (i=0; i<oldsize; i++)
      bdprintf("[%04d] allskeys[%zu]: %d %d %d\n", job->mynode, i, allskeys[i].color, allskeys[i].key, allskeys[i].rank);
    */

    comm_sort_splitdata(oldsize, allskeys);

    /*
    for (i=0; i<oldsize; i++)
      bdprintf("[%04d] allskeys[%zu]: %d %d %d\n", job->mynode, i, allskeys[i].color, allskeys[i].key, allskeys[i].rank);
    */

    k = 0;
    while (k < oldsize) {
      newlsize = 0;
      if (babel_is_local(oldcomm, allskeys[k].rank)) 
        newlsize++;
        
      for (j=k+1; j<oldsize; j++) {
        if (allskeys[j].color != allskeys[k].color)
          break;

        if (babel_is_local(oldcomm, allskeys[j].rank)) 
          newlsize++;
      }
      newsize = j-k;

      /*
      bdprintf("[%04d] color: %d, newlsize: %d, newsize: %d\n", 
          job->mynode, allskeys[k].color, newlsize, newsize);
      */

      if (newlsize > 0) {
        BD_GET_LOCK(job->comm_lock);
        if (job->nfreecomm == 0)
          comm_realloc_master(job);

        newcnum = job->commfreelist[--job->nfreecomm];
        BD_LET_LOCK(job->comm_lock);

        job->comms[newcnum] = newcomm = (bdmcomm_t *)gk_malloc(sizeof(bdmcomm_t), "comms[newcnum]");
        memset(newcomm, 0, sizeof(bdmcomm_t));

        /* intra-node component */
        newcomm->mcomm   = newcnum;
        newcomm->lsize   = newlsize;
        newcomm->counter = newlsize;
        newcomm->copid   = 0;

        /* intra<->inter component */
        newcomm->gsize  = newsize;
        newcomm->g2lmap = gk_imalloc(newcomm->gsize, "newcomm->g2lmap");
        newcomm->g2nmap = gk_imalloc(newcomm->gsize, "newcomm->g2nmap");
        newcomm->l2gmap = gk_imalloc(newcomm->lsize, "newcomm->l2gmap");
        newcomm->sranks = gk_imalloc(newcomm->lsize, "newcomm->sranks");

        /* create g2lmap and l2gmap */
        gk_iset(oldnnodes, 0, counts);
        for (j=0, i=0; i<newsize; i++) {
          newcomm->g2lmap[i] = counts[oldcomm->g2nmap[allskeys[k+i].rank]]++;
          if (babel_is_local(oldcomm, allskeys[k+i].rank)) {
            newcomm->l2gmap[j] = i;
            newcomm->sranks[j] = babel_get_srank(oldcomm, allskeys[k+i].rank);
            j++;
          }
        }

        /* determine the number of nodes in the new communicator */
        for (j=0, i=0; i<oldnnodes; i++) {
          if (counts[i] > 0)
            j++;
        }
        newcomm->nnodes = j;

        /* create slvdist/wnranks for the new communicator */
        newcomm->slvdist = gk_imalloc(newcomm->nnodes+1, "slvdist");
        newcomm->wnranks = gk_imalloc(newcomm->nnodes, "wnranks");
        for (j=0, i=0; i<oldnnodes; i++) {
          if (counts[i] > 0) {
            newcomm->slvdist[j] = counts[i];
            newcomm->wnranks[j] = oldcomm->wnranks[i];
            j++;
          }
        }
        MAKECSR(i, newcomm->nnodes, newcomm->slvdist);

        /*
        for (i=0; i<newcomm->nnodes+1; i++) {
          bdprintf("[%04d] color: %d, newcomm->slvdist[%zu] = %d\n",
              job->mynode, allskeys[k].color, i, newcomm->slvdist[i]);
        }
        for (i=0; i<newcomm->nnodes; i++) {
          bdprintf("[%04d] color: %d, newcomm->wnranks[%zu] = %d\n",
              job->mynode, allskeys[k].color, i, newcomm->wnranks[i]);
        }
        */

        /* create the g2nmap for the new communicator */
        MAKECSR(i, oldnnodes, counts);
        for (j=0, i=0; i<oldnnodes; i++) {
          if (counts[i+1]-counts[i] > 0)
            counts[i] = j++;
        }
        for (i=0; i<newsize; i++) 
          newcomm->g2nmap[i] = counts[oldcomm->g2nmap[allskeys[k+i].rank]];

        newcomm->mynode = counts[oldcomm->mynode];

        /*
        bdprintf("[%04d] color: %d, newcomm->nnodes: %d, newcomm->mynode: %d\n",
            job->mynode, allskeys[k].color, newcomm->nnodes, newcomm->mynode);
        */
      }

      /* inter-node component */
      if (newlsize > 0) 
        BDASSERT(MPI_Comm_split(oldcomm->mpi_comm, 1, 0, &(newcomm->mpi_comm)) 
            == MPI_SUCCESS);
      else 
        BDASSERT(MPI_Comm_split(oldcomm->mpi_comm, MPI_UNDEFINED, 0, &null_mpi_comm) 
            == MPI_SUCCESS);

      if (newlsize > 0) {
        MPI_Comm_size(newcomm->mpi_comm, &itmp);
        BDASSERT(itmp == newcomm->nnodes);
        MPI_Comm_rank(newcomm->mpi_comm, &itmp);
        BDASSERT(itmp == newcomm->mynode);

        /* get the commid from the master node */
        if (newcomm->mynode == 0) 
          newcomm->mpi_commid = mnode_get_next_commid(job);
        MPI_Bcast(&(newcomm->mpi_commid), 1, MPI_INT, 0, newcomm->mpi_comm);

        /*
        bdprintf("[%04d] color: %d, nnodes: %d, mynode: %d, mcomm: %d, size: %d, lsize: %d, newcomm->mpi_commid: %d\n",
            job->mynode, allskeys[k].color, 
            newcomm->nnodes, newcomm->mynode, newcomm->mcomm, 
            newcomm->slvdist[newcomm->nnodes], newcomm->lsize, 
            newcomm->mpi_commid);
        */

        /* lock */
        newcomm->rwlock = (pthread_rwlock_t *)gk_malloc(sizeof(pthread_rwlock_t), "comm->lock");
        BDASSERT(pthread_rwlock_init(newcomm->rwlock, NULL) == 0);

        /* slave component */
        scomm.mcomm = newcomm->mcomm;
        scomm.size  = newcomm->gsize;
        scomm.lsize = newcomm->lsize;
        scomm.nsize = newcomm->nnodes;
        scomm.nrank = newcomm->mynode;
        scomm.rrank = newcomm->l2gmap[0];

        /* unlock the processes involed in that sub-communicator */
        for (i=0; i<newlsize; i++) {
          scomm.rank  = newcomm->l2gmap[i];
          scomm.lrank = i;
          srank = newcomm->sranks[i];

          xfer_out_scb(job->scbs[srank], &scomm, sizeof(bdscomm_t), BDMPI_BYTE);
          slvpool_cunblock(job, srank);  
        }
      }

      k += newsize;
    }

    /* free the data that you allocated */
    gk_free((void **)&oldcomm->skeys, &allskeys, &counts, &displs, LTERM);
  }

  M_IFSET(BDMPI_DBG_IPCM, 
      bdprintf("[MSTR%04d.%04d] comm_split: mcomm: %d, counter: %d [exiting]\n", 
        job->mynode, msg->myrank, msg->mcomm, job->comms[msg->mcomm]->counter));

  BD_LET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

  gk_free((void **)&arg, LTERM);

  return NULL;
}


/*************************************************************************/
/*! Sorts the elements of splitdata_t in increasing color/key/rank order */
/*************************************************************************/
void comm_sort_splitdata(size_t n, splitdata_t *skeys)
{
#define splitdata_lt(a, b) \
  ((a)->color < (b)->color \
   || ((a)->color == (b)->color && (a)->key < (b)->key) \
   || ((a)->color == (b)->color && (a)->key == (b)->key && (a)->rank < (b)->rank))
GK_MKQSORT(splitdata_t, skeys, n, splitdata_lt);
#undef splitdata_lt
}
