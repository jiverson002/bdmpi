/*
 * comm.c
 *
 * Implements the various functions that deal with communicators.
 *
 * Started 4/4/2013
 * George
 *
 */


#include "bdmplib.h"


/*************************************************************************/
/* Returns the size of the communicator */
/*************************************************************************/
int bdmp_Comm_size(sjob_t *job, BDMPI_Comm comm, int *size)
{
  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }
  *size = comm->size;

  return BDMPI_SUCCESS;
}


/*************************************************************************/
/* Returns the lsize of the communicator */
/*************************************************************************/
int bdmp_Comm_lsize(sjob_t *job, BDMPI_Comm comm, int *lsize)
{
  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }
  *lsize = comm->lsize;

  return BDMPI_SUCCESS;
}


/*************************************************************************/
/* Returns the nsize of the communicator */
/*************************************************************************/
int bdmp_Comm_nsize(sjob_t *job, BDMPI_Comm comm, int *nsize)
{
  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }
  *nsize = comm->nsize;

  return BDMPI_SUCCESS;
}


/*************************************************************************/
/* Returns the rank of the calling process in the communicator */
/*************************************************************************/
int bdmp_Comm_rank(sjob_t *job, BDMPI_Comm comm, int *rank)
{
  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }
  *rank = comm->rank;

  return BDMPI_SUCCESS;
}


/*************************************************************************/
/* Returns the lrank of the calling process in the communicator */
/*************************************************************************/
int bdmp_Comm_lrank(sjob_t *job, BDMPI_Comm comm, int *lrank)
{
  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }
  *lrank = comm->lrank;

  return BDMPI_SUCCESS;
}


/*************************************************************************/
/* Returns the nrank of the calling process in the communicator */
/*************************************************************************/
int bdmp_Comm_nrank(sjob_t *job, BDMPI_Comm comm, int *nrank)
{
  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }
  *nrank = comm->nrank;

  return BDMPI_SUCCESS;
}


/*************************************************************************/
/* Returns the rrank of the calling process in the communicator */
/*************************************************************************/
int bdmp_Comm_rrank(sjob_t *job, BDMPI_Comm comm, int *rrank)
{
  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }
  *rrank = comm->rrank;

  return BDMPI_SUCCESS;
}


/*************************************************************************/
/* Duplicates a communicator */
/*************************************************************************/
int bdmp_Comm_dup(sjob_t *job, BDMPI_Comm comm, BDMPI_Comm *newcomm)
{
  size_t i;
  bdmsg_t msg, gomsg;

  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Comm_dup: entering: comm: %p\n", comm));

  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }

  memset(&msg, 0, sizeof(bdmsg_t));
  msg.msgtype = BDMPI_MSGTYPE_COMMDUP;
  msg.mcomm   = comm->mcomm;
  msg.myrank  = comm->rank;

  /* prepare to go to sleep */
  S_SB_IFSET(BDMPI_SB_SAVEALL) {
    if (job->jdesc->nr < job->jdesc->ns)
      sbma_mevictall();
  }

  /* notify the master that you entering a barrier */
  bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

  /* go to sleep... */
  BDMPL_SLEEP(job, gomsg, 1);

  /* allocate the new communicator */
  *newcomm = (bdscomm_t *)bd_malloc(sizeof(bdscomm_t), "newcomm");

  /* copy the new communicator info from the master */
  xfer_in_scb(job->scb, *newcomm, sizeof(bdscomm_t), BDMPI_BYTE);

  (*newcomm)->copid = 1;

  return BDMPI_SUCCESS;
}


/*************************************************************************/
/* Frees a communicator */
/*************************************************************************/
int bdmp_Comm_free(sjob_t *job, BDMPI_Comm *comm)
{
  bdmsg_t msg, gomsg;

  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Comm_dup: entering: comm: %p\n", *comm));

  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }

  /* prepare to go to sleep */
  S_SB_IFSET(BDMPI_SB_SAVEALL) {
    if (job->jdesc->nr < job->jdesc->ns)
      sbma_mevictall();
  }

  memset(&msg, 0, sizeof(bdmsg_t));
  msg.msgtype = BDMPI_MSGTYPE_COMMFREE;
  msg.mcomm   = (*comm)->mcomm;
  msg.myrank  = (*comm)->rank;

  /* notify the master that you entering a barrier */
  bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t));

  /* go to sleep... */
  BDMPL_SLEEP(job, gomsg, 1);

  /* remove the info associated with the old communicator */
  bd_free((void **)comm, LTERM);
  *comm = BDMPI_COMM_NULL;

  return BDMPI_SUCCESS;
}


/*************************************************************************/
/* Splits a communicator */
/*************************************************************************/
int bdmp_Comm_split(sjob_t *job, BDMPI_Comm comm, int color, int key,
          BDMPI_Comm *newcomm)
{
  bdmsg_t msg, gomsg;

  S_IFSET(BDMPI_DBG_IPCS, bdprintf("BDMPI_Comm_split: entering: comm: %p\n", comm));

  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }

  /* prepare to go to sleep */
  S_SB_IFSET(BDMPI_SB_SAVEALL) {
    if (job->jdesc->nr < job->jdesc->ns)
      sbma_mevictall();
  }

  *newcomm = BDMPI_COMM_NULL;

  memset(&msg, 0, sizeof(bdmsg_t));
  msg.msgtype = BDMPI_MSGTYPE_COMMSPLIT;
  msg.mcomm   = comm->mcomm;
  msg.myrank  = comm->rank;
  msg.source  = color;
  msg.dest    = key;

  /* notify the master that you entering a barrier */
  if (-1 == bdmq_send(job->reqMQ, &msg, sizeof(bdmsg_t))) {
    bdprintf("Failed on trying to send a req message: %s.\n",
      strerror(errno));
  }

  /* go to sleep... */
  BDMPL_SLEEP(job, gomsg, 1);

  /* create the new communicator */
  *newcomm = (bdscomm_t *)bd_malloc(sizeof(bdscomm_t), "newcomm");

  /* copy the new communicator info from the master */
  xfer_in_scb(job->scb, *newcomm, sizeof(bdscomm_t), BDMPI_BYTE);

  (*newcomm)->copid = 1;

  return BDMPI_SUCCESS;
}
