/*!
\file
\brief Various functions for dealing with message headers
\date Started 4/4/2013
\author George

TODO: You need to introduce a collective-operation-ID in order to eliminate
      the need for post and pre barriers.
*/

#include "bdmprun.h"



/*************************************************************************/
/*! Setup the per-slave FIFO message header queues */
/*************************************************************************/
void pending_setup(mjob_t *job)
{
  int i;
  pthread_mutexattr_t mtx_attr;

  job->psends = (header_t **)gk_malloc(sizeof(header_t *)*job->ns, "job->psends");
  for (i=0; i<job->ns; i++) 
    job->psends[i] = NULL;

  BDASSERT(pthread_mutexattr_init(&mtx_attr) == 0);
  BDASSERT(pthread_mutexattr_settype(&mtx_attr, PTHREAD_MUTEX_RECURSIVE) == 0);

  job->plocks = (pthread_mutex_t **)gk_malloc(sizeof(pthread_mutex_t *)*job->ns, "job->plocks");
  for (i=0; i<job->ns; i++) {
    job->plocks[i] = (pthread_mutex_t *)gk_malloc(sizeof(pthread_mutex_t), "plocks[i]");
    BDASSERT(pthread_mutex_init(job->plocks[i], &mtx_attr) == 0);
  }

  return;
}


/*************************************************************************/
/*! Cleanup the per-slave message headers */
/*************************************************************************/
void pending_cleanup(mjob_t *job)
{
  int i;
  header_t *curr, *next;

  for (i=0; i<job->ns; i++) {
    curr = job->psends[i];
    while (curr != NULL) {
      next = curr->next;
      pending_freeheader(job, &curr);
    }

    pthread_mutex_destroy(job->plocks[i]);
    gk_free((void **)&job->plocks[i], LTERM);
  }
  gk_free((void **)&job->psends, &job->plocks, LTERM);

  return;
}


/*************************************************************************/
/*! Free a header and associated info */
/*************************************************************************/
void pending_freeheader(mjob_t *job, header_t **r_hdr)
{
  header_t *hdr; 

  if (*r_hdr == NULL)
    return;

  hdr = *r_hdr;
  if (hdr->buf) {
    if (hdr->mlck)
      munlock(hdr->buf, hdr->len);

    gk_free((void **)&hdr->buf, LTERM);
  }
  gk_free((void **)r_hdr, LTERM);

  return;
}


/*************************************************************************/
/*! Adds a header corresponding to a send request. 
    \param msg is a send request. 
*/
/*************************************************************************/
void pending_addsend(mjob_t *job, bdmsg_t *msg, void *buf, size_t len)
{
  header_t *curr, *prev=NULL;

  /* get the global rank of the queue that will store the header */
  int qnum = babel_get_srank(job->comms[msg->mcomm], msg->dest); 
  BD_GET_LOCK(job->plocks[qnum]);

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] pending_addsend: from: %3d, to: %3d, tag: %3d, mcomm: %3d\n",
        job->mynode, msg->source, msg->dest, msg->tag, msg->mcomm));

  /* find the end of the queue */
  curr = job->psends[qnum];
  while (curr != NULL) {
    prev = curr;
    curr = curr->next;
  }

  curr = (header_t *)gk_malloc(sizeof(header_t), "curr");
  curr->msg  = *msg;
  curr->buf  = buf;
  curr->len  = len;
  curr->mlck = 0;
  curr->next = NULL;

  if (prev == NULL)
    job->psends[qnum] = curr;
  else
    prev->next = curr;

  BD_LET_LOCK(job->plocks[qnum]);

  return;
}

/*************************************************************************/
/*! Finds and extracts a pending send header. 
    \param msg is a recv message.
*/
/*************************************************************************/
header_t *pending_getsend(mjob_t *job, bdmsg_t *msg, int rmheader)
{
  header_t *curr, *prev=NULL;

  /* get the global rank of the queue that will store the header */
  int qnum = babel_get_srank(job->comms[msg->mcomm], msg->dest); 
  BD_GET_LOCK(job->plocks[qnum]);

  curr = job->psends[qnum];
  while (curr != NULL) {
    if (curr->msg.msgtype == BDMPI_MSGTYPE_SEND &&
        curr->msg.mcomm == msg->mcomm &&
        (msg->source == BDMPI_ANY_SOURCE || curr->msg.source == msg->source) &&  
        (msg->tag == BDMPI_ANY_TAG || curr->msg.tag == msg->tag)) { /* match */

      if (rmheader) { /* remove header from the list, if requested */
        if (prev == NULL)
          job->psends[qnum] = curr->next;
        else
          prev->next = curr->next;
      }

      break;
    }
    prev = curr;
    curr = curr->next;
  }

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] pending_getsend: from: %3d, to: %3d, tag: %3d, comm: %3d, status: %d\n",
        job->mynode, msg->source, msg->dest, msg->tag, msg->mcomm, (curr==NULL?0:1)));

  BD_LET_LOCK(job->plocks[qnum]);

  return curr;
}

/*************************************************************************/
/*! Locks/unlocks the appropriate psend link list. 
    \param msg is a send request. 
*/
/*************************************************************************/
void pending_locksend(mjob_t *job, bdmsg_t *msg)
{
  /* get the global rank of the queue that will store the header */
  int qnum = babel_get_srank(job->comms[msg->mcomm], msg->dest); 
  BD_GET_LOCK(job->plocks[qnum]);

  return;
}

void pending_unlocksend(mjob_t *job, bdmsg_t *msg)
{
  /* get the global rank of the queue that will store the header */
  int qnum = babel_get_srank(job->comms[msg->mcomm], msg->dest); 
  BD_LET_LOCK(job->plocks[qnum]);

  return;
}



/*************************************************************************/
/*! Adds a header corresponding to a bcast operation. 
    \param msg is a bcast request. 

    \note The headers are put in the psends[root]
*/
/*************************************************************************/
void pending_addbcast(mjob_t *job, bdmsg_t *msg, void *buf, size_t len, int icnt)
{
  header_t *curr, *prev=NULL;

  /* the header is always stored in the queue of lrank==0. */
  int qnum = job->comms[msg->mcomm]->sranks[0]; 
  BD_GET_LOCK(job->plocks[qnum]);

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] AddBcast: from: %3d, to: %3d, mcomm: %3d\n",
        job->mynode, msg->source, msg->dest, msg->mcomm));

  /* find the end of the queue for the root */
  curr = job->psends[qnum];
  while (curr != NULL) {
    prev = curr;
    curr = curr->next;
  }

  curr = (header_t *)gk_malloc(sizeof(header_t), "curr");
  curr->msg     = *msg;
  curr->buf     = buf;
  curr->len     = len;
  curr->mlck    = BDMPI_MLOCK_BCAST;
  curr->counter = icnt;
  curr->next    = NULL;

  if (curr->mlck)
    mlock(curr->buf, curr->len);

  if (prev == NULL)
    job->psends[qnum] = curr;
  else
    prev->next = curr;

  BD_LET_LOCK(job->plocks[qnum]);

  return;
}

/*************************************************************************/
/*! Finds and extracts a pending bcast header. 
    \param msg is a bcast message.

    \note The headers are located in the psends[root]
*/
/*************************************************************************/
header_t *pending_getbcast(mjob_t *job, bdmsg_t *msg, int countdown, 
              int *r_counter)
{
  header_t *curr, *prev=NULL;

  /* the header is always stored in the queue of lrank==0. */
  int qnum = job->comms[msg->mcomm]->sranks[0]; 
  BD_GET_LOCK(job->plocks[qnum]);

  curr = job->psends[qnum];
  while (curr != NULL) {
    if (curr->msg.msgtype == BDMPI_MSGTYPE_BCASTI &&
        curr->msg.copid == msg->copid &&
        curr->msg.mcomm == msg->mcomm)  { /* match */

      if (countdown) {
        if (--curr->counter == 0) { /* delete if counter became 0 */
          if (prev == NULL)
            job->psends[qnum] = curr->next;
          else
            prev->next = curr->next;
        }
      }

      break;
    }
    prev = curr;
    curr = curr->next;
  }

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] GetBcast: from: %3d, to: %3d, comm: %3d, status: %d\n",
        job->mynode, msg->source, msg->dest, msg->mcomm, (curr==NULL?0:1)));

  if (r_counter != NULL)
    *r_counter = curr->counter;

  BD_LET_LOCK(job->plocks[qnum]);

  return curr;
}



/*************************************************************************/
/*! Adds a header corresponding to a reduce operation. 
    \param msg is a reduce request. 

    \note The headers are put in the psends[root]
*/
/*************************************************************************/
void pending_addreduce(mjob_t *job, bdmsg_t *msg, void *buf, size_t len, int icnt)
{
  header_t *curr, *prev=NULL;

  /* the header is always stored in the queue of lrank==0. */
  int qnum = job->comms[msg->mcomm]->sranks[0]; 
  BD_GET_LOCK(job->plocks[qnum]);

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] AddReduce: from: %3d, to: %3d, mcomm: %3d\n",
        job->mynode, msg->source, msg->dest, msg->mcomm));

  /* find the end of the queue for the root */
  curr = job->psends[qnum];
  while (curr != NULL) {
    prev = curr;
    curr = curr->next;
  }

  curr = (header_t *)gk_malloc(sizeof(header_t), "curr");
  curr->msg     = *msg;
  curr->buf     = buf;
  curr->len     = len;
  curr->mlck    = BDMPI_MLOCK_REDUCE;
  curr->counter = icnt;
  curr->next    = NULL;

  if (curr->mlck)
    mlock(curr->buf, curr->len);

  if (prev == NULL)
    job->psends[qnum] = curr;
  else
    prev->next = curr;

  BD_LET_LOCK(job->plocks[qnum]);

  return;
}

/*************************************************************************/
/*! Finds and extracts a pending reduce header. 
    \param msg is a reduce message.

    \note The headers are located in the psends[root]
*/
/*************************************************************************/
header_t *pending_getreduce(mjob_t *job, bdmsg_t *msg)
{
  header_t *curr, *prev=NULL;

  /* the header is always stored in the queue of lrank==0. */
  int qnum = job->comms[msg->mcomm]->sranks[0]; 
  BD_GET_LOCK(job->plocks[qnum]);

  curr = job->psends[qnum];
  while (curr != NULL) {
    if (curr->msg.msgtype == BDMPI_MSGTYPE_REDUCEI &&
        curr->msg.copid == msg->copid &&
        curr->msg.mcomm == msg->mcomm) { /* match */
      if (--curr->counter == 0) { /* delete if counter became 0 */
        if (prev == NULL)
          job->psends[qnum] = curr->next;
        else
          prev->next = curr->next;
      }
      break;
    }
    prev = curr;
    curr = curr->next;
  }

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] GetReduce: from: %3d, to: %3d, mcomm: %3d, status: %d\n",
        job->mynode, msg->source, msg->dest, msg->mcomm, (curr==NULL?0:1)));

  BD_LET_LOCK(job->plocks[qnum]);

  return curr;
}



/*************************************************************************/
/*! Adds a header corresponding to a merge operation. 
    \param msg is a merge request. 

    \note The headers are put in the psends[root]
*/
/*************************************************************************/
void pending_addmerge(mjob_t *job, bdmsg_t *msg, void *buf, int icnt)
{
  header_t *curr, *prev=NULL;

  /* the header is always stored in the queue of lrank==0. */
  int qnum = job->comms[msg->mcomm]->sranks[0]; 
  BD_GET_LOCK(job->plocks[qnum]);

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] AddMerge: from: %3d, to: %3d, mcomm: %3d\n",
        job->mynode, msg->source, msg->dest, msg->mcomm));

  /* find the end of the queue for the root */
  curr = job->psends[qnum];
  while (curr != NULL) {
    prev = curr;
    curr = curr->next;
  }

  curr = (header_t *)gk_malloc(sizeof(header_t), "curr");
  curr->msg     = *msg;
  curr->buf     = buf;
  curr->len     = 0;
  curr->mlck    = 0;
  curr->counter = icnt;
  curr->next    = NULL;

  if (prev == NULL)
    job->psends[qnum] = curr;
  else
    prev->next = curr;

  BD_LET_LOCK(job->plocks[qnum]);

  return;
}

/*************************************************************************/
/*! Finds and extracts a pending merge header. 
    \param msg is a merge message.

    \note The headers are located in the psends[root]
*/
/*************************************************************************/
header_t *pending_getmerge(mjob_t *job, bdmsg_t *msg)
{
  header_t *curr, *prev=NULL;

  /* the header is always stored in the queue of lrank==0. */
  int qnum = job->comms[msg->mcomm]->sranks[0]; 
  BD_GET_LOCK(job->plocks[qnum]);

  curr = job->psends[qnum];
  while (curr != NULL) {
    if (curr->msg.msgtype == BDMPI_MSGTYPE_MERGEI &&
        curr->msg.copid == msg->copid &&
        curr->msg.mcomm == msg->mcomm) { /* match */
      if (--curr->counter == 0) { /* delete if counter became 0 */
        if (prev == NULL)
          job->psends[qnum] = curr->next;
        else
          prev->next = curr->next;
      }
      break;
    }
    prev = curr;
    curr = curr->next;
  }

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] GetMerge: from: %3d, to: %3d, mcomm: %3d, status: %d\n",
        job->mynode, msg->source, msg->dest, msg->mcomm, (curr==NULL?0:1)));

  BD_LET_LOCK(job->plocks[qnum]);

  return curr;
}



/*************************************************************************/
/*! Adds a header corresponding to an allreduce operation. 
    \param msg is an allreduce request. 

    \note The headers are put in the psends[0]
*/
/*************************************************************************/
void pending_addallreduce(mjob_t *job, bdmsg_t *msg, void *buf, size_t len)
{
  header_t *curr, *prev=NULL;

  /* the header is always stored in the queue of lrank==0. */
  int qnum = job->comms[msg->mcomm]->sranks[0]; 
  BD_GET_LOCK(job->plocks[qnum]);

  /* find the end of the queue for the root */
  curr = job->psends[qnum];
  while (curr != NULL) {
    prev = curr;
    curr = curr->next;
  }

  curr = (header_t *)gk_malloc(sizeof(header_t), "curr");
  curr->msg     = *msg;
  curr->buf     = buf;
  curr->len     = len;
  curr->mlck    = BDMPI_MLOCK_REDUCE;
  curr->counter = 2*job->comms[msg->mcomm]->lsize; 
  curr->next    = NULL;

  if (curr->mlck)
    mlock(curr->buf, curr->len);

  if (prev == NULL)
    job->psends[qnum] = curr;
  else
    prev->next = curr;

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] pending_addallreduce: myrank: %3d, mcomm: %3d, curr: %p\n",
        job->mynode, msg->myrank, msg->mcomm, (void *)curr));

  BD_LET_LOCK(job->plocks[qnum]);

  return;
}

/*************************************************************************/
/*! Finds and extracts a pending allreduce header. 
    \param msg is an allreduce message.

    \note The headers are located in the psends[0]
*/
/*************************************************************************/
header_t *pending_getallreduce(mjob_t *job, bdmsg_t *msg, int *r_counter)
{
  header_t *curr, *prev=NULL;

  /* the header is always stored in the queue of lrank==0. */
  int qnum = job->comms[msg->mcomm]->sranks[0]; 
  BD_GET_LOCK(job->plocks[qnum]);

  curr = job->psends[qnum];
  while (curr != NULL) {
    if (curr->msg.msgtype == BDMPI_MSGTYPE_ALLREDUCEI &&
        curr->msg.copid == msg->copid &&
        curr->msg.mcomm == msg->mcomm) { /* match */
      if (--curr->counter == 0) { /* delete if counter became 0 */
        if (prev == NULL)
          job->psends[qnum] = curr->next;
        else
          prev->next = curr->next;
      }
      break;
    }
    prev = curr;
    curr = curr->next;
  }

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] pending_getallreduce: myrank: %3d, mcomm: %3d, curr: %p, status: %d\n",
        job->mynode, msg->myrank, msg->mcomm, (void *)curr, (curr==NULL?0:1)));

  if (r_counter != NULL)
    *r_counter = curr->counter;

  BD_LET_LOCK(job->plocks[qnum]);

  return curr;
}



/*************************************************************************/
/*! Adds a header corresponding to an alltoall operation. 
    \param msg is a alltoall request. 

    \note The headers are put in the psends[msg->dest]
*/
/*************************************************************************/
void pending_addalltoall(mjob_t *job, bdmsg_t *msg, void *buf, size_t len)
{
  header_t *curr, *prev=NULL;

  /* get the global rank of the queue that will store the header */
  int qnum = babel_get_srank(job->comms[msg->mcomm], msg->dest);
  BD_GET_LOCK(job->plocks[qnum]);

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] AddAlltoall: from: %3d, to: %3d, mcomm: %3d\n",
        job->mynode, msg->myrank, msg->dest, msg->mcomm));

  /* find the end of the queue for the root */
  curr = job->psends[qnum];
  while (curr != NULL) {
    prev = curr;
    curr = curr->next;
  }

  curr = (header_t *)gk_malloc(sizeof(header_t), "curr");
  curr->msg  = *msg;
  curr->buf  = buf;
  curr->len  = len;
  curr->mlck = (len == 0 ? 0 : BDMPI_MLOCK_ALLTOALL);
  curr->next = NULL;

  if (curr->mlck)
    mlock(curr->buf, curr->len);

  if (prev == NULL)
    job->psends[qnum] = curr;
  else
    prev->next = curr;

  BD_LET_LOCK(job->plocks[qnum]);

  return;
}

/*************************************************************************/
/*! Finds and extracts a pending alltoall header 
    \param msg is a send message.
*/
/*************************************************************************/
header_t *pending_getalltoall(mjob_t *job, bdmsg_t *msg)
{
  header_t *curr, *prev=NULL;

  /* get the global rank of the queue that will store the header */
  int qnum = babel_get_srank(job->comms[msg->mcomm], msg->myrank);
  BD_GET_LOCK(job->plocks[qnum]);

  curr = job->psends[qnum];
  while (curr != NULL) {
    if (curr->msg.msgtype == BDMPI_MSGTYPE_ALLTOALLI &&
        msg->copid == curr->msg.copid &&
        msg->mcomm == curr->msg.mcomm) { /* match */
      if (prev == NULL)
        job->psends[qnum] = curr->next;
      else
        prev->next = curr->next;
      break;
    }
    prev = curr;
    curr = curr->next;
  }

  BD_LET_LOCK(job->plocks[qnum]);

  return curr;
}



/*************************************************************************/
/*! Adds a header corresponding to a scatter operation. 
    \param msg is a scatter request. 

    \note The headers are put in the psends[msg->dest]
*/
/*************************************************************************/
void pending_addscatter(mjob_t *job, bdmsg_t *msg, void *buf, size_t len)
{
  header_t *curr, *prev=NULL;

  /* get the global rank of the queue that will store the header */
  int qnum = babel_get_srank(job->comms[msg->mcomm], msg->dest);
  BD_GET_LOCK(job->plocks[qnum]);

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] AddScatter: from: %3d, to: %3d, mcomm: %3d\n",
        job->mynode, msg->myrank, msg->dest, msg->mcomm));

  /* find the end of the queue for the root */
  curr = job->psends[qnum];
  while (curr != NULL) {
    prev = curr;
    curr = curr->next;
  }

  curr = (header_t *)gk_malloc(sizeof(header_t), "curr");
  curr->msg  = *msg;
  curr->buf  = buf;
  curr->len  = len;
  curr->mlck = (len == 0 ? 0 : BDMPI_MLOCK_SCATTER);
  curr->next = NULL;

  if (curr->mlck)
    mlock(curr->buf, curr->len);

  if (prev == NULL)
    job->psends[qnum] = curr;
  else
    prev->next = curr;

  BD_LET_LOCK(job->plocks[qnum]);

  return;
}

/*************************************************************************/
/*! Finds and extracts a pending scatter header 
    \param msg is a send message.
*/
/*************************************************************************/
header_t *pending_getscatter(mjob_t *job, bdmsg_t *msg)
{
  header_t *curr, *prev=NULL;

  /* get the global rank of the queue that will store the header */
  int qnum = babel_get_srank(job->comms[msg->mcomm], msg->myrank);
  BD_GET_LOCK(job->plocks[qnum]);

  curr = job->psends[qnum];
  while (curr != NULL) {
    if (curr->msg.msgtype == BDMPI_MSGTYPE_SCATTERI &&
        msg->copid == curr->msg.copid &&
        msg->mcomm == curr->msg.mcomm &&
        msg->source == curr->msg.source) { /* match */
      if (prev == NULL)
        job->psends[qnum] = curr->next;
      else
        prev->next = curr->next;
      break;
    }
    prev = curr;
    curr = curr->next;
  }

  BD_LET_LOCK(job->plocks[qnum]);

  return curr;
}



/*************************************************************************/
/*! Adds a header corresponding to a gather operation. 
    \param msg is a gather request. 

    \note The headers are put in the psends[msg->dest]
*/
/*************************************************************************/
void pending_addgather(mjob_t *job, bdmsg_t *msg, void *buf, size_t len)
{
  header_t *curr, *prev=NULL;

  /* get the global rank of the queue that will store the header */
  int qnum = babel_get_srank(job->comms[msg->mcomm], msg->dest);
  BD_GET_LOCK(job->plocks[qnum]);

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] AddGather: from: %3d, to: %3d, mcomm: %3d\n",
        job->mynode, msg->myrank, msg->dest, msg->mcomm));

  /* find the end of the queue for the root */
  curr = job->psends[qnum];
  while (curr != NULL) {
    prev = curr;
    curr = curr->next;
  }

  curr = (header_t *)gk_malloc(sizeof(header_t), "curr");
  curr->msg  = *msg;
  curr->buf  = buf;
  curr->len  = len;
  curr->mlck = (len == 0 ? 0 : BDMPI_MLOCK_GATHER);
  curr->next = NULL;

  if (curr->mlck)
    mlock(curr->buf, curr->len);

  if (prev == NULL)
    job->psends[qnum] = curr;
  else
    prev->next = curr;

  BD_LET_LOCK(job->plocks[qnum]);

  return;
}

/*************************************************************************/
/*! Finds and extracts a pending gather header 
    \param msg is a send message.
*/
/*************************************************************************/
header_t *pending_getgather(mjob_t *job, bdmsg_t *msg)
{
  header_t *curr, *prev=NULL;

  /* get the global rank of the queue that will store the header */
  int qnum = babel_get_srank(job->comms[msg->mcomm], msg->myrank);
  BD_GET_LOCK(job->plocks[qnum]);

  curr = job->psends[qnum];
  while (curr != NULL) {
    if (curr->msg.msgtype == BDMPI_MSGTYPE_GATHERI &&
        msg->copid == curr->msg.copid &&
        msg->mcomm == curr->msg.mcomm) { /* match */
      if (prev == NULL)
        job->psends[qnum] = curr->next;
      else
        prev->next = curr->next;
      break;
    }
    prev = curr;
    curr = curr->next;
  }

  BD_LET_LOCK(job->plocks[qnum]);

  return curr;
}



/*************************************************************************/
/*! Adds a header corresponding to an allgather operation. 
    \param msg is a allgather request. 

    \note The headers are put in the psends[msg->myrank]
*/
/*************************************************************************/
void pending_addallgather(mjob_t *job, bdmsg_t *msg, void *buf, size_t len)
{
  header_t *curr, *prev=NULL;

  /* the header is always stored in the queue of lrank==0. */
  int qnum = job->comms[msg->mcomm]->sranks[0]; 
  BD_GET_LOCK(job->plocks[qnum]);

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] AddAllgather: myrank: %3d, source: %3d, mcomm: %3d\n",
        job->mynode, msg->myrank, msg->source, msg->mcomm));

  /* find the end of the queue for the root */
  curr = job->psends[qnum];
  while (curr != NULL) {
    prev = curr;
    curr = curr->next;
  }

  curr = (header_t *)gk_malloc(sizeof(header_t), "curr");
  curr->msg     = *msg;
  curr->buf     = buf;
  curr->len     = len;
  curr->mlck    = (len == 0 ? 0 : BDMPI_MLOCK_ALLGATHER);
  curr->counter = job->comms[msg->mcomm]->lsize;
  curr->next    = NULL;

  if (curr->mlck)
    mlock(curr->buf, curr->len);

  if (prev == NULL)
    job->psends[qnum] = curr;
  else
    prev->next = curr;

  BD_LET_LOCK(job->plocks[qnum]);

  return;
}


/*************************************************************************/
/*! Finds and extracts a pending allgather header. 
    \param msg is a allgather message.

    \note The headers are located in the psends[msg->source]
*/
/*************************************************************************/
header_t *pending_getallgather(mjob_t *job, bdmsg_t *msg, int *r_counter)
{
  header_t *curr, *prev=NULL;

  /* the header is always stored in the queue of lrank==0. */
  int qnum = job->comms[msg->mcomm]->sranks[0]; 
  BD_GET_LOCK(job->plocks[qnum]);

  curr = job->psends[qnum];
  while (curr != NULL) {
    if (curr->msg.msgtype == BDMPI_MSGTYPE_ALLGATHERI &&
        curr->msg.copid == msg->copid &&
        curr->msg.mcomm == msg->mcomm &&
        curr->msg.source == msg->source) { /* match */

      if (r_counter != NULL) {
        if (--curr->counter == 0) { /* delete if counter became 0 */
          if (prev == NULL)
            job->psends[qnum] = curr->next;
          else
            prev->next = curr->next;
        }
      }
      break;
    }
    prev = curr;
    curr = curr->next;
  }

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] GetAllgather: from: %3d, comm: %3d, status: %d\n",
        job->mynode, msg->source, msg->mcomm, (curr==NULL?0:1)));

  if (r_counter != NULL)
    *r_counter = curr->counter;
  BD_LET_LOCK(job->plocks[qnum]);

  return curr;
}




