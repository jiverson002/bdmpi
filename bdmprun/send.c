/*!
\file
\brief Various functions from performing p2p send operations.
\date Started 4/6/2013
\author George
*/


#include "bdmprun.h"


/*************************************************************************/
/*! Response to a BDMPI_Send/BDMPI_Isend.
    Copies the message from the slave, saves it on the disk, and creates
    a header entry for it.
*/
/*************************************************************************/
void *mstr_send(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  bdmcomm_t *comm;
  int fd, srank, dest_node, orig_mcomm, done=1;
  char *buf=NULL;
  size_t rsize, msize;
  char *fname;

  orig_mcomm = msg->mcomm;  /* save it for later */

  BD_GET_RDLOCK(job->comms[orig_mcomm]->rwlock); /* lock communicator */

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_send: source: %d, dest: %d [entering]\n",
        job->mynode, msg->myrank, msg->source, msg->dest));

  /* hook to the key info */
  comm  = job->comms[msg->mcomm];
  srank = babel_get_srank(comm, msg->myrank);

  dest_node = babel_get_node(comm, msg->dest);

  msize = bdmp_msize(msg->count, msg->datatype);

  if (msg->fnum == -1) { /* in memory message */
    buf = gk_cmalloc(msize, "mstr_send: buf");
    xfer_in_scb(job->scbs[srank], buf, msg->count, msg->datatype);

    M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d] mstr_send: source: %d, dest: %d [done xfer_in_scb]\n",
          job->mynode, msg->source, msg->dest));
  }

  if (babel_is_local(comm, msg->dest)) { /* dest is local */
    M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d:%04d] [local] mstr_send: source: %d, dest: %d [dest_node: %d]\n",
        job->mynode, comm->mynode, msg->source, msg->dest, dest_node));

    pending_locksend(job, msg);
    pending_addsend(job, msg, buf, (buf == NULL ? 0 : msize));
    slvpool_munblock(job, babel_get_srank(comm, msg->dest));
    pending_unlocksend(job, msg);
  }
  else { /* destination is remote */
    M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d:%04d] [remote] mstr_send: source: %d, dest: %d [dest_node: l:%d, w:%d]\n",
        job->mynode, comm->mynode, msg->source, msg->dest, dest_node,
        comm->wnranks[dest_node]));

    BD_GET_LOCK(job->comm_lock);
    msg->mpi_tag = job->next_mpi_tag;
    job->next_mpi_tag += 2;
    BD_LET_LOCK(job->comm_lock);

    msg->mcomm = comm->mpi_commid; /* replace local ID with a comm-consistent ID */

    /* send the message header using the global node number of wcomm */
    BDASSERT(MPI_Send(msg, sizeof(bdmsg_t), MPI_BYTE, comm->wnranks[dest_node],
                      BDMPI_HDR_TAG, job->mpi_wcomm)
             == MPI_SUCCESS);

    if (msg->fnum == -1) { /* data is already in buf */
      BDASSERT(MPI_Send(buf, msize, MPI_BYTE, dest_node, msg->mpi_tag, comm->mpi_comm)
               == MPI_SUCCESS);
    }
    else { /* data resides on disk */
      BDASSERT(asprintf(&fname, "%s/%zd", job->jdesc->wdir, msg->fnum) != -1);
      BDASSERT((rsize = gk_getfsize(fname)) != -1);
      BDASSERT(rsize == msize);
      BDASSERT((fd = open(fname, O_RDONLY)) != -1);

      buf = gk_cmalloc(job->mmsize, "mstr_send: mpi buf");
      do {
        BDASSERT((rsize = read(fd, buf, job->mmsize)) > 0);
        BDASSERT(MPI_Send(buf, rsize, MPI_BYTE, dest_node, msg->mpi_tag, comm->mpi_comm)
                 == MPI_SUCCESS);
        msize -= rsize;
      } while (msize > 0);

      close(fd);
      unlink(fname);
      free(fname);
    }

    /* wait for the other node to tell you that it is done receiving */
    BDASSERT(MPI_Recv(&done, 1, MPI_INT, dest_node, msg->mpi_tag+1, comm->mpi_comm,
                      MPI_STATUS_IGNORE)
        == MPI_SUCCESS);

    gk_free((void **)&buf, LTERM);
  }

  /* Tell the slave that we are done. This is to enforce that the sends are
     delivered in the order that were sent. */
  bdmq_send(job->c2sMQs[srank], &done, sizeof(int));

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_send: source: %d, dest: %d [exiting]\n",
        job->mynode, msg->myrank, msg->source, msg->dest));

  BD_LET_RDLOCK(job->comms[orig_mcomm]->rwlock); /* unlock communicator */

  gk_free((void **)&arg, LTERM);

  return NULL;
}


/*************************************************************************/
/*! Response to a BDMPI_Send/BDMPI_Isend on a remote node.
    Uses MPI to receive the data from the source slave, saves it to
    memory/on disk, and creates a header entry for it.
*/
/*************************************************************************/
void *mstr_send_remote(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  bdmcomm_t *comm;
  int fd, source_node, count, done=1;
  char *buf=NULL;
  size_t rsize, msize;
  char *fname;
  MPI_Status status;

  /* incoming msg->mcomm actually stores mpi_commid, so we need to map it
     to this node's corresponding mcomm */
  msg->mcomm = babel_get_my_mcomm(job, msg->mcomm);

  BD_GET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_send_remote: source: %d, dest: %d [entering]\n",
      job->mynode, msg->myrank, msg->source, msg->dest));

  /* hook to the key info */
  comm = job->comms[msg->mcomm];

  BDASSERT(babel_is_local(comm, msg->dest));

  /* the node/master that is the source */
  source_node = babel_get_node(comm, msg->source);

  /* the size of the message */
  msize = bdmp_msize(msg->count, msg->datatype);

  /* lock it here in order to ensure message serializability. */
  pending_locksend(job, msg);

  if (msg->fnum == -1) { /* in memory message */
    buf = gk_cmalloc(msize, "mstr_send_remote: buf");
    BDASSERT(MPI_Recv(buf, msize, MPI_BYTE, source_node, msg->mpi_tag,
                 comm->mpi_comm, MPI_STATUS_IGNORE)
        == MPI_SUCCESS);
  }
  else { /* data resides on disk */
    BD_GET_LOCK(job->comm_lock);
    msg->fnum = xfer_getfnum();
    BD_LET_LOCK(job->comm_lock);

    BDASSERT(asprintf(&fname, "%s/%zd", job->jdesc->wdir, msg->fnum) != -1);
    BDASSERT((fd = open(fname, O_CREAT|O_WRONLY|O_TRUNC, S_IRUSR|S_IWUSR)) != -1);

    buf = gk_cmalloc(job->mmsize, "mstr_send: mpi buf");
    do {
      BDASSERT(MPI_Recv(buf, job->mmsize, MPI_BYTE, source_node, msg->mpi_tag,
                   comm->mpi_comm, &status)
          == MPI_SUCCESS);

      BDASSERT(MPI_Get_count(&status, MPI_BYTE, &count)
          == MPI_SUCCESS);
      rsize = count;
      BDASSERT(rsize == job->mmsize || rsize == msize);

      BDASSERT(gk_write(fd, buf, rsize) == rsize);
      msize -= rsize;
    } while (msize > 0);

    close(fd);
    free(fname);

    gk_free((void **)&buf, LTERM);
  }

#if 0
  /* Notify remote master that a receive request has been completed. */

  mmsg.mcomm = commid;
  mmsg.dest  = msg->source;
  mmsg.type  = BDMPI_MSGTYPE_RECVD;

  /* send the message header using the global node number of wcomm */
  BDASSERT(MPI_Send(mmsg, sizeof(bdmsg_t), MPI_BYTE,
                    comm->wnranks[source_node], BDMPI_HDR_TAG,
                    job->mpi_wcomm)
           == MPI_SUCCESS);
#endif

  pending_addsend(job, msg, buf, (buf == NULL ? 0 : msize));
  slvpool_munblock(job, babel_get_srank(comm, msg->dest));

  pending_unlocksend(job, msg);

  /* tell the sender that you are done with the receives */
  BDASSERT(MPI_Send(&done, 1, MPI_INT, source_node, msg->mpi_tag+1, comm->mpi_comm)
      == MPI_SUCCESS);

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mstr_send_remote: source: %d, dest: %d [exiting]\n",
      job->mynode, msg->myrank, msg->source, msg->dest));

  BD_LET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

  gk_free((void **)&arg, LTERM);

  return NULL;
}
