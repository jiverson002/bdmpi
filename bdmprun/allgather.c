/*!
\file
\brief Various functions for allgather operations.
\date Started 4/6/2013
\author George
*/


#include "bdmprun.h"


/*************************************************************************/
/*! Response to a BDMPI_Allgather. - Send part
    Protocol:
       Blocks the process.
       Decreases the counter associated with the counter.
       Copies the data
       If counter becomes 0, then moves all processes to runnable state
       and sets the counter back to the size of that communicator.

*/
/*************************************************************************/
void *mstr_allgather_send(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg), nmsg;
  header_t *hdr;
  int i, srank, root_nrank, sleeping, fd;
  size_t rsize, msize;
  bdmcomm_t *comm;
  char *buf=NULL, *fname;

  BD_GET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  /* hook to the key info */
  comm  = job->comms[msg->mcomm];
  srank = babel_get_srank(comm, msg->myrank);

  /* deal with copid */
  if (comm->counter == comm->lsize)
    comm->copid++;
  if (bdmq_send(job->c2sMQs[srank], &(comm->copid), sizeof(int)) == -1)
    bdprintf("Failed to send a message to %d: %s\n", srank, strerror(errno));
  msg->copid = comm->copid;

  /* receive the specific info about this message */
  xfer_in_scb(job->scbs[srank], msg, sizeof(bdmsg_t), BDMPI_BYTE);

  if (msg->fnum == -1) {
    /* allocate memory and get the data */
    buf = gk_cmalloc(msg->count*bdmp_sizeof(msg->datatype), "allgather: buf");
    xfer_in_scb(job->scbs[srank], buf, msg->count, msg->datatype);
  }

  pending_addallgather(job, msg, buf,
      (buf == NULL ? 0 : bdmp_msize(msg->count, msg->datatype)));

  /* block the slave, once you get the ok from the slave */
  xfer_in_scb(job->scbs[srank], &sleeping, sizeof(int), BDMPI_BYTE);
  slvpool_cblock(job, srank);


  /* do some work, once all local slaves are done */
  if (--comm->counter == 0) {
    if (comm->nnodes > 1) { /* if non-local, get all masters involved */
      BD_LET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

      /* perform a sequence of broadcasts to get the data to all masters */
      for (i=0; i<comm->gsize; i++) {
        root_nrank = babel_get_node(comm, i);

        if (babel_is_local(comm, i)) {
          msg->source = i;
          if ((hdr = pending_getallgather(job, msg, NULL)) == NULL)
            slvpool_abort(1, "mstr_allgather_send: could not locate pending_allgather [%d %d %d]\n",
                msg->source, msg->dest, (int)msg->mcomm);

          nmsg = hdr->msg; /* copy the info into 'msg' */

          /* send information about the message to everybody */
          BDASSERT(MPI_Bcast(&nmsg, sizeof(bdmsg_t), MPI_BYTE, root_nrank, comm->mpi_comm)
                   == MPI_SUCCESS);

          msize = bdmp_msize(nmsg.count, nmsg.datatype);

          /* send the actual data */
          if (nmsg.fnum == -1) {  /* in-memory storage */
            BDASSERT(MPI_Bcast(hdr->buf, msize, MPI_BYTE, root_nrank, comm->mpi_comm)
                     == MPI_SUCCESS);
          }
          else { /* disk-based storage */
            BDASSERT(asprintf(&fname, "%s/%zd", job->jdesc->wdir, nmsg.fnum) != -1);
            BDASSERT((rsize = gk_getfsize(fname)) != -1);
            BDASSERT(rsize == msize);
            BDASSERT((fd = open(fname, O_RDONLY)) != -1);

            buf = gk_cmalloc(job->mmsize, "mstr_allgather_send: mpi buf");
            do {
              BDASSERT((rsize = read(fd, buf, job->mmsize)) > 0);
              BDASSERT(MPI_Bcast(buf, rsize, MPI_BYTE, root_nrank, comm->mpi_comm)
                       == MPI_SUCCESS);
              msize -= rsize;
            } while (msize > 0);

            free(fname);
            gk_free((void **)&buf, LTERM);
          }
        }
        else {
          /* everybody receives the information about the message */
          BDASSERT(MPI_Bcast(&nmsg, sizeof(bdmsg_t), MPI_BYTE, root_nrank, comm->mpi_comm)
                   == MPI_SUCCESS);

          msize = bdmp_msize(nmsg.count, nmsg.datatype);

          nmsg.mcomm = msg->mcomm;  /* switch it to the locally known # ('msg' was
                                       a message from that master about the same
                                       collective operation */

          /* get the actual data */
          if (nmsg.fnum == -1) {  /* in-memory storage */
            /* allocate memory and get the data */
            buf = gk_cmalloc(msize, "allgather: buf");

            BDASSERT(MPI_Bcast(buf, msize, MPI_BYTE, root_nrank, comm->mpi_comm)
                     == MPI_SUCCESS);
          }
          else { /* disk-based storage */
            BD_GET_LOCK(job->comm_lock);
            nmsg.fnum = xfer_getfnum();
            BD_LET_LOCK(job->comm_lock);

            BDASSERT(asprintf(&fname, "%s/%zd", job->jdesc->wdir, nmsg.fnum) != -1);
            BDASSERT((fd = open(fname, O_CREAT|O_WRONLY|O_TRUNC, S_IRUSR|S_IWUSR)) != -1);

            buf = gk_cmalloc(job->mmsize, "mstr_allgather_send: mpi buf");
            do {
              rsize = (msize > job->mmsize ? job->mmsize : msize);
              BDASSERT(MPI_Bcast(buf, rsize, MPI_BYTE, root_nrank, comm->mpi_comm)
                       == MPI_SUCCESS);
              BDASSERT(gk_write(fd, buf, rsize) == rsize);
              msize -= rsize;
            } while (msize > 0);

            close(fd);
            free(fname);

            gk_free((void **)&buf, LTERM);
          }

          pending_addallgather(job, &nmsg, buf,
              (buf == NULL ? 0 : bdmp_msize(nmsg.count, nmsg.datatype)));
        }
      }

      BD_GET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */
    }

    /* all processes have called the allgather, unblock them in order to get to
       the second step */
    comm->counter = comm->lsize;
    for (i=0; i<comm->lsize; i++)
      slvpool_cunblock(job, comm->sranks[i]);
  }

  BD_LET_WRLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

  gk_free((void **)&arg, LTERM);

  return NULL;
}


/*************************************************************************/
/*! Response to a BDMPI_Allgather - Recv part.
    Protocol:
      - Finds the headers of the allgathers
      - Reads and copies the data to the slave.
    Note:
      - By construction the header of the allgather should be there! If not,
        then something went very bad.
*/
/*************************************************************************/
void *mstr_allgather_recv(void *arg)
{
  mjob_t *job = ((ptarg_t *)arg)->job;
  bdmsg_t *msg = &(((ptarg_t *)arg)->msg);
  int k, response, srank, counter;
  header_t *hdr;
  bdmcomm_t *comm;
  ssize_t fnum;

  BD_GET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* lock communicator */

  /* hook to the key info */
  comm  = job->comms[msg->mcomm];
  srank = babel_get_srank(comm, msg->myrank);

  /* go and copy the data received from all the slaves to that slave */
  for (k=0; k<comm->gsize; k++) {
    /* get the header */
    msg->source = k;
    if ((hdr = pending_getallgather(job, msg, &counter)) == NULL)
      slvpool_abort(1, "mstr_allgather_recv: could not locate pending_allgather [%d %d %d]\n",
          msg->source, msg->dest, (int)msg->mcomm);

    /* send the header first */
    xfer_out_scb(job->scbs[srank], &hdr->msg, sizeof(bdmsg_t), BDMPI_BYTE);

    /* send the actual data to slave */
    if ((fnum = hdr->msg.fnum) == -1)
      xfer_out_scb(job->scbs[srank], hdr->buf, hdr->msg.count, hdr->msg.datatype);

    /* wait for the ACK from the slave */
    bdmq_recv(job->c2mMQs[srank], &response, sizeof(int));

    if (counter == 0) { /* switch from RD to WR locks */
      BD_LET_RDLOCK(job->comms[msg->mcomm]->rwlock);
      BD_GET_WRLOCK(job->comms[msg->mcomm]->rwlock);

      pending_freeheader(job, &hdr);

      if (fnum != -1)
        xfer_unlink(fnum);

      BD_LET_WRLOCK(job->comms[msg->mcomm]->rwlock);
      BD_GET_RDLOCK(job->comms[msg->mcomm]->rwlock);
    }
  }

  BD_LET_RDLOCK(job->comms[msg->mcomm]->rwlock); /* unlock communicator */

  gk_free((void **)&arg, LTERM);

  return NULL;
}
