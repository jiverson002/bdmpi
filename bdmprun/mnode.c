/*!
\file
\brief Various functions for the tasks associated with the master node.
\date Started 6/9/2013
\author George
*/


#include "bdmprun.h"



/*************************************************************************/
/*! Services any requests to the root master from the other masters. */
/*************************************************************************/
void mnode_service_xcomms(mjob_t *job)
{
  int flag;
  MPI_Message message;
  MPI_Status status;
  bdmsg_t msg;

  for (;;) {
    BDASSERT(MPI_Improbe(MPI_ANY_SOURCE, BDMPI_MRQ_TAG, job->mpi_wcomm, 
                 &flag, &message, &status)
        == MPI_SUCCESS);
  
    if (!flag) 
      break;
    
  
    BDASSERT(MPI_Mrecv(&msg, sizeof(bdmsg_t), MPI_BYTE, &message, &status)
        == MPI_SUCCESS);
  
    switch (msg.msgtype) {
      case BDMPI_MSGTYPE_CID:
        BD_GET_LOCK(job->comm_lock);
        msg.tag = job->next_mpi_commid++;
        BD_LET_LOCK(job->comm_lock);

        BDASSERT(MPI_Send(&msg, sizeof(bdmsg_t), MPI_BYTE, msg.myrank, 
                     BDMPI_MRS_TAG, job->mpi_wcomm) 
                 == MPI_SUCCESS); 
        break;
      
      default:
        slvpool_abort(1, "Got message: %d\n", msg.msgtype);
        return;
    }
  } 

  return;
}


/*************************************************************************/
/*! Queries the master node to get the next commid. */
/*************************************************************************/
int mnode_get_next_commid(mjob_t *job)
{
  int mpi_commid;
  bdmsg_t msg;
  MPI_Status status;

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mnode_get_next_commid: [entering]\n", 
        job->mynode));

  if (job->mynode == 0) {
    BD_GET_LOCK(job->comm_lock);
    mpi_commid = job->next_mpi_commid++;
    BD_LET_LOCK(job->comm_lock);
  }
  else {
    msg.msgtype = BDMPI_MSGTYPE_CID;
    msg.myrank  = job->mynode;

    /* send the request message to the master node */
    BDASSERT(MPI_Send(&msg, sizeof(bdmsg_t), MPI_BYTE, 0, BDMPI_MRQ_TAG, job->mpi_wcomm) 
             == MPI_SUCCESS); 

    /* wait for the response */
    BDASSERT(MPI_Recv(&msg, sizeof(bdmsg_t), MPI_BYTE, 0, BDMPI_MRS_TAG, job->mpi_wcomm,
                 &status)
        == MPI_SUCCESS);

    mpi_commid = msg.tag;
  }

  M_IFSET(BDMPI_DBG_IPCM, bdprintf("[MSTR%04d.%04d] mnode_get_next_commid: [exiting]\n", 
        job->mynode));

  return mpi_commid;
}
