/*!
\file
\brief Contails various functions for dealing with intra/inter communicators and 
       manipulations of their ranks.
\date Started 6/5/2013
\author George
*/

#include "bdmprun.h"



/*************************************************************************/
/*! Returns the slave rank of a global rank. It aborts if the global rank 
    does not correspond to a rank mapped to the local node. */
/*************************************************************************/
int babel_get_srank(bdmcomm_t *comm, int grank)
{
  if (grank < 0 || grank >= comm->gsize)
    slvpool_abort(1, "babel_get_srank: grank: %d is out of range [0:%d)\n", grank, comm->gsize);

  if (comm->g2nmap[grank] != comm->mynode)
    slvpool_abort(1, "babel_get_srank: grank: %d is not mapped on this node; mynode: %d, g2nmap: %d.\n",
        grank, comm->mynode, comm->g2nmap[grank]);

  return comm->sranks[comm->g2lmap[grank]];
}


/*************************************************************************/
/*! Returns the local rank of a global rank. It aborts if the global 
    rank does not correspond to a rank within the intra communicator. */
/*************************************************************************/
int babel_get_lrank(bdmcomm_t *comm, int grank)
{
  if (grank < 0 || grank >= comm->gsize)
    slvpool_abort(1, "babel_get_lrank: grank: %d is out of range [0:%d)\n", grank, comm->gsize);

  if (comm->g2nmap[grank] != comm->mynode)
    slvpool_abort(1, "babel_get_lrank: grank: %d is not mapped on this node; mynode: %d, g2nmap: %d.\n",
        grank, comm->mynode, comm->g2nmap[grank]);

  return comm->g2lmap[grank];
}


/*************************************************************************/
/*! Returns the global rank of a local rank. It aborts if the local 
    rank does not correspond to a rank within the intra communicator. */
/*************************************************************************/
int babel_get_grank(bdmcomm_t *comm, int lrank)
{
  if (lrank < 0 || lrank >= comm->lsize)
    slvpool_abort(1, "babel_get_grank: lrank: %d is out of range [0:%d)\n", lrank, comm->lsize);

  return comm->l2gmap[lrank];
}


/*************************************************************************/
/*! Returns the node number corresponding to the global rank. */
/*************************************************************************/
int babel_get_node(bdmcomm_t *comm, int grank)
{
  if (grank < 0 || grank >= comm->gsize)
    slvpool_abort(1, "babel_get_node: grank: %d is out of range [0:%d)\n", grank, comm->gsize);

  return comm->g2nmap[grank];
}


/*************************************************************************/
/*! Returns true if the global rank corresponds to a local slave. */
/*************************************************************************/
int babel_is_local(bdmcomm_t *comm, int grank)
{
  if (grank < 0 || grank >= comm->gsize)
    slvpool_abort(1, "babel_is_local: grank: %d is out of range [0:%d)\n", grank, comm->gsize);

  return (comm->g2nmap[grank] == comm->mynode);
}


/*************************************************************************/
/*! Returns the index into job->comms[] of the communicator with a given
    mpi_commid. */
/*************************************************************************/
int babel_get_my_mcomm(mjob_t *job, int mpi_commid)
{
  int i, mcomm=-1;

  if (mpi_commid == BDMPI_COMM_WORLD)
    return BDMPI_COMM_WORLD;  
  if (mpi_commid == BDMPI_COMM_NODE)
    return BDMPI_COMM_NODE;  
  
  BD_GET_LOCK(job->comm_lock);
  for (i=job->maxncomm-1; i>0; i--) {
    if (job->comms[i] != NULL && job->comms[i]->mpi_commid == mpi_commid) {
      mcomm = job->comms[i]->mcomm;
      break;
    }
  }
  BD_LET_LOCK(job->comm_lock);

  BDASSERT(mcomm != -1);

  return mcomm;
}
