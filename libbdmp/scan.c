/*!
\file
\brief Implements the scan operation.
\date Started 9/2/2013
\author George/Shaden
*/

#include "bdmplib.h"


/*************************************************************************/
/*! Performs BDMPI_Scan() */
/*************************************************************************/
int bdmp_Scan(sjob_t *job, void *sendbuf, void *recvbuf, size_t count, 
        BDMPI_Datatype datatype, BDMPI_Op op, BDMPI_Comm comm)
{
  int i, k, npes, mype, d, partner, mask, tag, ierror=BDMPI_SUCCESS;
  size_t msize;
  char *abuf, *tbuf;

  /* some error checking */
  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }
  if (!datatype_isvalid(datatype)) {
    fprintf(stderr, "The sendtype is invalid.\n");
    return BDMPI_ERR_TYPE;
  }

  npes = comm->size;
  mype = comm->rank;
  
  tag = (++comm->copid)*BDMPL_COPID_MULT + BDMPL_ISCAN_TAG;

  msize = bdmp_msize(count, datatype);

  abuf = bd_malloc(msize, "abuf: bdmp_Scan");
  tbuf = bd_malloc(msize, "tbuf: bdmp_Scan");

  memcpy(abuf, sendbuf, msize);
  memcpy(recvbuf, sendbuf, msize);

  /* the # of dimensions of the enclosing hcube */
  d = gk_log2(2*npes-1);

  /* perform the d send/recv exchanges */
  for (mask=1, i=0; i<d; i++, mask=(mask<<1)) {
    partner = mype ^ mask;
    if (partner < npes) {
      if ((ierror = BDMPI_Sendrecv(abuf, msize, BDMPI_BYTE, partner, 
                        tag, tbuf, msize, BDMPI_BYTE, partner, tag, 
                        comm, BDMPI_STATUS_IGNORE)) != BDMPI_SUCCESS) 
        break; /* Sendrecv failed, jump out of loop */

      scan_op(abuf, tbuf, recvbuf, (partner < mype), count, datatype, op);
    }
  }

  bd_free((void**)&abuf, (void**)&tbuf, LTERM);

  return ierror;
}


/*************************************************************************/
/*! Performs BDMPI_Exscan() */
/*************************************************************************/
int bdmp_Exscan(sjob_t *job, void *sendbuf, void *recvbuf, size_t count, 
        BDMPI_Datatype datatype, BDMPI_Op op, BDMPI_Comm comm)
{
  int i, k, npes, mype, d, partner, mask, tag, ierror=BDMPI_SUCCESS;
  size_t msize;
  char *abuf, *tbuf;

  /* some error checking */
  if (comm == BDMPI_COMM_NULL) {
    fprintf(stderr, "Undefined communicator.\n");
    return BDMPI_ERR_COMM;
  }
  if (!datatype_isvalid(datatype)) {
    fprintf(stderr, "The sendtype is invalid.\n");
    return BDMPI_ERR_TYPE;
  }

  npes = comm->size;
  mype = comm->rank;
  
  tag = (++comm->copid)*BDMPL_COPID_MULT + BDMPL_ESCAN_TAG;

  msize = bdmp_msize(count, datatype);

  abuf = bd_malloc(msize, "abuf: bdmp_Scan");
  tbuf = bd_malloc(msize, "tbuf: bdmp_Scan");

  memcpy(abuf, sendbuf, msize);

  scan_init_op(recvbuf, count, datatype, op);

  /* the # of dimensions of the enclosing hcube */
  d = gk_log2(npes+npes-1);

  /* perform the d send/recv exchanges */
  for (mask=1, i=0; i<d; i++, mask=(mask<<1)) {
    partner = mype ^ mask;
    if (partner < npes) {
      if ((ierror = BDMPI_Sendrecv(abuf, msize, BDMPI_BYTE, partner, 
                        tag, tbuf, msize, BDMPI_BYTE, partner, tag, 
                        comm, BDMPI_STATUS_IGNORE)) != BDMPI_SUCCESS) 
        break; /* Sendrecv failed, jump out of loop */

      scan_op(abuf, tbuf, recvbuf, (partner < mype), count, datatype, op);
    }
  }

  bd_free((void**)&abuf, (void**)&tbuf, LTERM);

  return ierror;
}


