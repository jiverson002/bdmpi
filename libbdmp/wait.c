/*!
\file
\brief Implements the BDMPI_Wait function.
\date Started 4/14/2013
\author George
*/

#include "bdmplib.h"


/*************************************************************************/
/* Performs a wait operation */
/*************************************************************************/
int bdmp_Wait(sjob_t *job, BDMPI_Request *r_request, BDMPI_Status *status)
{
  BDMPI_Request request = *r_request;
  int ierror;

  S_IFSET(BDMPI_DBG_IPCS, 
      bdprintf("BDMPI_Wait: Testing a non-blocking request [goMQlen: %d]\n", bdmq_length(job->goMQ)));

  if (request == BDMPI_REQUEST_NULL) {
    return BDMPI_SUCCESS;
  }

  if (request->state != BDMPI_INPROGRESS) {
    if (status != BDMPI_STATUS_IGNORE)
      *status = request->status;
  }
  else { /* this happens only for BDMPI_REQUEST_IRECV */
    BDASSERT(request->type == BDMPI_REQUEST_IRECV);
    ierror = bdmp_Recv(job, request->buf, request->status.count, 
                  request->status.datatype, request->status.BDMPI_SOURCE, 
                  request->status.BDMPI_TAG, request->status.comm, status);
  }

  bd_free((void **)r_request, LTERM);
  *r_request = BDMPI_REQUEST_NULL;

  return BDMPI_SUCCESS;
}
