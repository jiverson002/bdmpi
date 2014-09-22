/*!
\file
\brief Implements the BDMPI_Test function.
\date Started 4/14/2013
\author George
*/

#include "bdmplib.h"


/*************************************************************************/
/* Performs a test operation */
/*************************************************************************/
int bdmp_Test(sjob_t *job, BDMPI_Request *r_request, int *flag, BDMPI_Status *status)
{
  int ierror;
  BDMPI_Request request = *r_request;
  BDMPI_Request nrequest = BDMPI_REQUEST_NULL;

  S_IFSET(BDMPI_DBG_IPCS, 
      bdprintf("BDMPI_Test: Testing a non-blocking request [goMQlen: %d]\n", bdmq_length(job->goMQ)));

  if (request->state == BDMPI_INPROGRESS) {
    BDASSERT(request->type == BDMPI_REQUEST_IRECV);

    ierror = bdmp_Irecv(job, request->buf, request->status.count, 
                 request->status.datatype, request->status.BDMPI_SOURCE, 
                 request->status.BDMPI_TAG, request->status.comm, &nrequest);

    gk_free((void **)r_request, LTERM);

    if (nrequest->state == BDMPI_SUCCESS) { /* async operation finished */
      if (status != BDMPI_STATUS_IGNORE) 
        *status = nrequest->status;

      gk_free((void **)&nrequest, LTERM);
      *r_request = BDMPI_REQUEST_NULL;
      *flag = 1;
    }
    else { /* async operation did not finish */
      *r_request = nrequest;
      *flag = 0;
    }
  }
  else {
    *flag = 1;

    if (status != BDMPI_STATUS_IGNORE) 
      *status  = request->status;
    
    gk_free((void **)r_request, LTERM);
    *r_request = BDMPI_REQUEST_NULL;
  }

  return BDMPI_SUCCESS;
}

