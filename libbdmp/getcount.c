/*!
\file
\brief Implements the BDMPI_Get_count() function.
\date Started 4/12/2013
\author George
*/

#include "bdmplib.h"


/*************************************************************************/
/* Returns the counts of data items of type datatype received. 
   TODO: The implementation of this is not fully compliant. */
/*************************************************************************/
int bdmp_Get_count(sjob_t *job, BDMPI_Status *status, BDMPI_Datatype datatype, 
          size_t *count)
{
  *count = status->count;

  return BDMPI_SUCCESS;
}

