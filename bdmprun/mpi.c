/*!
\file
\brief Various functions from converting between BDMP constructs to the
       corresponding MPI constructs.
\date Started 6/7/2013
\author George
*/


#include "bdmprun.h"


/*************************************************************************/
/*! Converts a BDMPI_Datatype into an MPI_Datatype */
/*************************************************************************/
MPI_Datatype mpi_dt(BDMPI_Datatype datatype)
{
  switch (datatype) {
    case BDMPI_CHAR:
      return MPI_CHAR;

    case BDMPI_SIGNED_CHAR:
      return MPI_SIGNED_CHAR;

    case BDMPI_UNSIGNED_CHAR:
      return MPI_UNSIGNED_CHAR;

    case BDMPI_BYTE:
      return MPI_BYTE;

    case BDMPI_WCHAR:
      return MPI_WCHAR;

    case BDMPI_SHORT:
      return MPI_SHORT;

    case BDMPI_UNSIGNED_SHORT:
      return MPI_UNSIGNED_SHORT;

    case BDMPI_INT:
      return MPI_INT;

    case BDMPI_UNSIGNED:
      return MPI_UNSIGNED;

    case BDMPI_LONG:
      return MPI_LONG;

    case BDMPI_UNSIGNED_LONG:
      return MPI_UNSIGNED_LONG;

    case BDMPI_LONG_LONG_INT:
      return MPI_LONG_LONG_INT;

    case BDMPI_UNSIGNED_LONG_LONG:
      return MPI_UNSIGNED_LONG_LONG;

    case BDMPI_INT8_T:
      return MPI_INT8_T;

    case BDMPI_UINT8_T:
      return MPI_UINT8_T;

    case BDMPI_INT16_T:
      return MPI_INT16_T;

    case BDMPI_UINT16_T:
      return MPI_UINT16_T;

    case BDMPI_INT32_T:
      return MPI_INT32_T;

    case BDMPI_UINT32_T:
      return MPI_UINT32_T;

    case BDMPI_INT64_T:
      return MPI_INT64_T;

    case BDMPI_UINT64_T:
      return MPI_UINT64_T;

    case BDMPI_SIZE_T:
      if (sizeof(size_t) == 4)
        return MPI_UINT32_T;
      else
        return MPI_UINT64_T;

    case BDMPI_SSIZE_T:
      if (sizeof(ssize_t) == 4)
        return MPI_INT32_T;
      else
        return MPI_INT64_T;

    case BDMPI_FLOAT:
      return MPI_FLOAT;

    case BDMPI_DOUBLE:
      return MPI_DOUBLE;

    case BDMPI_FLOAT_INT:
      return MPI_FLOAT_INT;

    case BDMPI_DOUBLE_INT:
      return MPI_DOUBLE_INT;

    case BDMPI_SHORT_INT:
      return MPI_SHORT_INT;

    case BDMPI_LONG_INT:
      return MPI_LONG_INT;

    case BDMPI_2INT:
      return MPI_2INT;

    default:
      slvpool_abort(1, "Undefined datatype: %d.\n", datatype);
  }

  return 0;
}


/*************************************************************************/
/*! Converts a BDMPI_Op into an MPI_Op */
/*************************************************************************/
MPI_Op mpi_op(BDMPI_Op op)
{
  switch (op) {
    case BDMPI_MAX:
      return MPI_MAX;

    case BDMPI_MIN:
      return MPI_MIN;

    case BDMPI_SUM:
      return MPI_SUM;

    case BDMPI_PROD:
      return MPI_PROD;

    case BDMPI_LAND:
      return MPI_LAND;

    case BDMPI_BAND:
      return MPI_BAND;

    case BDMPI_LOR:
      return MPI_LOR;

    case BDMPI_BOR:
      return MPI_BOR;

    case BDMPI_LXOR:
      return MPI_LXOR;

    case BDMPI_BXOR:
      return MPI_BXOR;

    case BDMPI_MAXLOC:
      return MPI_MAXLOC;

    case BDMPI_MINLOC:
      return MPI_MINLOC;

    default:
      slvpool_abort(1, "Undefined op.\n");
  }

  return 0;
}
