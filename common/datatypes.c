/*!
\file
\brief Implements various functions dealing with datatypes

\date Started 4/5/2013
\author George

\version\verbatim $Id: tch.c 6117 2012-09-25 18:14:47Z karypis $ \endverbatim
*/


#include "common.h"


/*************************************************************************/
/* Returns the sizeof the BDMPI_Datatypes */
/*************************************************************************/
size_t bdmp_sizeof(BDMPI_Datatype datatype)
{
  switch (datatype) {
      case BDMPI_CHAR:
      case BDMPI_SIGNED_CHAR:
      case BDMPI_UNSIGNED_CHAR:
      case BDMPI_BYTE:
        return sizeof(char);

      case BDMPI_WCHAR:
        return sizeof(wchar_t);

      case BDMPI_SHORT:
      case BDMPI_UNSIGNED_SHORT:
        return sizeof(short);

      case BDMPI_INT:
      case BDMPI_UNSIGNED:
        return sizeof(int);

      case BDMPI_LONG:
      case BDMPI_UNSIGNED_LONG:
        return sizeof(long);

      case BDMPI_LONG_LONG_INT:
      case BDMPI_UNSIGNED_LONG_LONG:
        return sizeof(long long int);

      case BDMPI_INT8_T:
      case BDMPI_UINT8_T:
        return sizeof(int8_t);

      case BDMPI_INT16_T:
      case BDMPI_UINT16_T:
        return sizeof(int16_t);

      case BDMPI_INT32_T:
      case BDMPI_UINT32_T:
        return sizeof(int32_t);

      case BDMPI_INT64_T:
      case BDMPI_UINT64_T:
        return sizeof(int64_t);

      case BDMPI_SIZE_T:
        return sizeof(size_t);

      case BDMPI_SSIZE_T:
        return sizeof(ssize_t);

      case BDMPI_FLOAT:
        return sizeof(float);

      case BDMPI_DOUBLE:
        return sizeof(double);

      case BDMPI_FLOAT_INT:
        return sizeof(bdvlp_fi_t);

      case BDMPI_DOUBLE_INT:
        return sizeof(bdvlp_di_t);

      case BDMPI_LONG_INT:
        return sizeof(bdvlp_li_t);

      case BDMPI_SHORT_INT:
        return sizeof(bdvlp_si_t);

      case BDMPI_2INT:
        return sizeof(bdvlp_ii_t);

      default:
        errexit("[%5d] +Undefined datatype: %d.\n", (int)getpid(), (int)datatype);
  }

  return 0;
}


/*************************************************************************/
/* Returns the total size of a message */
/*************************************************************************/
size_t bdmp_msize(size_t count, BDMPI_Datatype datatype)
{
  return count*bdmp_sizeof(datatype);
}


/*************************************************************************/
/* Returns 0/1 if the specified datatype is supported */
/*************************************************************************/
int datatype_isvalid(BDMPI_Datatype datatype)
{
  switch (datatype) {
    case BDMPI_CHAR:
    case BDMPI_SIGNED_CHAR:
    case BDMPI_UNSIGNED_CHAR:
    case BDMPI_BYTE:
    case BDMPI_WCHAR:
    case BDMPI_SHORT:
    case BDMPI_UNSIGNED_SHORT:
    case BDMPI_INT:
    case BDMPI_UNSIGNED:
    case BDMPI_LONG:
    case BDMPI_UNSIGNED_LONG:
    case BDMPI_LONG_LONG_INT:
    case BDMPI_UNSIGNED_LONG_LONG:
    case BDMPI_INT8_T:
    case BDMPI_UINT8_T:
    case BDMPI_INT16_T:
    case BDMPI_UINT16_T:
    case BDMPI_INT32_T:
    case BDMPI_UINT32_T:
    case BDMPI_INT64_T:
    case BDMPI_UINT64_T:
    case BDMPI_SIZE_T:
    case BDMPI_SSIZE_T:
    case BDMPI_FLOAT:
    case BDMPI_DOUBLE:
    case BDMPI_FLOAT_INT:
    case BDMPI_DOUBLE_INT:
    case BDMPI_LONG_INT:
    case BDMPI_SHORT_INT:
    case BDMPI_2INT:
      return 1;
    default:
      return 0;
  }
}
