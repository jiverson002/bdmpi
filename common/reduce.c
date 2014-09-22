/*!
\file
\brief Implements the reduce operation.
\date Started 4/15/2013
\author George
*/

#include "common.h"

#define REDUCE_APPLY_OP1(i, a, b, count, op) \
  do {\
    switch (op) {\
      case BDMPI_MIN:\
        for (i=0; i<count; i++) \
          (a)[i] = ((a)[i]<(b)[i] ? (a)[i] : (b)[i]);\
        break;\
      case BDMPI_MAX:\
        for (i=0; i<count; i++) \
          (a)[i] = ((a)[i]>(b)[i] ? (a)[i] : (b)[i]);\
        break;\
      case BDMPI_SUM:\
        for (i=0; i<count; i++) \
          (a)[i] += (b)[i];\
        break;\
      case BDMPI_PROD:\
        for (i=0; i<count; i++) \
          (a)[i] *= (b)[i];\
        break;\
      case BDMPI_LAND:\
        for (i=0; i<count; i++) \
          (a)[i] = (a)[i] && (b)[i];\
        break;\
      case BDMPI_BAND:\
        for (i=0; i<count; i++) \
          (a)[i] &= (b)[i];\
        break;\
      case BDMPI_LOR:\
        for (i=0; i<count; i++) \
          (a)[i] = (a)[i] || (b)[i];\
        break;\
      case BDMPI_BOR:\
        for (i=0; i<count; i++) \
          (a)[i] |= (b)[i];\
        break;\
      case BDMPI_BXOR:\
        for (i=0; i<count; i++) \
          (a)[i] ^= (b)[i];\
        break;\
      default:\
        errexit("BDMPI_Op: %d not implemented\n", op);\
    }\
  } while (0)

#define REDUCE_APPLY_OP2(i, a, b, count, op) \
  do {\
    switch (op) {\
      case BDMPI_MIN:\
        for (i=0; i<count; i++) \
          (a)[i] = ((a)[i]<(b)[i] ? (a)[i] : (b)[i]);\
        break;\
      case BDMPI_MAX:\
        for (i=0; i<count; i++) \
          (a)[i] = ((a)[i]>(b)[i] ? (a)[i] : (b)[i]);\
        break;\
      case BDMPI_SUM:\
        for (i=0; i<count; i++) \
          (a)[i] += (b)[i];\
        break;\
      case BDMPI_PROD:\
        for (i=0; i<count; i++) \
          (a)[i] *= (b)[i];\
        break;\
      default:\
        errexit("BDMPI_Op: %d not implemented\n", op);\
    }\
  } while (0)

#define REDUCE_APPLY_OP3(i, a, b, count, op) \
  do {\
    switch (op) {\
      case BDMPI_MINLOC:\
        for (i=0; i<count; i++) \
          (a)[i] = ((a)[i].val<(b)[i].val ? (a)[i] : (b)[i]);\
        break;\
      case BDMPI_MAXLOC:\
        for (i=0; i<count; i++) \
          (a)[i] = ((a)[i].val>(b)[i].val ? (a)[i] : (b)[i]);\
        break;\
      default:\
        errexit("BDMPI_Op: %d not implemented\n", op);\
    }\
  } while (0)


/*************************************************************************/
/*! Performs the (a = a OP b) element-wise reduction operation.
 */
/*************************************************************************/
void reduce_op(void *a, void *b, size_t count, BDMPI_Datatype datatype, 
         BDMPI_Op op)
{
  size_t i;

  switch (datatype) {
    case BDMPI_CHAR:
    case BDMPI_BYTE:
      REDUCE_APPLY_OP1(i, (char *)a, (char *)b, count, op);
      break;
    case BDMPI_SIGNED_CHAR:
      REDUCE_APPLY_OP1(i, (signed char *)a, (signed char *)b, count, op);
      break;
    case BDMPI_UNSIGNED_CHAR:
      REDUCE_APPLY_OP1(i, (unsigned char *)a, (unsigned char *)b, count, op);
      break;
    case BDMPI_SHORT:
      REDUCE_APPLY_OP1(i, (short *)a, (short *)b, count, op);
      break;
    case BDMPI_UNSIGNED_SHORT:
      REDUCE_APPLY_OP1(i, (unsigned short *)a, (unsigned short *)b, count, op);
      break;
    case BDMPI_INT:
      REDUCE_APPLY_OP1(i, (int *)a, (int *)b, count, op);
      break;
    case BDMPI_UNSIGNED:
      REDUCE_APPLY_OP1(i, (unsigned int *)a, (unsigned int *)b, count, op);
      break;
    case BDMPI_LONG:
      REDUCE_APPLY_OP1(i, (long *)a, (long *)b, count, op);
      break;
    case BDMPI_UNSIGNED_LONG:
      REDUCE_APPLY_OP1(i, (unsigned long *)a, (unsigned long *)b, count, op);
      break;
    case BDMPI_LONG_LONG_INT:
      REDUCE_APPLY_OP1(i, (long long int *)a, (long long int *)b, count, op);
      break;
    case BDMPI_UNSIGNED_LONG_LONG:
      REDUCE_APPLY_OP1(i, (unsigned long long *)a, (unsigned long long *)b, count, op);
      break;
    case BDMPI_INT8_T:
      REDUCE_APPLY_OP1(i, (int8_t *)a, (int8_t *)b, count, op);
      break;
    case BDMPI_INT16_T:
      REDUCE_APPLY_OP1(i, (int16_t *)a, (int16_t *)b, count, op);
      break;
    case BDMPI_INT32_T:
      REDUCE_APPLY_OP1(i, (int32_t *)a, (int32_t *)b, count, op);
      break;
    case BDMPI_INT64_T:
      REDUCE_APPLY_OP1(i, (int64_t *)a, (int64_t *)b, count, op);
      break;
    case BDMPI_UINT8_T:
      REDUCE_APPLY_OP1(i, (uint8_t *)a, (uint8_t *)b, count, op);
      break;
    case BDMPI_UINT16_T:
      REDUCE_APPLY_OP1(i, (uint16_t *)a, (uint16_t *)b, count, op);
      break;
    case BDMPI_UINT32_T:
      REDUCE_APPLY_OP1(i, (uint32_t *)a, (uint32_t *)b, count, op);
      break;
    case BDMPI_UINT64_T:
      REDUCE_APPLY_OP1(i, (uint64_t *)a, (uint64_t *)b, count, op);
      break;
    case BDMPI_SIZE_T:
      REDUCE_APPLY_OP1(i, (size_t *)a, (size_t *)b, count, op);
      break;
    case BDMPI_SSIZE_T:
      REDUCE_APPLY_OP1(i, (ssize_t *)a, (ssize_t *)b, count, op);
      break;
    case BDMPI_FLOAT:
      REDUCE_APPLY_OP2(i, (float *)a, (float *)b, count, op);
      break;
    case BDMPI_DOUBLE:
      REDUCE_APPLY_OP2(i, (double *)a, (double *)b, count, op);
      break;
    case BDMPI_FLOAT_INT:
      REDUCE_APPLY_OP3(i, (bdvlp_fi_t *)a, (bdvlp_fi_t *)b, count, op);
      break;
    case BDMPI_DOUBLE_INT:
      REDUCE_APPLY_OP3(i, (bdvlp_di_t *)a, (bdvlp_di_t *)b, count, op);
      break;
    case BDMPI_LONG_INT:
      REDUCE_APPLY_OP3(i, (bdvlp_li_t *)a, (bdvlp_li_t *)b, count, op);
      break;
    case BDMPI_SHORT_INT:
      REDUCE_APPLY_OP3(i, (bdvlp_si_t *)a, (bdvlp_si_t *)b, count, op);
      break;
    case BDMPI_2INT:
      REDUCE_APPLY_OP3(i, (bdvlp_ii_t *)a, (bdvlp_ii_t *)b, count, op);
      break;
    default:
      errexit("Reduction has not implemented for datatype %d\n", datatype);
  }
}


/*************************************************************************/
/* Returns 0/1 if the specified op is supported */
/*************************************************************************/
int op_isvalid(BDMPI_Op op)
{
  switch (op) {
    case BDMPI_MAX:
    case BDMPI_MIN:          
    case BDMPI_SUM:
    case BDMPI_PROD:
    case BDMPI_LAND:
    case BDMPI_BAND:
    case BDMPI_LOR:
    case BDMPI_BOR:
    case BDMPI_LXOR:
    case BDMPI_BXOR:
    case BDMPI_MAXLOC:
    case BDMPI_MINLOC:
      return 1;
    default:
      return 0;
  }
}

