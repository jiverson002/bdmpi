/*!
\file
\brief Implements the computations associated with the scan operation.
\date Started 9/4/2013
\author George
*/

#include "common.h"

#define SCAN_INIT_OP1(i, type, a, count, op, min, max) \
  do {\
    switch (op) {\
      case BDMPI_MIN:\
        for (i=0; i<count; i++) ((type *)a)[i] = (max);\
        break;\
      case BDMPI_MAX:\
        for (i=0; i<count; i++) ((type *)a)[i] = (min);\
        break;\
      case BDMPI_SUM:\
      case BDMPI_LOR:\
      case BDMPI_BOR:\
        for (i=0; i<count; i++) ((type *)a)[i] = 0;\
        break;\
      case BDMPI_PROD:\
        for (i=0; i<count; i++) ((type *)a)[i] = 1;\
        break;\
      case BDMPI_LAND:\
      case BDMPI_BAND:\
      case BDMPI_BXOR:\
        for (i=0; i<count; i++) ((type *)a)[i] = (type)0xffffffffffffffff;\
        break;\
      default:\
        errexit("BDMPI_Op: %d not implemented\n", op);\
    }\
  } while (0)

#define SCAN_INIT_OP2(i, type, a, count, op, min, max) \
  do {\
    switch (op) {\
      case BDMPI_MIN:\
        for (i=0; i<count; i++) ((type *)a)[i] = (max);\
        break;\
      case BDMPI_MAX:\
        for (i=0; i<count; i++) ((type *)a)[i] = (min);\
        break;\
      case BDMPI_SUM:\
        for (i=0; i<count; i++) ((type *)a)[i] = 0;\
        break;\
      case BDMPI_PROD:\
        for (i=0; i<count; i++) ((type *)a)[i] = 1;\
        break;\
      default:\
        errexit("BDMPI_Op: %d not implemented\n", op);\
    }\
  } while (0)


#define SCAN_INIT_OP3(i, type, a, count, op, min, max) \
  do {\
    switch (op) {\
      case BDMPI_MINLOC:\
        for (i=0; i<count; i++) {\
          ((type *)a)[i].val = (max);\
          ((type *)a)[i].loc = -1;\
        }\
        break;\
      case BDMPI_MAXLOC:\
        for (i=0; i<count; i++) {\
          ((type *)a)[i].val = (min);\
          ((type *)a)[i].loc = -1;\
        }\
        break;\
      default:\
        errexit("BDMPI_Op: %d not implemented\n", op);\
    }\
  } while (0)


#define SCAN_APPLY_OP1(i, a, b, r, a2r, count, op) \
  do {\
    switch (op) {\
      case BDMPI_MIN:\
        for (i=0; i<count; i++) { \
          (a)[i] = ((a)[i]<(b)[i] ? (a)[i] : (b)[i]);\
          if (a2r) \
            (r)[i] = ((r)[i]<(b)[i] ? (r)[i] : (b)[i]);\
        }\
        break;\
      case BDMPI_MAX:\
        for (i=0; i<count; i++) { \
          (a)[i] = ((a)[i]>(b)[i] ? (a)[i] : (b)[i]);\
          if (a2r) \
            (r)[i] = ((r)[i]>(b)[i] ? (r)[i] : (b)[i]);\
        }\
        break;\
      case BDMPI_SUM:\
        for (i=0; i<count; i++) { \
          (a)[i] += (b)[i];\
          if (a2r) \
            (r)[i] += (b)[i];\
        }\
        break;\
      case BDMPI_PROD:\
        for (i=0; i<count; i++) { \
          (a)[i] *= (b)[i];\
          if (a2r) \
            (r)[i] *= (b)[i];\
        }\
        break;\
      case BDMPI_LAND:\
        for (i=0; i<count; i++) { \
          (a)[i] = (a)[i] && (b)[i];\
          if (a2r) \
            (r)[i] = (r)[i] && (b)[i];\
        }\
        break;\
      case BDMPI_BAND:\
        for (i=0; i<count; i++) { \
          (a)[i] &= (b)[i];\
          if (a2r) \
            (r)[i] &= (b)[i];\
        }\
        break;\
      case BDMPI_LOR:\
        for (i=0; i<count; i++) { \
          (a)[i] = (a)[i] || (b)[i];\
          if (a2r) \
            (r)[i] = (r)[i] || (b)[i];\
        }\
        break;\
      case BDMPI_BOR:\
        for (i=0; i<count; i++) { \
          (a)[i] |= (b)[i];\
          if (a2r) \
            (r)[i] |= (b)[i];\
        }\
        break;\
      case BDMPI_BXOR:\
        for (i=0; i<count; i++) { \
          (a)[i] ^= (b)[i];\
          if (a2r) \
            (r)[i] ^= (b)[i];\
        }\
        break;\
      default:\
        errexit("BDMPI_Op: %d not implemented\n", op);\
    }\
  } while (0)


#define SCAN_APPLY_OP2(i, a, b, r, a2r, count, op) \
  do {\
    switch (op) {\
      case BDMPI_MIN:\
        for (i=0; i<count; i++) { \
          (a)[i] = ((a)[i]<(b)[i] ? (a)[i] : (b)[i]);\
          if (a2r) \
            (r)[i] = ((r)[i]<(b)[i] ? (r)[i] : (b)[i]);\
        }\
        break;\
      case BDMPI_MAX:\
        for (i=0; i<count; i++) { \
          (a)[i] = ((a)[i]>(b)[i] ? (a)[i] : (b)[i]);\
          if (a2r) \
            (r)[i] = ((r)[i]>(b)[i] ? (r)[i] : (b)[i]);\
        }\
        break;\
      case BDMPI_SUM:\
        for (i=0; i<count; i++) { \
          (a)[i] += (b)[i];\
          if (a2r) \
            (r)[i] += (b)[i];\
        }\
        break;\
      case BDMPI_PROD:\
        for (i=0; i<count; i++) { \
          (a)[i] *= (b)[i];\
          if (a2r) \
            (r)[i] *= (b)[i];\
        }\
        break;\
      default:\
        errexit("BDMPI_Op: %d not implemented\n", op);\
    }\
  } while (0)


#define SCAN_APPLY_OP3(i, a, b, r, a2r, count, op) \
  do {\
    switch (op) {\
      case BDMPI_MINLOC:\
        for (i=0; i<count; i++) { \
          (a)[i] = ((a)[i].val<(b)[i].val ? (a)[i] : (b)[i]);\
          if (a2r) \
            (r)[i] = ((r)[i].val<(b)[i].val ? (r)[i] : (b)[i]);\
        }\
        break;\
      case BDMPI_MAXLOC:\
        for (i=0; i<count; i++) { \
          (a)[i] = ((a)[i].val>(b)[i].val ? (a)[i] : (b)[i]);\
          if (a2r) \
            (r)[i] = ((r)[i].val>(b)[i].val ? (r)[i] : (b)[i]);\
        }\
        break;\
      default:\
        errexit("BDMPI_Op: %d not implemented\n", op);\
    }\
  } while (0)



/*************************************************************************/
/*! Initializes the results array for a scan operation. */
/*************************************************************************/
void scan_init_op(void *a, size_t count, BDMPI_Datatype datatype, BDMPI_Op op)
{
  size_t i;

  switch (datatype) {
    case BDMPI_BYTE:
      SCAN_INIT_OP1(i, unsigned char, a, count, op, 0, UCHAR_MAX);
      break;

    case BDMPI_CHAR:
      SCAN_INIT_OP1(i, unsigned char, a, count, op, SCHAR_MIN, SCHAR_MAX);
      break;

    case BDMPI_SIGNED_CHAR:
      SCAN_INIT_OP1(i, signed char, a, count, op, SCHAR_MIN, SCHAR_MAX);
      break;

    case BDMPI_UNSIGNED_CHAR:
      SCAN_INIT_OP1(i, unsigned char, a, count, op, 0, UCHAR_MAX);
      break;

    case BDMPI_SHORT:
      SCAN_INIT_OP1(i, short, a, count, op, SHRT_MIN, SHRT_MAX);
      break;

    case BDMPI_UNSIGNED_SHORT:
      SCAN_INIT_OP1(i, unsigned short, a, count, op, 0, USHRT_MAX);
      break;

    case BDMPI_INT:
      SCAN_INIT_OP1(i, int, a, count, op, INT_MIN, INT_MAX);
      break;

    case BDMPI_UNSIGNED:
      SCAN_INIT_OP1(i, unsigned int, a, count, op, 0, UINT_MAX);
      break;

    case BDMPI_LONG:
      SCAN_INIT_OP1(i, long, a, count, op, LONG_MIN, LONG_MAX);
      break;

    case BDMPI_UNSIGNED_LONG:
      SCAN_INIT_OP1(i, unsigned long, a, count, op, 0, ULONG_MAX);
      break;

    case BDMPI_LONG_LONG_INT:
      SCAN_INIT_OP1(i, long long int, a, count, op, LLONG_MIN, LLONG_MAX);
      break;

    case BDMPI_UNSIGNED_LONG_LONG:
      SCAN_INIT_OP1(i, unsigned long long, a, count, op, 0, ULLONG_MAX);
      break;

    case BDMPI_INT8_T:
      SCAN_INIT_OP1(i, int8_t, a, count, op, INT8_MIN, INT8_MAX);
      break;

    case BDMPI_INT16_T:
      SCAN_INIT_OP1(i, int16_t, a, count, op, INT16_MIN, INT16_MAX);
      break;

    case BDMPI_INT32_T:
      SCAN_INIT_OP1(i, int32_t, a, count, op, INT32_MIN, INT32_MAX);
      break;

    case BDMPI_INT64_T:
      SCAN_INIT_OP1(i, int64_t, a, count, op, INT64_MIN, INT64_MAX);
      break;

    case BDMPI_UINT8_T:
      SCAN_INIT_OP1(i, uint8_t, a, count, op, 0, UINT8_MAX);
      break;

    case BDMPI_UINT16_T:
      SCAN_INIT_OP1(i, uint16_t, a, count, op, 0, UINT16_MAX);
      break;

    case BDMPI_UINT32_T:
      SCAN_INIT_OP1(i, uint32_t, a, count, op, 0, UINT32_MAX);
      break;

    case BDMPI_UINT64_T:
      SCAN_INIT_OP1(i, uint64_t, a, count, op, 0, UINT64_MAX);
      break;

    case BDMPI_SIZE_T:
      SCAN_INIT_OP1(i, size_t, a, count, op, 0, SIZE_MAX);
      break;

    case BDMPI_SSIZE_T:
      SCAN_INIT_OP1(i, ssize_t, a, count, op, INT64_MIN, INT64_MAX);
      break;

    case BDMPI_FLOAT:
      SCAN_INIT_OP2(i, float, a, count, op, FLT_MIN, FLT_MAX);
      break;

    case BDMPI_DOUBLE:
      SCAN_INIT_OP2(i, double, a, count, op, DBL_MIN, DBL_MAX);
      break;

    case BDMPI_FLOAT_INT:
      SCAN_INIT_OP3(i, bdvlp_fi_t, a, count, op, FLT_MIN, FLT_MAX);
      break;

    case BDMPI_DOUBLE_INT:
      SCAN_INIT_OP3(i, bdvlp_di_t, a, count, op, DBL_MIN, DBL_MAX);
      break;

    case BDMPI_LONG_INT:
      SCAN_INIT_OP3(i, bdvlp_li_t, a, count, op, LONG_MIN, LONG_MAX);
      break;

    case BDMPI_SHORT_INT:
      SCAN_INIT_OP3(i, bdvlp_si_t, a, count, op, SHRT_MIN, SHRT_MAX);
      break;

    case BDMPI_2INT:
      SCAN_INIT_OP3(i, bdvlp_ii_t, a, count, op, INT_MIN, INT_MAX);
      break;

    default:
      errexit("Scan has not implemented for datatype %d\n", datatype);
  }
}


/*************************************************************************/
/*! Performs the computations associated pairwise scan. */
/*************************************************************************/
void scan_op(void *a, void *b, void *r, int add2r, size_t count, 
         BDMPI_Datatype datatype, BDMPI_Op op)
{
  size_t i;

  switch (datatype) {
    case BDMPI_CHAR:
    case BDMPI_BYTE:
      SCAN_APPLY_OP1(i, (char *)a, (char *)b, (char *)r, 
          add2r, count, op);
      break;

    case BDMPI_SIGNED_CHAR:
      SCAN_APPLY_OP1(i, (signed char *)a, (signed char *)b, (signed char *)r, 
          add2r, count, op);
      break;

    case BDMPI_UNSIGNED_CHAR:
      SCAN_APPLY_OP1(i, (unsigned char *)a, (unsigned char *)b, (unsigned char *)r, 
          add2r, count, op);
      break;

    case BDMPI_SHORT:
      SCAN_APPLY_OP1(i, (short *)a, (short *)b, (short *)r, 
          add2r, count, op);
      break;

    case BDMPI_UNSIGNED_SHORT:
      SCAN_APPLY_OP1(i, (unsigned short *)a, (unsigned short *)b, (unsigned short *)r, 
          add2r, count, op);
      break;

    case BDMPI_INT:
      SCAN_APPLY_OP1(i, (int *)a, (int *)b, (int *)r, 
          add2r, count, op);
      break;

    case BDMPI_UNSIGNED:
      SCAN_APPLY_OP1(i, (unsigned int *)a, (unsigned int *)b, (unsigned int *)r, 
          add2r, count, op);
      break;

    case BDMPI_LONG:
      SCAN_APPLY_OP1(i, (long *)a, (long *)b, (long *)r, 
          add2r, count, op);
      break;

    case BDMPI_UNSIGNED_LONG:
      SCAN_APPLY_OP1(i, (unsigned long *)a, (unsigned long *)b, (unsigned long *)r, 
          add2r, count, op);
      break;

    case BDMPI_LONG_LONG_INT:
      SCAN_APPLY_OP1(i, (long long int *)a, (long long int *)b, (long long int *)r, 
          add2r, count, op);
      break;

    case BDMPI_UNSIGNED_LONG_LONG:
      SCAN_APPLY_OP1(i, (unsigned long long *)a, (unsigned long long *)b, (unsigned long long *)r, 
          add2r, count, op);
      break;

    case BDMPI_INT8_T:
      SCAN_APPLY_OP1(i, (int8_t *)a, (int8_t *)b, (int8_t *)r, 
          add2r, count, op);
      break;

    case BDMPI_INT16_T:
      SCAN_APPLY_OP1(i, (int16_t *)a, (int16_t *)b, (int16_t *)r, 
          add2r, count, op);
      break;

    case BDMPI_INT32_T:
      SCAN_APPLY_OP1(i, (int32_t *)a, (int32_t *)b, (int32_t *)r, 
          add2r, count, op);
      break;

    case BDMPI_INT64_T:
      SCAN_APPLY_OP1(i, (int64_t *)a, (int64_t *)b, (int64_t *)r, 
          add2r, count, op);
      break;

    case BDMPI_UINT8_T:
      SCAN_APPLY_OP1(i, (uint8_t *)a, (uint8_t *)b, (uint8_t *)r, 
          add2r, count, op);
      break;

    case BDMPI_UINT16_T:
      SCAN_APPLY_OP1(i, (uint16_t *)a, (uint16_t *)b, (uint16_t *)r, 
          add2r, count, op);
      break;

    case BDMPI_UINT32_T:
      SCAN_APPLY_OP1(i, (uint32_t *)a, (uint32_t *)b, (uint32_t *)r, 
          add2r, count, op);
      break;

    case BDMPI_UINT64_T:
      SCAN_APPLY_OP1(i, (uint64_t *)a, (uint64_t *)b, (uint64_t *)r, 
          add2r, count, op);
      break;

    case BDMPI_SIZE_T:
      SCAN_APPLY_OP1(i, (size_t *)a, (size_t *)b, (size_t *)r, 
          add2r, count, op);
      break;

    case BDMPI_SSIZE_T:
      SCAN_APPLY_OP1(i, (ssize_t *)a, (ssize_t *)b, (ssize_t *)r, 
          add2r, count, op);
      break;

    case BDMPI_FLOAT:
      SCAN_APPLY_OP2(i, (float *)a, (float *)b, (float *)r, 
          add2r, count, op);
      break;

    case BDMPI_DOUBLE:
      SCAN_APPLY_OP2(i, (double *)a, (double *)b, (double *)r, 
          add2r, count, op);
      break;

    case BDMPI_FLOAT_INT:
      SCAN_APPLY_OP3(i, (bdvlp_fi_t *)a, (bdvlp_fi_t *)b, (bdvlp_fi_t *)r, 
          add2r, count, op);
      break;

    case BDMPI_DOUBLE_INT:
      SCAN_APPLY_OP3(i, (bdvlp_di_t *)a, (bdvlp_di_t *)b, (bdvlp_di_t *)r, 
          add2r, count, op);
      break;

    case BDMPI_LONG_INT:
      SCAN_APPLY_OP3(i, (bdvlp_li_t *)a, (bdvlp_li_t *)b, (bdvlp_li_t *)r, 
          add2r, count, op);
      break;

    case BDMPI_SHORT_INT:
      SCAN_APPLY_OP3(i, (bdvlp_si_t *)a, (bdvlp_si_t *)b, (bdvlp_si_t *)r, 
          add2r, count, op);
      break;

    case BDMPI_2INT:
      SCAN_APPLY_OP3(i, (bdvlp_ii_t *)a, (bdvlp_ii_t *)b, (bdvlp_ii_t *)r, 
          add2r, count, op);
      break;

    default:
      errexit("Scan has not implemented for datatype %d\n", datatype);
  }
}


