!     /* -*- Mode: Fortran; -*- */
!
!     /* 
!     \file bdmpif.h
!     \brief returning functions and constants
!
!     (C) 2013 Regents of the University of Minnesota
!
!     \date Started 11/25/13
!     \author Dominique
!     */
!


!     /*******************************************************************/
!     /* constants */
!     /*******************************************************************/

!     /* communicators */
      INTEGER BDMPI_COMM_WORLD
      PARAMETER (BDMPI_COMM_WORLD=0)
      INTEGER BDMPI_COMM_SELF
      PARAMETER (BDMPI_COMM_SELF=1)
      INTEGER BDMPI_COMM_CWORLD
      PARAMETER (BDMPI_COMM_CWORLD=2)
      INTEGER BDMPI_COMM_NODE
      PARAMETER (BDMPI_COMM_NODE=3)

!     /* types */
      INTEGER BDMPI_DATATYPE_NULL
      PARAMETER (BDMPI_DATATYPE_NULL=-1)
      INTEGER BDMPI_CHAR
      PARAMETER (BDMPI_CHAR=1001)
      INTEGER BDMPI_SIGNED_CHAR
      PARAMETER (BDMPI_SIGNED_CHAR=1002)
      INTEGER BDMPI_UNSIGNED_CHAR
      PARAMETER (BDMPI_UNSIGNED_CHAR=1003)
      INTEGER BDMPI_BYTE
      PARAMETER (BDMPI_BYTE=1004)
      INTEGER BDMPI_WCHAR
      PARAMETER (BDMPI_WCHAR=1005)
      INTEGER BDMPI_SHORT
      PARAMETER (BDMPI_SHORT=1006)
      INTEGER BDMPI_UNSIGNED_SHORT
      PARAMETER (BDMPI_UNSIGNED_SHORT=1007)
      INTEGER BDMPI_INT
      PARAMETER (BDMPI_INT=1008)
      INTEGER BDMPI_UNSIGNED
      PARAMETER (BDMPI_UNSIGNED=1009)
      INTEGER BDMPI_LONG
      PARAMETER (BDMPI_LONG=1010)
      INTEGER BDMPI_UNSIGNED_LONG
      PARAMETER (BDMPI_UNSIGNED_LONG=1011)
      INTEGER BDMPI_LONG_LONG_INT
      PARAMETER (BDMPI_LONG_LONG_INT=1012)
      INTEGER BDMPI_UNSIGNED_LONG_LONG
      PARAMETER (BDMPI_UNSIGNED_LONG_LONG=1013)
      INTEGER BDMPI_INT8_T
      PARAMETER (BDMPI_INT8_T=1014)
      INTEGER BDMPI_UINT8_T
      PARAMETER (BDMPI_UINT8_T=1015)
      INTEGER BDMPI_INT16_T
      PARAMETER (BDMPI_INT16_T=1016)
      INTEGER BDMPI_UINT16_T
      PARAMETER (BDMPI_UINT16_T=1017)
      INTEGER BDMPI_INT32_T
      PARAMETER (BDMPI_INT32_T=1018)
      INTEGER BDMPI_UINT32_T
      PARAMETER (BDMPI_UINT32_T=1019)
      INTEGER BDMPI_INT64_T
      PARAMETER (BDMPI_INT64_T=1020)
      INTEGER BDMPI_UINT64_T
      PARAMETER (BDMPI_UINT64_T=1021)
      INTEGER BDMPI_SIZE_T
      PARAMETER (BDMPI_SIZE_T=1022)
      INTEGER BDMPI_SSIZE_T
      PARAMETER (BDMPI_SSIZE_T=1023)
      INTEGER BDMPI_FLOAT
      PARAMETER (BDMPI_FLOAT=1024)
      INTEGER BDMPI_DOUBLE
      PARAMETER (BDMPI_DOUBLE=1025)
      INTEGER BDMPI_DOUBLE_PRECISION
      PARAMETER (BDMPI_DOUBLE_PRECISION=1025)
      INTEGER BDMPI_FLOAT_INT
      PARAMETER (BDMPI_FLOAT_INT=1026)
      INTEGER BDMPI_DOUBLE_INT
      PARAMETER (BDMPI_DOUBLE_INT=1027)
      INTEGER BDMPI_LONG_INT
      PARAMETER (BDMPI_LONG_INT=1028)
      INTEGER BDMPI_SHORT_INT
      PARAMETER (BDMPI_SHORT_INT=1029)
      INTEGER BDMPI_2INT
      PARAMETER (BDMPI_2INT=1030)

!     /* operations */
      INTEGER BDMPI_OP_NULL
      PARAMETER (BDMPI_OP_NULL=-1)
      INTEGER BDMPI_MAX
      PARAMETER (BDMPI_MAX=1001)
      INTEGER BDMPI_MIN
      PARAMETER (BDMPI_MIN=1002)
      INTEGER BDMPI_SUM
      PARAMETER (BDMPI_SUM=1003)
      INTEGER BDMPI_PROD
      PARAMETER (BDMPI_PROD=1004)
      INTEGER BDMPI_LAND
      PARAMETER (BDMPI_LAND=1005)
      INTEGER BDMPI_BAND
      PARAMETER (BDMPI_BAND=1006)
      INTEGER BDMPI_LOR
      PARAMETER (BDMPI_LOR=1007)
      INTEGER BDMPI_BOR
      PARAMETER (BDMPI_BOR=1008)
      INTEGER BDMPI_LXOR
      PARAMETER (BDMPI_LXOR=1009)
      INTEGER BDMPI_BXOR
      PARAMETER (BDMPI_BXOR=1010)
      INTEGER BDMPI_MAXLOC
      PARAMETER (BDMPI_MAXLOC=1011)
      INTEGER BDMPI_MINLOC
      PARAMETER (BDMPI_MINLOC=1012)

!     /* status */
      INTEGER BDMPI_STATUS_SIZE
      PARAMETER (BDMPI_STATUS_SIZE=48)
      INTEGER BDMPI_STATUS_IGNORE
      PARAMETER (BDMPI_STATUS_IGNORE=-1)

!     /* processes */
      INTEGER BDMPI_PROC_NULL
      PARAMETER (BDMPI_PROC_NULL=-1)

!     /* requests */
      INTEGER BDMPI_REQUEST_NULL
      PARAMETER (BDMPI_REQUEST_NULL=-1)

!     /* return codes */
      INTEGER BDMPI_SUCCESS
      PARAMETER (BDMPI_SUCCESS=1001)
      INTEGER BDMPI_INPROGRESS
      PARAMETER (BDMPI_INPROGRESS=1002)
      INTEGER BDMPI_ERR_UNKNOWN
      PARAMETER (BDMPI_ERR_UNKNOWN=1003)
      INTEGER BDMPI_ERR_COMM
      PARAMETER (BDMPI_ERR_COMM=1004)
      INTEGER BDMPI_ERR_ARG
      PARAMETER (BDMPI_ERR_ARG=1005)
      INTEGER BDMPI_ERR_BUFFER
      PARAMETER (BDMPI_ERR_BUFFER=1006)
      INTEGER BDMPI_ERR_COUNT
      PARAMETER (BDMPI_ERR_COUNT=1007)
      INTEGER BDMPI_ERR_TYPE
      PARAMETER (BDMPI_ERR_TYPE=1008)
      INTEGER BDMPI_ERR_TAG
      PARAMETER (BDMPI_ERR_TAG=1009)
      INTEGER BDMPI_ERR_ROOT
      PARAMETER (BDMPI_ERR_ROOT=1010)
      INTEGER BDMPI_ERR_RANK
      PARAMETER (BDMPI_ERR_RANK=1011)
      INTEGER BDMPI_ERR_TRUNCATE
      PARAMETER (BDMPI_ERR_TRUNCATE=1012)
      INTEGER BDMPI_ERR_OP
      PARAMETER (BDMPI_ERR_OP=1013)
      
!     /****************************************************************/
!     /* returning functions */
!     /****************************************************************/
      EXTERNAL BDMPI_WTIME
      REAL*8 BDMPI_WTIME

