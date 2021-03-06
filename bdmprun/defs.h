/*
 * defs.h
 *
 * This file contains various constant definitions
 *
 * Started 3/31/13
 * George
 *
 */

/* Constants for command-line arguments */
#define BDMPRUN_CMD_NS            1
#define BDMPRUN_CMD_NR            2
#define BDMPRUN_CMD_WDIR          3
#define BDMPRUN_CMD_SMSIZE        4
#define BDMPRUN_CMD_IMSIZE        5
#define BDMPRUN_CMD_MMSIZE        6
#define BDMPRUN_CMD_NOLOCKMEM     7
#define BDMPRUN_CMD_DBGLVL        100
#define BDMPRUN_CMD_SBMA          200
#define BDMPRUN_CMD_PGSIZE        201
#define BDMPRUN_CMD_RMSIZE        202
#define BDMPRUN_CMD_SAVALL        203
#define BDMPRUN_CMD_HELP          1000


/* Various defaults */
#define BDMPRUN_DEFAULT_NS            1
#define BDMPRUN_DEFAULT_NR            1
#define BDMPRUN_DEFAULT_NC            1
#define BDMPRUN_DEFAULT_WDIR          "/tmp/bdmpi"
#define BDMPRUN_DEFAULT_SBOPTS        BDMPI_SB_DISCARD
#define BDMPRUN_DEFAULT_IMSIZE        4
#define BDMPRUN_DEFAULT_SMSIZE        20
#define BDMPRUN_DEFAULT_MMSIZE        32
#define BDMPRUN_DEFAULT_PGSIZE        4
#define BDMPRUN_DEFAULT_RMSIZE        786432  /* 3072 MiB */
#define BDMPRUN_DEFAULT_SBMA          "none"
#define BDMPRUN_DEFAULT_LOCKMEM       1
#define BDMPRUN_DEFAULT_DBGLVL        0


/* Wakeup selection schemes */
#define BDMPRUN_WAKEUP_FIRST    1
#define BDMPRUN_WAKEUP_LAST     2
#define BDMPRUN_WAKEUP_FIFO     3
#define BDMPRUN_WAKEUP_LIFO     4
#define BDMPRUN_WAKEUP_VRSS     5
#define BDMPRUN_WAKEUP_PEND     6

/* MPI message tags (should be < 100) */
#define BDMPI_HDR_TAG            10  /* header for remote send */
#define BDMPI_ABORT_TAG          15  /* slvpool_abort request */
#define BDMPI_MRQ_TAG            20  /* requests to master node */
#define BDMPI_MRS_TAG            21  /* response from master node */

/* Various behavior controling constants */
#define BDMPI_MLOCK_BCAST        1
#define BDMPI_MLOCK_REDUCE       1
#define BDMPI_MLOCK_ALLGATHER    1
#define BDMPI_MLOCK_ALLTOALL     0
#define BDMPI_MLOCK_GATHER       0
#define BDMPI_MLOCK_SCATTER      0

/* Predefinied communicators */
#define BDMPI_COMM_WORLD         ((int)0)
#define BDMPI_COMM_NODE          ((int)1)
