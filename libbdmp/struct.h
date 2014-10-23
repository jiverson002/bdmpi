/*
 * struct.h
 *
 * This file contains data structures 
 *
 * Started 3/31/13
 * George
 *
 * $Id: struct.h 13900 2013-03-24 15:27:07Z karypis $
 */

#ifndef _STRUCT_H_
#define _STRUCT_H_


/*************************************************************************/
/*! This data structure stores all information associated with a bdmp 
    execution state and environment */
/*************************************************************************/
typedef struct {
  int rank;             /*!< my global rank */
  int lrank;            /*!< my global intra rank */

  /* process IDs */
  pid_t mypid;          /*!< The pid of the master */


  /* scheduling parameters */
  bdmq_t *goMQ;         /*!< Message queue to listed for go commands */

  /* shared memory objects */
  bdsm_t *globalSM;     /*!< The global SMR */
  bdscb_t *scb;         /*!< The communication buffer for master-slave communication */

  /* message queues for communication */
  bdmq_t *reqMQ;        /*!< Message queue for sending requests to the master */
  bdmq_t *c2sMQ;        /*!< Message queue for 1-1 communication master => slave */
  bdmq_t *c2mMQ;        /*!< Message queue for 1-1 communication slave => master */

  /* locks for IP synchronization */
  bdlock_t *mlockMX;    /*!< Mutex for mlock() system calls */
  bdlock_t *criticalMX; /*!< Mutex for critical sections */

  /* various other job-related constructs */
  bdjdesc_t *jdesc;     /*!< The job description to be allocated in the global 
                             shared memory */
  pid_t *spids;         /*!< The pids of the slave processes. This will
                             be allocated in the global SMR. */

  /* other */
  size_t smallmsg;      /*!< The size of a message in datatype to be considered
                             small and thus stored in memory */

} sjob_t;


#endif 
