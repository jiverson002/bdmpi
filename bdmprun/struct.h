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
/*! Data structure for storing a message header */
/*************************************************************************/
typedef struct header_s {
  bdmsg_t msg;            /*!< most of the message header info */
  char *buf;              /*!< pointer to the in-memory data */
  size_t len;             /*!< the size of the buffer */
  int mlck;               /*!< 0/1 if the memory in buf needs to be mlocked */
  int counter;            /*!< count-down counter for delivery */
  struct header_s *next;  /*!< the next header */
} header_t;


/*************************************************************************/
/*! A structure storing color/key/rank information for comm_split */
/*************************************************************************/
typedef struct {
  int color;            /*!< The color which determines the group membership */
  int key;              /*!< The key that determines the rank */
  int rank;             /*!< The initial rank */
} splitdata_t;


/*************************************************************************/
/*! The communicator on the master */
/*************************************************************************/
typedef struct {
  /* Inter-node component */
  MPI_Comm mpi_comm;  /*!< The MPI communicator */
  int nnodes;         /*!< The number of masters (MPI nodes) */
  int mynode;         /*!< The master's rank (MPI rank) */
  int *wnranks;       /*!< The ranks of the nodes in job->mpi_wcomm (i.e., global ranks) */
  int mpi_commid;     /*!< A unique communicator ID among the nodes involved */

  /* Intra-node component */
  int mcomm;            /*!< The master's ID of the communicator */
  int lsize;            /*!< The number of local processes in the communicator */
  int counter;          /*!< A counter for collective operations */
  int copid;            /*!< The next available collective operation ID */
  int *sranks;          /*!< The ranks of the slave processes within the node */

  /* Inter<=>Intra maps;
     grank: the BDMP rank of a process.
     lrank: the rank of a process within the slaves controlled by the master.
  */
  int gsize;            /*!< the total number of BDMP ranks */
  int *slvdist;         /*!< The distribution of the slaves */
  int *g2nmap;          /*!< map of BDMP ranks to nodes [nn*ns] */
  int *g2lmap;          /*!< map of BDMP ranks to lranks [nn*ns] */
  int *l2gmap;          /*!< map of lranks to BDMP ranks [ns] */

  /* Auxiliary data */
  pthread_rwlock_t *rwlock; /*!< Used to control concurrent access */
  splitdata_t *skeys;       /*!< The information storing color/key/rank for comm_split */
} bdmcomm_t;



/*************************************************************************/
/*! This data structure stores all information associated with a bdmp
    execution state and environment */
/*************************************************************************/
typedef struct {
  /* command-line parameters */
  int ns;               /*!< The number of slave processes to fork() */
  int nr_input;         /*!< The maximum number of running processes (-nr) */
  int nc;               /*!< The maximum number of slaves in a critical section */
  int sbopts;           /*!< The sb library options */
  size_t smsize;        /*!< The size of the shared comm buffer (-smsize*pagesize) */
  size_t imsize;        /*!< The maximum size of a message for in-memory buffering */
  size_t mmsize;        /*!< The maximum buffer size of inter-node p2p communication */
  size_t sbsize;        /*!< The minimum size for sbmalloc() */
  size_t pgsize;        /*!< The number of system pages per sb_malloc() page */
  size_t rmsize;        /*!< The maximum number of aggregate systems pages resident */
  int lockmem;          /*!< Specifies if the master will be locking its bcast/reduce
                             buffers */
  int dbglvl;           /*!< The dbglvl of the execution */
  char *iwdir;          /*!< The workding directory for storing various files (input) */
  char *exefile;        /*!< The name of the executable to be run */
  char **exeargv;       /*!< The command-line arguments of the executable */

  /* memory structure */
  size_t memrss;        /*!< Current resident set size for slaves on node */
  size_t memmax;        /*!< Maximum amount of memory available on system */
  size_t * slvrss;      /*!< Current resident set size for slaves on node */
  size_t * slvtot;      /*!< Total memory allocated for slaves on node */

  /* process IDs */
  pid_t mpid;           /*!< The pid of the master */

  /* scheduling parameters */
  bdmq_t **goMQs;       /*!< Message queues to signal slaves that they can run */

  /* shared memory objects */
  bdsm_t *globalSM;     /*!< The global shared memory region */

  /* shared communication buffers */
  bdscb_t **scbs;       /*!< The per slave shared communication buffer */

  /* message queues for communication */
  bdmq_t *reqMQ;        /*!< Slaves to master request queue */
  bdmq_t **c2sMQs;      /*!< Message queues for 1-1 communication master => slave */
  bdmq_t **c2mMQs;      /*!< Message queues for 1-1 communication slave => master */

  /* locks for IP synchronization */
  bdlock_t *mlockMX;    /*!< Mutex for mlock() system calls */
  bdlock_t *criticalMX; /*!< Mutex for critical sections */

  /* headers for the various pending messages */
  header_t **psends;          /*!< The per slave link-list of message headers for
                                   pending sends */
  pthread_mutex_t **plocks;   /*!< The locks for the above link-lists */


  /* various other job-related constructs */
  bdjdesc_t *jdesc;     /*!< The job description to be allocated in the global
                             shared memory */
  pid_t *spids;         /*!< The pids of the slave processes. This will
                             be allocated in the global SMR. */

  /* lists of running, runnable, and blocked slaves */
  int nalive, *alivelist, *alivemap;          /*!< The slave processes still existing */
  int nrunnable, *runnablelist, *runnablemap; /*!< The runnable slaves */
  int nrunning, *runninglist, *runningmap;    /*!< The slaves that are currently running */
  int nmblocked, *mblockedlist, *mblockedmap; /*!< The slaves that are currently blocked
                                                   waiting for messages */
  int ncblocked, *cblockedlist, *cblockedmap; /*!< The slaves that are currently blocked
                                                   due to collective operation */
  int *blockedts;                             /*!< Timestamp of when it was blocked */

  /* controls the initialization/termination of the slave pool */
  int njoined;          /*!< The number of slaves that have joined the pool */
  int nR;               /*!< The maximum number of running processes */


  /* communicators */
  int maxncomm;              /*!< The size of comms */
  int nfreecomm;             /*!< The number of communicators */
  int *commfreelist;         /*!< The list of free communicators */
  bdmcomm_t **comms;         /*!< An array of communicators */

  /* timers */
  double totalTmr;           /*!< total time spent by the master */
  double routeTmr;           /*!< time spent by the router */
  double sendTmr;            /*!< send comm timers */
  double recvTmr;            /*!< recv comm timers */
  double colTmr;             /*!< collective comm timers */
  double barrierTmr;         /*!< barrier timers */
  double commTmr;            /*!< communicator timers */
  double aux1Tmr;            /*!< aux1 timers */
  double aux2Tmr;            /*!< aux2 timers */
  double aux3Tmr;            /*!< aux3 timers */

  /* MPI info */
  MPI_Comm mpi_wcomm;   /*!< The world-communicator */
  int nnodes;           /*!< The number of MPI nodes */
  int mynode;           /*!< My rank within the MPI nodes */
  int *slvdist;         /*!< The node-based distribution of slave processes */
  int next_mpi_commid;  /*!< The next unique communicator ID */
  int next_mpi_tag;     /*!< The next tag for the message */

  /* pthread info */
  pthread_mutex_t *schedule_lock; /*!< Used for manipulating the scheduling structures */
  pthread_mutex_t *comm_lock;     /*!< Used for manipulating the communicator structures */
} mjob_t;


/*************************************************************************/
/*! Data structure for storing a pthread-function argument */
/*************************************************************************/
typedef struct {
  mjob_t *job;    /*!< the master's job info structure */
  bdmsg_t msg;    /*!< most of the message header info */
  char *buf;      /*!< pointer to the in-memory data */
} ptarg_t;


/*************************************************************************/
/*! This structure is used for storing the data of a merge */
/*************************************************************************/
typedef struct {
  size_t len;
  int *ids;
  void *buf;
} mergeinfo_t;

#endif
