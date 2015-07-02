/*
 * bdipc.h
 *
 * Header info for bdmp's wrapping of unix IPC constructs 
 *
 * Started 3/31/13
 * George
 */

#ifndef _BDIPC_H_
#define _BDIPC_H_

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <semaphore.h>
#include <mqueue.h>


/*************************************************************************/
/*! This data structure is a wrapper for message queues */
/*************************************************************************/
typedef struct {
  char *name;           /*!< The name of the message queue */
  size_t msgsize;       /*!< The maximum size of the allowed message */
  void *buf;            /*!< The buffer for the messages */
  mqd_t mqdes;          /*!< The message queue descriptor */
} bdmq_t;


/*************************************************************************/
/*! This data structure provides inter-process locks via semaphores */
/*************************************************************************/
typedef struct {
  char *name;           /*!< The name of the semaphore */
  sem_t *sem;           /*!< The semaphore descriptor */
} bdlock_t;


/*************************************************************************/
/*! This data structure stores basic information for a shared memory mmap
    allocation */
/*************************************************************************/
typedef struct {
  char *smname;         /*!< The name of the shared memory region */
  size_t size;          /*!< The allocation size */
  size_t off;           /*!< Offset of first free word */
  int fd;               /*!< File descriptor of the shared memory region */
  void *mem;            /*!< Pointer to the start of the allocation */

  char *semname;        /*!< The name of the semaphore */
  sem_t *sem;           /*!< The semaphore controling access */
} bdsm_t;



/*************************************************************************/
/*! This data structure represents a shared communication buffer with
    empty/full synchronization constructs */
/*************************************************************************/
typedef struct {
  char *smname;         /*!< The name of the shared memory region */
  size_t size;          /*!< The allocation size */
  int fd;               /*!< File descriptor of the shared memory region */
  void *buf;            /*!< Pointer to the start of the allocation */

  char *esemname;       /*!< The name of the semaphore controlling access */
  sem_t *esem;          /*!< The semaphore for checking if the buffer is empty */
  char *fsemname;       /*!< The name of the semaphore controlling access */
  sem_t *fsem;          /*!< The semaphore for checking if the buffer is full */
} bdscb_t;



/*************************************************************************/
/* Prototypes */
/*************************************************************************/
/* bdscb.c */
bdscb_t *bdscb_create(size_t size, char *tag, int num);
bdscb_t *bdscb_open(size_t size, char *tag, int num);
void bdscb_close(bdscb_t *scb);
void bdscb_destroy(bdscb_t *scb);
int bdscb_wait_empty(bdscb_t *scb);
int bdscb_wait_full(bdscb_t *scb);
int bdscb_post_empty(bdscb_t *scb);
int bdscb_post_full(bdscb_t *scb);

/* bdsm.c */
bdsm_t *bdsm_create(size_t size, char *tag, int num);
bdsm_t *bdsm_open(size_t size, char *tag, int num);
void bdsm_close(bdsm_t *sm);
void bdsm_destroy(bdsm_t *sm);
int bdsm_lock(bdsm_t *sm);
int bdsm_unlock(bdsm_t *sm);
void *bdsm_malloc(bdsm_t *sm, size_t nbytes, char *msg);

/* bdmq.c */
bdmq_t *bdmq_create(char *tag, int num);
bdmq_t *bdmq_open(char *tag, int num);
void bdmq_close(bdmq_t *mq);
void bdmq_destroy(bdmq_t *mq);
ssize_t bdmq_send(bdmq_t *mq, void *buf, size_t size);
ssize_t bdmq_recv(bdmq_t *mq, void *buf, size_t size);
ssize_t bdmq_timedrecv(bdmq_t *mq, void *buf, size_t size, long dt);
void bdmq_register_function(bdmq_t *mq, void (*function)(union sigval));
int bdmq_length(bdmq_t *mq);

/* bdlock.c */
bdlock_t *bdlock_create(char *tag, int num, int value);
bdlock_t *bdlock_open(char *tag, int num);
void bdlock_close(bdlock_t *lock);
void bdlock_destroy(bdlock_t *lock);
int bdlock_lock(bdlock_t *lock);
int bdlock_unlock(bdlock_t *lock);


#endif 
