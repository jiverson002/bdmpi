/*
 * bdscb.c
 *
 * The BDMP's functions for providing a shared communication buffer with 
 *  read/write syncronization
 *
 * Started 4/4/2013
 * George
 *
 */

#include "common.h"


/*************************************************************************/
/*! Creates the neccessary components of a shared communication buffer.

    \param size is the number of bytes to allocate (it should be a multiple
           of PAGESIZE
    \param tag is a string to be used as a prefix of the SMS
    \param num is used to uniquely identify the SMS been allocated 
           (eg., based on master or slave pid)
    \return the allocated shared communication buffer.
*/
/*************************************************************************/
bdscb_t *bdscb_create(size_t size, char *tag, int num)
{
  char name[256];
  bdscb_t *scb=NULL;

  scb = (bdscb_t *)bd_malloc(sizeof(bdscb_t), "bdscb_create: scb");

  /* allocate the shared memory, set its size, and map it */
  scb->size = size;

  sprintf(name, "/scbsm%s%d", tag, num);
  scb->smname = bd_strdup(name);
  shm_unlink(name);

  scb->fd = shm_open(name, O_CREAT|O_RDWR|O_TRUNC, S_IRUSR|S_IWUSR);
  if (scb->fd == -1)
    errexit("Failed on shm_open(scb->fd): %s\n", strerror(errno));

  if (ftruncate(scb->fd, size) == -1)
    errexit("Failed on ftruncate(scb->fd): %s\n", strerror(errno));
  
  scb->buf = mmap(NULL, size, PROT_READ|PROT_WRITE, MAP_SHARED, scb->fd, 0);
  if (scb->buf == MAP_FAILED)
    errexit("Failed on mmap(scb->buf): %s\n", strerror(errno));


  /* allocate the empty semaphore */
  sprintf(name, "/scbesem%s%d", tag, num);  
  scb->esemname = bd_strdup(name);
  sem_unlink(name);

  scb->esem = sem_open(name, O_CREAT, S_IRUSR|S_IWUSR, 1);
  if (scb->esem == SEM_FAILED)
    errexit("Failed on sem_open(scb->esem): %s\n", strerror(errno));


  /* allocate the full semaphore */
  sprintf(name, "/scbfsem%s%d", tag, num);  
  scb->fsemname = bd_strdup(name);
  sem_unlink(name);

  scb->fsem = sem_open(name, O_CREAT, S_IRUSR|S_IWUSR, 0);
  if (scb->fsem == SEM_FAILED)
    errexit("Failed on sem_open(scb->fsem): %s\n", strerror(errno));


  return scb;
}


/*************************************************************************/
/*! Opens the neccessary components of a shared communication buffer.
    Note that this call does not create the various IPC constructs and 
    it will fail if they have not already been created.

    \param size is the number of bytes to allocate (it should be a multiple
           of PAGESIZE
    \param tag is a string to be used as a prefix of the SMS
    \param num is used to uniquely identify the SMS been allocated 
           (eg., based on master or slave pid)
    \return the allocated shared communication buffer.
*/
/*************************************************************************/
bdscb_t *bdscb_open(size_t size, char *tag, int num)
{
  char name[256];
  bdscb_t *scb;

  scb = (bdscb_t *)bd_malloc(sizeof(bdscb_t), "bdscb_create: scb");

  /* allocate the shared memory object and map it */
  scb->size = size;

  sprintf(name, "/scbsm%s%d", tag, num);
  scb->smname = bd_strdup(name);

  scb->fd = shm_open(name, O_RDWR, S_IRUSR|S_IWUSR);
  if (scb->fd == -1) 
    errexit("Failed on shm_open(scb->fd): %s\n", strerror(errno));

  /*  I do not think this is needed.
  if (ftruncate(scb->fd, size) == -1)
    errexit("Failed on ftruncate(scb->fd): %s\n", strerror(errno));
  */

  scb->buf = mmap(NULL, size, PROT_READ|PROT_WRITE, MAP_SHARED, scb->fd, 0);
  if (scb->buf == MAP_FAILED)
    errexit("Failed on mmap(scb->buf): %s\n", strerror(errno));


  /* allocate the empty semaphore */
  sprintf(name, "/scbesem%s%d", tag, num);  
  scb->esemname = bd_strdup(name);

  scb->esem = sem_open(name, 0);
  if (scb->esem == SEM_FAILED)
    errexit("Failed on sem_open(scb->esem): %s\n", strerror(errno));


  /* allocate the full semaphore */
  sprintf(name, "/scbfsem%s%d", tag, num);  
  scb->fsemname = bd_strdup(name);

  scb->fsem = sem_open(name, 0);
  if (scb->fsem == SEM_FAILED)
    errexit("Failed on sem_open(scb->fsem): %s\n", strerror(errno));

  return scb;
}


/*************************************************************************/
/*! Closes the shared communication buffer, but it does not destroy it. 
    This should be called from the slaves.
*/
/*************************************************************************/
void bdscb_close(bdscb_t *scb)
{
  if (sem_close(scb->esem) == -1)
    errexit("Failed on sem_close(scb->esem): %s\n", strerror(errno));
  if (sem_close(scb->fsem) == -1)
    errexit("Failed on sem_close(scb->fsem): %s\n", strerror(errno));
  if (munmap(scb->buf, scb->size) == -1)
    errexit("Failed on munmap(scb->mem): %s\n", strerror(errno));
  if (close(scb->fd) == -1)
    errexit("Failed on close(scb->fd): %s\n", strerror(errno));

  bd_free((void **)&scb->smname, &scb->esemname, &scb->fsemname, &scb, LTERM);
}


/*************************************************************************/
/*! Closes and destroys the shared communication buffer. 
    This should be called from the master.
*/
/*************************************************************************/
void bdscb_destroy(bdscb_t *scb)
{
  if (sem_close(scb->esem) == -1)
    errexit("Failed on sem_close(scb->esem): %s\n", strerror(errno));
  if (sem_unlink(scb->esemname) == -1)
    errexit("Failed on sem_unlink(scb->esemname): %s\n", strerror(errno));

  if (sem_close(scb->fsem) == -1)
    errexit("Failed on sem_close(scb->fsem): %s\n", strerror(errno));
  if (sem_unlink(scb->fsemname) == -1)
    errexit("Failed on sem_unlink(scb->fsemname): %s\n", strerror(errno));

  if (munmap(scb->buf, scb->size) == -1)
    errexit("Failed on munmap(scb->mem): %s\n", strerror(errno));
  if (close(scb->fd) == -1)
    errexit("Failed on close(scb->fd): %s\n", strerror(errno));
  if (shm_unlink(scb->smname) == -1)
    errexit("Failed on shm_unlink(scb->smname): %s\n", strerror(errno));

  bd_free((void **)&scb->smname, &scb->esemname, &scb->fsemname, &scb, LTERM);
}


/*************************************************************************/
/*! Blocks until the buffer is empty */
/*************************************************************************/
int bdscb_wait_empty(bdscb_t *scb)
{
  int ret;
  for (;;) {
    ret = sem_wait(scb->esem);
    if (-1 == ret) {
      if (EINTR == errno)
        errno = 0;
      else
        return -1;
    }
    else {
      return 0;
    }
  }
  //return sem_wait(scb->esem);
}


/*************************************************************************/
/*! Blocks until the buffer is full */
/*************************************************************************/
int bdscb_wait_full(bdscb_t *scb)
{
  int ret;
  for (;;) {
    ret = sem_wait(scb->fsem);
    if (-1 == ret) {
      if (EINTR == errno)
        errno = 0;
      else
        return -1;
    }
    else {
      return 0;
    }
  }
  //return sem_wait(scb->fsem);
}


/*************************************************************************/
/*! Signals that a buffer is empty */
/*************************************************************************/
int bdscb_post_empty(bdscb_t *scb)
{
  return sem_post(scb->esem);
}


/*************************************************************************/
/*! Signals that a buffer is full */
/*************************************************************************/
int bdscb_post_full(bdscb_t *scb)
{
  return sem_post(scb->fsem);
}
