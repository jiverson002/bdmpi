/*
 * bdsm.c
 *
 * The BDMP's functions for providing shared memory functionality.
 *
 * Started 3/31/2013
 * George
 *
 */

#include "common.h"


/*************************************************************************/
/*! Creates the neccessary components of a shared memory (SMS) segment 
    used by BDMP.

    \param size is the number of bytes to allocate (it should be a multiple
           of PAGESIZE
    \param tag is a string to be used as a prefix of the SMS
    \param num is used to uniquely identify the SMS been allocated 
           (eg., based on master or slave pid)
    \return the allocated shared memory region and associated semaphore.
*/
/*************************************************************************/
bdsm_t *bdsm_create(size_t size, char *tag, int num)
{
  char name[256];
  bdsm_t *sm;

  sm = (bdsm_t *)bd_malloc(sizeof(bdsm_t), "bdsm_create: sm");

  /* allocate the shared memory, size it, and map it */
  sm->size = size;
  sm->off  = 0;

  sprintf(name, "/sm%s%d", tag, num);
  sm->smname = bd_strdup(name);
  shm_unlink(name);

  sm->fd = shm_open(name, O_CREAT|O_RDWR|O_TRUNC, S_IRUSR|S_IWUSR);
  if (sm->fd == -1)
    errexit("Failed on shm_open(sm->fd): %s\n", strerror(errno));
  
  if (ftruncate(sm->fd, size) == -1)
    errexit("Failed on ftruncate(sm->fd): %s\n", strerror(errno));

  sm->mem = mmap(NULL, size, PROT_READ|PROT_WRITE, MAP_SHARED, sm->fd, 0);
  if (sm->mem == MAP_FAILED)
    errexit("Failed on mmap(sm->mem): %s\n", strerror(errno));

  /* allocate the associated semaphore */
  sprintf(name, "/sem%s%d", tag, num);
  sm->semname = bd_strdup(name);
  sem_unlink(name);

  sm->sem = sem_open(name, O_CREAT, S_IRUSR|S_IWUSR, 1);
  if (sm->sem == SEM_FAILED)
    errexit("Failed on sem_open(sm->sem): %s\n", strerror(errno));

  return sm;
}


/*************************************************************************/
/*! Opens the neccessary components of a shared memory (SMS) segment used 
    by BDMP.

    \param size is the number of bytes to allocate (it should be a multiple
           of PAGESIZE
    \param tag is a string to be used as a prefix of the SMS
    \param num is used to uniquely identify the SMS been allocated 
           (eg., based on master or slave pid)
    \return the allocated shared memory region and associated semaphore.
*/
/*************************************************************************/
bdsm_t *bdsm_open(size_t size, char *tag, int num)
{
  char name[256];
  bdsm_t *sm;

  sm = (bdsm_t *)bd_malloc(sizeof(bdsm_t), "bdsm_create: sm");

  /* allocate the shared memory, size it, and map it */
  sm->size = size;
  sm->off  = 0;

  sprintf(name, "/sm%s%d", tag, num);
  sm->smname = bd_strdup(name);

  sm->fd = shm_open(name, O_RDWR, 0);
  if (sm->fd == -1)
    errexit("Failed on shm_open(sm->fd): %s\n", strerror(errno));

  /* I do not think this is needed
  if (ftruncate(sm->fd, size) == -1)
    errexit("Failed on ftruncate(sm->fd): %s\n", strerror(errno));
  */

  sm->mem = mmap(NULL, size, PROT_READ|PROT_WRITE, MAP_SHARED, sm->fd, 0);
  if (sm->mem == MAP_FAILED)
    errexit("Failed on mmap(sm->mem): %s\n", strerror(errno));

  /* allocate the associated semaphore */
  sprintf(name, "/sem%s%d", tag, num);
  sm->semname = bd_strdup(name);

  sm->sem = sem_open(name, 0);
  if (sm->sem == SEM_FAILED)
    errexit("Failed on sem_open(sm->sem): %s\n", strerror(errno));

  return sm;
}


/*************************************************************************/
/*! Closes the shared communication region.
    This should be called from the slaves.
*/
/*************************************************************************/
void bdsm_close(bdsm_t *sm)
{
  if (sem_close(sm->sem) == -1)
    errexit("Failed on sem_close(sm->sem): %s\n", strerror(errno));
  if (munmap(sm->mem, sm->size) == -1)
    errexit("Failed on munmap(sm->mem): %s\n", strerror(errno));
  if (close(sm->fd) == -1)
    errexit("Failed on close(sm->fd): %s\n", strerror(errno));

  bd_free((void **)&sm->smname, &sm->semname, &sm, LTERM);
}


/*************************************************************************/
/*! Closes and destroys the shared communication region.
    This should be called from the master.
*/
/*************************************************************************/
void bdsm_destroy(bdsm_t *sm)
{
  if (sem_close(sm->sem) == -1)
    errexit("Failed on sem_close(sm->sem): %s\n", strerror(errno));
  if (sem_unlink(sm->semname) == -1)
    errexit("Failed on sem_unlink(sm->semname): %s\n", strerror(errno));

  if (munmap(sm->mem, sm->size) == -1)
    errexit("Failed on munmap(sm->mem): %s\n", strerror(errno));
  if (close(sm->fd) == -1)
    errexit("Failed on close(sm->fd): %s\n", strerror(errno));
  if (shm_unlink(sm->smname) == -1)
    errexit("Failed on shm_unlink(sm->smname): %s\n", strerror(errno));

  bd_free((void **)&sm->smname, &sm->semname, &sm, LTERM);
}


/*************************************************************************/
/*! Acquires the access semaphore */
/*************************************************************************/
int bdsm_lock(bdsm_t *sm)
{
  return sem_wait(sm->sem);
}


/*************************************************************************/
/*! Releases the access semaphore */
/*************************************************************************/
int bdsm_unlock(bdsm_t *sm)
{
  return sem_post(sm->sem);
}


/*************************************************************************/
/*! Allocates memory from the SMR */
/*************************************************************************/
void *bdsm_malloc(bdsm_t *sm, size_t nbytes, char *msg)
{
  void *ptr = NULL;

  /* pad to make pointers 8-byte aligned */
  nbytes += (nbytes%8 == 0 ? 0 : 8 - nbytes%8);

  if (sm->off+nbytes < sm->size) {
    ptr = ((char *)sm->mem)+sm->off;
    sm->off += nbytes;
  }
  else {
    fprintf(stderr, "bdsm_malloc: could not allocate %zu bytes for [%s]\n", nbytes, msg);
  }

  return ptr;
}
