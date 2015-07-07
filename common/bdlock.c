/*
 * bdlock.c
 *
 * The BDMP's functions for providing inter-process locks via semaphores.
 *
 * Started 3/31/2013
 * George
 *
 */

#include "common.h"


/*************************************************************************/
/*! Creates the neccessary components of an IP lock used by BDMP.

    \param tag is a string to be used as a prefix of the SMS
    \param num is used to uniquely identify the SMS been allocated 
           (eg., based on master or slave pid)
    \return the created lock.
*/
/*************************************************************************/
bdlock_t *bdlock_create(char *tag, int num, int value)
{
  char name[256];
  bdlock_t *lock;

  lock = (bdlock_t *)bd_malloc(sizeof(bdlock_t), "bdlock_create: lock");

  sprintf(name, "/lk%s%d", tag, num);
  lock->name = bd_strdup(name);
  sem_unlink(name);

  lock->sem = sem_open(name, O_CREAT, S_IRUSR|S_IWUSR, value);
  if (lock->sem == SEM_FAILED)
    errexit("Failed on sem_open(lock->sem): %s\n", strerror(errno));

  return lock;
}


/*************************************************************************/
/*! Opens the neccessary components of an IP lock used by BDMP.

    \param tag is a string to be used as a prefix of the SMS
    \param num is used to uniquely identify the SMS been allocated 
           (eg., based on master or slave pid)
    \return the created lock.
*/
/*************************************************************************/
bdlock_t *bdlock_open(char *tag, int num)
{
  char name[256];
  bdlock_t *lock;

  lock = (bdlock_t *)bd_malloc(sizeof(bdlock_t), "bdlock_create: lock");

  sprintf(name, "/lk%s%d", tag, num);
  lock->name = bd_strdup(name);

  lock->sem = sem_open(name, 0);
  if (lock->sem == SEM_FAILED)
    errexit("Failed on sem_open(lock->sem): %s\n", strerror(errno));

  return lock;
}


/*************************************************************************/
/*! Closes an IP lock.
    This should be called from the slaves.
*/
/*************************************************************************/
void bdlock_close(bdlock_t *lock)
{
  if (sem_close(lock->sem) == -1)
    errexit("Failed on sem_close(lock->sem): %s\n", strerror(errno));

  bd_free((void **)&lock->name, &lock, LTERM);
}


/*************************************************************************/
/*! Closes and destroys an IP lock.
    This should be called from the master.
*/
/*************************************************************************/
void bdlock_destroy(bdlock_t *lock)
{
  if (sem_close(lock->sem) == -1)
    errexit("Failed on sem_close(lock->sem): %s\n", strerror(errno));
  if (sem_unlink(lock->name) == -1)
    errexit("Failed on sem_unlink(lock->name): %s\n", strerror(errno));

  bd_free((void **)&lock->name, &lock, LTERM);
}


/*************************************************************************/
/*! Acquires the lock */
/*************************************************************************/
int bdlock_lock(bdlock_t *lock)
{
  int ret;
  for (;;) {
    ret = sem_wait(lock->sem);
    if (-1 == ret && EINTR != errno)
      break;
    else if (-1 != ret)
      break;
  }
  return ret;
  //return sem_wait(lock->sem);
}


/*************************************************************************/
/*! Releases the lock */
/*************************************************************************/
int bdlock_unlock(bdlock_t *lock)
{
  return sem_post(lock->sem);
}
