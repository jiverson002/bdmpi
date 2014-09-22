/*!
\file
\brief Various utility macros
\date Started 6/6/13
\author George
*/

#ifdef XXX
#define BD_GET_LOCK(lock)\
  BDASSERT(pthread_mutex_lock(lock) == 0)

#define BD_LET_LOCK(lock)\
  BDASSERT(pthread_mutex_unlock(lock) == 0)

#define BD_GET_RDLOCK(lock)\
  BDASSERT(pthread_rwlock_rdlock(lock) == 0)

#define BD_GET_WRLOCK(lock)\
  BDASSERT(pthread_rwlock_wrlock(lock) == 0)

#define BD_LET_RDLOCK(lock)\
  BDASSERT(pthread_rwlock_unlock(lock) == 0)

#define BD_LET_WRLOCK(lock)\
  BDASSERT(pthread_rwlock_unlock(lock) == 0)

#define BD_LET_RWLOCK(lock)\
  BDASSERT(pthread_rwlock_unlock(lock) == 0)

#endif
