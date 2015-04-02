#ifndef _GNU_SOURCE
# define _GNU_SOURCE
#endif

#include <dlfcn.h>  /* dlsym */
#include <signal.h> /* raise */
#include <stdarg.h> /* va_* */
#include <stdio.h>  /* *printf */
#include <stdint.h> /* int*_t */
#include <string.h> /* memmove */
#include <unistd.h> /* ssize_t */


/****************************************************************************/
/* need to provide temporary calloc function for dlsym */
/****************************************************************************/
#define HOOK_INIT(func)                                                     \
do {                                                                        \
  if (NULL == libc_##func)                                                  \
    *((void **) &libc_##func) = dlsym(RTLD_NEXT, #func);                    \
} while (0)


/*************************************************************************/
/*! This function is my wrapper around malloc that provides the following
    enhancements over malloc:
    * It always allocates one byte of memory, even if 0 bytes are requested.
      This is to ensure that checks of returned values do not lead to NULL
      due to 0 bytes requested.
    * It zeros-out the memory that is allocated. This is for a quick init
      of the underlying datastructures.
*/
/**************************************************************************/
void *bd_malloc(size_t nbytes, char *msg)
{
  static void* (*libc_malloc)(size_t)=NULL;

  void *ptr=NULL;

  HOOK_INIT(malloc);

  if (nbytes == 0)
    nbytes++;  /* Force mallocs to actually allocate some memory */

  ptr = (void *)libc_malloc(nbytes);

  if (ptr == NULL) {
    fprintf(stderr, "***Memory allocation failed for %s. Requested size: %zu "
      "bytes", msg, nbytes);
    raise(SIGKILL);
    return NULL;
  }

  return ptr;
}


/*************************************************************************
* This function is my wrapper around free, allows multiple pointers
**************************************************************************/
void bd_free(void **ptr1,...)
{
  static void (*libc_free)(void*)=NULL;
  va_list plist;
  void **ptr;

  HOOK_INIT(free);

  if (*ptr1 != NULL) {
    libc_free(*ptr1);
  }
  *ptr1 = NULL;

  va_start(plist, ptr1);
  while ((ptr = va_arg(plist, void **)) != (void**)0) {
    if (*ptr != NULL) {
      libc_free(*ptr);
    }
    *ptr = NULL;
  }
  va_end(plist);
}
