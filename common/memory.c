#ifndef _GNU_SOURCE
# define _GNU_SOURCE
#endif

#include <dlfcn.h>  /* dlsym */
#include <signal.h> /* raise */
#include <stdarg.h> /* va_* */
#include <stdint.h> /* int*_t */
#include <stdio.h>  /* *printf */
#include <string.h> /* memmove */
#include <unistd.h> /* ssize_t */


extern void * __libc_malloc(size_t);
extern void * __libc_free(void*);


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
  void *ptr=NULL;

  if (nbytes == 0)
    nbytes++;  /* Force mallocs to actually allocate some memory */

  ptr = (void *)__libc_malloc(nbytes);

  if (ptr == NULL) {
    fprintf(stderr, "***Memory allocation failed for %s. Requested size: %zu "
      "bytes", msg, nbytes);
    raise(SIGKILL);
    return NULL;
  }

  memset(ptr, 0, nbytes);

  return ptr;
}


/*************************************************************************/
/* This function is my wrapper around free, allows multiple pointers */
/*************************************************************************/
void bd_free(void **ptr1,...)
{
  va_list plist;
  void **ptr;

  if (*ptr1 != NULL) {
    __libc_free(*ptr1);
  }
  *ptr1 = NULL;

  va_start(plist, ptr1);
  while ((ptr = va_arg(plist, void **)) != (void**)0) {
    if (*ptr != NULL) {
      __libc_free(*ptr);
    }
    *ptr = NULL;
  }
  va_end(plist);
}


/************************************************************************/
/*! \brief Duplicates a string

This function is a replacement for C's standard <em>strdup()</em> function.
The key differences between the two are that gk_strdup():
  - uses the dynamic memory allocation routines of \e GKlib. 
  - it correctly handles NULL input strings.

The string that is returned must be freed by gk_free().

\param orgstr is the string that will be duplicated.
\returns A pointer to the newly created string.
\sa gk_free()
*/
/*************************************************************************/
char *bd_strdup(char *orgstr)
{
  int len;
  char *str=NULL;

  if (orgstr != NULL) {
    len = strlen(orgstr)+1;
    str = bd_malloc(len*sizeof(char), "gk_strdup: str");
    strcpy(str, orgstr);
  }

  return str;
}
