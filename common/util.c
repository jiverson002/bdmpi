/*!
\file
\brief Implements various utility functions

\date Started 5/4/2013
\author George

\version\verbatim $Id: tch.c 6117 2012-09-25 18:14:47Z karypis $ \endverbatim
*/

#include "common.h"


/*************************************************************************/
/*! A more flexible version of madvise */
/*************************************************************************/
int bdmp_madvise(char *ptr, size_t size, int advise)
{
  ptrdiff_t addr = (ptrdiff_t)ptr, naddr;

  naddr = addr + 4095;
  naddr = (naddr>>12)<<12;

  size = size - (size_t)(naddr-addr);
  //printf("%p %p %zu\n", (void *)ptr, (void *)naddr, size);

  return madvise((void *)naddr, size, advise);
}


