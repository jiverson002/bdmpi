/*!
\file
\brief Implements various functions for reporting/debugging
\date Started 4/4/2013
\author George
*/

#include "common.h"


/*! \brief mutex controling output to stdout from the client/server */
static pthread_mutex_t printf_mtx = PTHREAD_MUTEX_INITIALIZER;


/*************************************************************************/
/*! This function uses a mutex to control concurrent writes */
/*************************************************************************/
void bdprintf(char *f_str,...)
{
  va_list argp;
  char hostname[9], string[16384];

  /* print it to the temporary string */
  va_start(argp, f_str);
  vsnprintf(string, 16383, f_str, argp);
  va_end(argp);
  string[16383] = '\0';

  /* output it to stderr */
  BD_GET_LOCK(&printf_mtx);
  gethostname(hostname, 9);
  fprintf(stderr, "[%8s:%6d]%s", hostname, (int)getpid(), string);
  fflush(stderr);
  BD_LET_LOCK(&printf_mtx);

}
