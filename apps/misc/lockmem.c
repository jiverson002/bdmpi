/*!
\file
\brief A simple program that allocates a locks a user specified number of bytes.
\date Started 10/1/2013
\author George
*/

#include <GKlib.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <sys/mman.h>



/*************************************************************************/
/*! Let the game begin! */
/*************************************************************************/
int main(int argc, char *argv[])
{
  int nbytes, sleeptime;
  char *buf;

  setbuf(stdout, NULL);
  setbuf(stderr, NULL);

  if (argc != 3) {
    printf("Usage %s nbytes sleeptime\n", argv[0]);
    exit(0);
  }

  nbytes    = atoi(argv[1]);
  sleeptime = atoi(argv[2]);

  printf("Allocating %d bytes\n", nbytes);

  buf = gk_cmalloc(nbytes, "buf");
  mlock(buf, nbytes);

  printf("Sleeping for %d seconds\n", sleeptime);
  sleep(sleeptime);

  free(buf);
}
