/*!
\file
\brief A virtual-memory based quicksort program for comparison purposes.
\date Started 4/23/2013
\author George
*/


#include <GKlib.h>
#include <bdmpi.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>


/**************************************************************************/
/**************************************************************************/
int main(int argc, char **argv)
{
  size_t i, n, npes;
  int *elmnts;
  double tmr=0.0;

  if (argc < 2) {
    fprintf(stderr, "Usage: %s npes nelems\n", argv[0]);
    return EXIT_FAILURE;
  }

  npes = atoi(argv[1]);
  n    = atoi(argv[2]) * npes;

  printf("n: %zu\n", n);

  srand(5);

  gk_startwctimer(tmr);
  elmnts = gk_imalloc(n, "elmnts");
  for (i=0; i<n; i++) 
    elmnts[i] = rand();

  gk_isorti(n, elmnts);

  for (i=1; i<n; i++) {
    if (elmnts[i-1] > elmnts[i])
      printf("Error at index %zu: %d %d\n", i, elmnts[i-1], elmnts[i]);
  }

  gk_stopwctimer(tmr);

  printf("Done: time: %8.3lf\n", gk_getwctimer(tmr));

  return EXIT_SUCCESS;
}



