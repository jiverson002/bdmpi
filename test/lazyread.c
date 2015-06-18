#include <GKlib.h>
#include <bdmpi.h>
#include <stdio.h>
#include <stdlib.h>

int
main(int argc, char* argv[])
{
  int i, n, np, rank, pagesize;
  int ret=0, gret=0;
  long sum, gsum;
  MPI_Request req;
  MPI_Status status;
  int * arr=NULL;
  double tmr0, tmr1, tmr2, cur, max;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  gk_clearwctimer(tmr0);
  gk_clearwctimer(tmr1);
  gk_clearwctimer(tmr2);

  n = 1024*1024*16;
  pagesize = sysconf(_SC_PAGESIZE);

  if (NULL == (arr=(int*)malloc(n*sizeof(int)))) {
    ret = -1;
    goto DONE;
  }

  /* populate every entry in array */
  gk_startwctimer(tmr0);
  for (i=0; i<n; ++i)
    arr[i] = rank+1;
  gk_stopwctimer(tmr0);

  /* force slaves to write memory */
  MPI_Barrier(MPI_COMM_WORLD);

  /* do some work reading one entry from every page */
  gk_startwctimer(tmr1);
  for (sum=0,i=0; i<n/pagesize; ++i)
    sum += arr[i*pagesize];
  gk_stopwctimer(tmr1);

  /* force slaves to write memory */
  MPI_Barrier(MPI_COMM_WORLD);

  /* do some work reading the same number of entires from consecutive pages */
  gk_startwctimer(tmr2);
  for (i=0; i<n/pagesize; ++i)
    sum += arr[i];
  gk_stopwctimer(tmr2);

  /* force slaves to write memory */
  MPI_Barrier(MPI_COMM_WORLD);

  for (i=0; i<n; ++i) {
    if (rank+1 != arr[i]) {
      printf("FAIL.1\n");
      goto FAIL;
    }
  }

  free(arr);

DONE:
  MPI_Reduce(&sum, &gsum, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  if (0 == rank && 2*(n/pagesize)*(rank+1) != gsum)
    printf("FAIL.2\n");

  cur = gk_getwctimer(tmr0);
  MPI_Reduce(&cur, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (0 == rank && 0.0 < max)
    printf(" work0Tmr:  %10.4lf\n", max);

  cur = gk_getwctimer(tmr1);
  MPI_Reduce(&cur, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (0 == rank && 0.0 < max)
    printf(" work1Tmr:  %10.4lf\n", max);

  cur = gk_getwctimer(tmr2);
  MPI_Reduce(&cur, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (0 == rank && 0.0 < max)
    printf(" work2Tmr:  %10.4lf\n", max);

  MPI_Reduce(&ret, &gret, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
  if (0 == rank)
    printf(" return: %d\n", gret);

FAIL:
  MPI_Finalize();

  return EXIT_SUCCESS;
}
