/*!
\file
\brief A program to test Reduce/Allreduce 
\date Started 4/16/2013
\author George
*/

#include <GKlib.h>
#include <bdmpi.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>


#define N 50000

/*************************************************************************/
/*! Let the game begin! */
/*************************************************************************/
int main(int argc, char *argv[])
{
  int i, j, myrank, npes;
  pid_t pid, ppid;
  int sbuf[N], rbuf[N];
  MPI_Request request;

  setbuf(stdout, NULL);
  setbuf(stderr, NULL);

  pid  = getpid();
  ppid = getppid();

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  printf("[%d:%d] npes: %d, myrank: %d\n", (int)ppid, (int)pid, npes, myrank);

  srand(10*myrank);

  if (myrank == 0)
    printf("Starting reduce tests...\n");
  MPI_Barrier(MPI_COMM_WORLD);

  for (i=0; i<npes; i++) {
    for (j=0; j<N; j++)
      sbuf[j] = RandomInRange(128);
    printf("Root %4d: Reporter: %4d: buf: %5d %5d %5d\n", i, myrank, sbuf[0], sbuf[N/2], sbuf[N-1]);

    MPI_Reduce(sbuf, rbuf, N, MPI_INT, MPI_MIN, i, MPI_COMM_WORLD);

    if (myrank == i)
      printf("Root %4d: Reporter: %4d: minreduce: %5d %5d %5d\n", i, myrank, rbuf[0], rbuf[N/2], rbuf[N-1]);
  }

  for (i=0; i<npes; i++) {
    for (j=0; j<N; j++)
      sbuf[j] = RandomInRange(128);
    printf("Root %4d: Reporter: %4d: buf: %5d %5d %5d\n", i, myrank, sbuf[0], sbuf[N/2], sbuf[N-1]);

    MPI_Reduce(sbuf, rbuf, N, MPI_INT, MPI_MAX, i, MPI_COMM_WORLD);

    if (myrank == i)
      printf("Root %4d: Reporter: %4d: maxreduce: %5d %5d %5d\n", i, myrank, rbuf[0], rbuf[N/2], rbuf[N-1]);
  }

  for (i=0; i<npes; i++) {
    for (j=0; j<N; j++)
      sbuf[j] = RandomInRange(128);
    printf("Root %4d: Reporter: %4d: buf: %5d %5d %5d\n", i, myrank, sbuf[0], sbuf[N/2], sbuf[N-1]);

    MPI_Reduce(sbuf, rbuf, N, MPI_INT, MPI_SUM, i, MPI_COMM_WORLD);

    if (myrank == i)
      printf("Root %4d: Reporter: %4d: sumreduce: %5d %5d %5d\n", i, myrank, rbuf[0], rbuf[N/2], rbuf[N-1]);
  }


  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    printf("Starting allreduce tests...\n");
  MPI_Barrier(MPI_COMM_WORLD);

  for (j=0; j<N; j++)
    sbuf[j] = RandomInRange(128);
  printf("Reporter: %4d: buf: %5d %5d %5d\n", myrank, sbuf[0], sbuf[N/2], sbuf[N-1]);

  MPI_Allreduce(sbuf, rbuf, N, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  printf("Reporter: %4d: minreduce: %5d %5d %5d\n", myrank, rbuf[0], rbuf[N/2], rbuf[N-1]);
  //MPI_Barrier(MPI_COMM_WORLD);

  MPI_Allreduce(sbuf, rbuf, N, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  printf("Reporter: %4d: maxreduce: %5d %5d %5d\n", myrank, rbuf[0], rbuf[N/2], rbuf[N-1]);
  //MPI_Barrier(MPI_COMM_WORLD);

  MPI_Allreduce(sbuf, rbuf, N, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  printf("Reporter: %4d: sumreduce: %5d %5d %5d\n", myrank, rbuf[0], rbuf[N/2], rbuf[N-1]);

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    printf("Done.\n");

  MPI_Finalize();

  return EXIT_SUCCESS;
}
