/*!
\file
\brief A program to test Gather 
\date Started 4/20/2013
\author George
*/

#include <GKlib.h>
#include <bdmpi.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#define N 300000


/*************************************************************************/
/*! Let the game begin! */
/*************************************************************************/
int main(int argc, char *argv[])
{
  int i, k, myrank, npes, n;
  pid_t pid, ppid;
  int *sendbuf, *recvbuf;

  n = N;

  pid  = getpid();
  ppid = getppid();

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    printf("Testing gather...\n");

  sendbuf = gk_imalloc(n*npes, "sendbuf");
  recvbuf = gk_imalloc(n*npes, "recvbuf");

  for (i=0; i<n; i++)
    sendbuf[i] = i+myrank;

  MPI_Gather(sendbuf, n, MPI_INT, recvbuf, n, MPI_INT, npes/2, MPI_COMM_WORLD);

  if (myrank == npes/2) {
    for (k=0; k<npes; k++) {
      for (i=0; i<n; i++) {
        if (recvbuf[k*n+i] != i+k)
          printf("[%3d]Error: recvbuf[%d]: got %d instead of %d\n", myrank, k*n+i, recvbuf[k*n+i], i+k);
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    printf("Done.\n");

  gk_free((void **)&sendbuf, &recvbuf, LTERM);

  MPI_Finalize();

  return EXIT_SUCCESS;
}
