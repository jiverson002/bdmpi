/*!
\file
\brief A program to test Scatter 
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
  int i, k, myrank, npes, n=N;
  pid_t pid, ppid;
  int *sendbuf, *recvbuf;

  pid  = getpid();
  ppid = getppid();

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if (myrank == 0)
    printf("Testing scatter...\n");
  MPI_Barrier(MPI_COMM_WORLD);

  sendbuf = gk_imalloc(n*npes, "sendbuf");
  recvbuf = gk_imalloc(n*npes, "recvbuf");

  if (myrank == npes/2) {
    for (i=0; i<n*npes; i++)
      sendbuf[i] = i;
  }

  MPI_Scatter(sendbuf, n, MPI_INT, recvbuf, n, MPI_INT, npes/2, MPI_COMM_WORLD);

  for (i=0; i<n; i++) {
    if (recvbuf[i] != myrank*n+i)
      printf("[%3d]Error: recvbuf[%d]: got %d instead of %d\n", myrank, i, recvbuf[i], myrank*n+i);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    printf("Done.\n");

  gk_free((void **)&sendbuf, &recvbuf, LTERM);

  MPI_Finalize();

  return EXIT_SUCCESS;
}
