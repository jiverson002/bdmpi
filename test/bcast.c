/*!
\file
\brief A program to test Bcast
\date Started 4/5/2013
\author George
*/

#include <GKlib.h>
#include <bdmpi.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>



/*************************************************************************/
/*! Let the game begin! */
/*************************************************************************/
int main(int argc, char *argv[])
{
  int i, k, myrank, npes, n=10000, trial;
  pid_t pid, ppid;
  char buf[16384];
  int *sendbuf, *recvbuf;

  pid  = getpid();
  ppid = getppid();

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  printf("[%d:%d] npes: %d, myrank: %d\n", (int)ppid, (int)pid, npes, myrank);

  if (myrank == 0)
    printf("Testing bcast...\n");
  MPI_Barrier(MPI_COMM_WORLD);

  for (i=0; i<npes; i++) {
    if (myrank == i) 
      sprintf(buf, "Root is %4d [%d]", i, pid);

    MPI_Bcast(buf, 128, MPI_CHAR, i, MPI_COMM_WORLD);

    if (myrank == (i+1)%npes)
      printf("Root %4d: Reporter: %4d: msg: %s\n", i, myrank, buf);
    if (myrank == (npes+i-1)%npes)
      printf("Root %4d: Reporter: %4d: msg: %s\n", i, myrank, buf);
  }

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    printf("Done.\n");

  MPI_Finalize();

  return EXIT_SUCCESS;
}
