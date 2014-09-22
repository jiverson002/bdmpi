#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <bdmpi.h>

int main(int argc, char *argv[])
{
  int i, myrank, npes;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  for (i=0; i<npes; i++) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == i) {
      printf("[%5d:%5d] Hello from rank: %2d out of %2d\n", 
          (int)getppid(), (int)getpid(), myrank, npes);
      fflush(stdout);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();

  /* It is important to exit with a EXIT_SUCCESS status */
  return EXIT_SUCCESS; 
}

