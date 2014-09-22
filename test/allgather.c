/*!
\file
\brief A program to test Allgather 
\date Started 4/5/2013
\author George
*/

#include <GKlib.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <bdmpi.h>



/*************************************************************************/
/*! Let the game begin! */
/*************************************************************************/
int main(int argc, char *argv[])
{
  int i, j, k, myrank, npes, n=10000, trial;
  pid_t pid, ppid;
  char buf[16384];
  int *sendbuf, *recvbuf, *counts, *displs;

  pid  = getpid();
  ppid = getppid();

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  printf("[%d:%d] npes: %d, myrank: %d\n", (int)ppid, (int)pid, npes, myrank);
  fflush(stdout);


  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    printf("Allgather...\n");

  sendbuf = gk_imalloc(n, "sendbuf");
  recvbuf = gk_imalloc(n*npes, "recvbuf");

  for (trial=0; trial<10; trial++) {
    if (myrank == 0) {
      for (i=0; i<n; i++)
        sendbuf[i] = i;
    }
    MPI_Bcast(sendbuf, n, MPI_INT, 0, MPI_COMM_WORLD);

    /* test 1 */
    MPI_Allgather(sendbuf, n, MPI_INT, recvbuf, n, MPI_INT, MPI_COMM_WORLD);
  
    for (k=0; k<npes; k++) {
      for (i=0; i<n; i++) {
        if (recvbuf[k*n+i] != i)
          printf("[%3d]Error: recvbuf[%d]: got %d instead of %d\n", myrank, k*n+i, recvbuf[k*n+i], i);
      }
    }
  
    /* test 2 */
    for (i=0; i<n; i++)
      sendbuf[i] = i+myrank;
  
    MPI_Allgather(sendbuf, n, MPI_INT, recvbuf, n, MPI_INT, MPI_COMM_WORLD);
  
    for (k=0; k<npes; k++) {
      for (i=0; i<n; i++) {
        if (recvbuf[k*n+i] != i+k)
          printf("[%3d]Error: recvbuf[%d]: got %d instead of %d\n", myrank, k*n+i, recvbuf[k*n+i], i+k);
      }
    }
  }

  gk_free((void **)&sendbuf, &recvbuf, LTERM);
  if (myrank == 0)
    printf("Done.\n");

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    printf("Allgatherv...\n");
  fflush(stdout);

  sendbuf = gk_imalloc(n+myrank, "sendbuf");
  recvbuf = gk_imalloc((n+npes)*npes, "recvbuf");
  counts  = gk_imalloc(npes, "counts");
  displs  = gk_imalloc(npes, "displs");

  for (trial=0; trial<10; trial++) {
    gk_iset((n+npes)*npes, -1, recvbuf);
    for (i=0; i<n+myrank; i++)
      sendbuf[i] = i+myrank;
    for (i=0; i<npes; i++)
      counts[i] = n+i;
    for (i=0; i<npes; i++) {
      displs[i] = i*n;
      for (j=0; j<i; j++)
        displs[i] += j;
    }
  
    MPI_Allgatherv(sendbuf, n+myrank, MPI_INT, recvbuf, counts, displs, MPI_INT, 
        MPI_COMM_WORLD);
  
    for (k=0; k<npes; k++) {
      for (i=0; i<n+k; i++) {
        if (recvbuf[displs[k]+i] != i+k)
          printf("[%3d]Error: recvbuf[%d]: got %d instead of %d\n", 
              myrank, displs[k]+i, recvbuf[displs[k]+i], i+k);
      }
    }
  }

  gk_free((void **)&sendbuf, &recvbuf, &counts, &displs, LTERM);


  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    printf("Done.\n");

  MPI_Finalize();

  return EXIT_SUCCESS;
}
