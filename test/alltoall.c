/*!
\file
\brief A program to test Alltoall
\date Started 4/20/2013
\author George
*/

#include <GKlib.h>
#include <bdmpi.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#define N 5000


/*************************************************************************/
/*! Let the game begin! */
/*************************************************************************/
int main(int argc, char *argv[])
{
  int iter, i, k, myrank, npes, n;
  pid_t pid, ppid;
  int *sendbuf, *recvbuf;
  size_t *sendcounts, *recvcounts, *sdispls, *rdispls;

  setbuf(stdout, NULL);
  setbuf(stderr, NULL);

  n = N;

  pid  = getpid();
  ppid = getppid();

  BDMPI_Init(&argc, &argv);

  BDMPI_Comm_size(BDMPI_COMM_WORLD, &npes);
  BDMPI_Comm_rank(BDMPI_COMM_WORLD, &myrank);


  BDMPI_Barrier(BDMPI_COMM_WORLD);
  if (myrank == 0)
    printf("Testing alltoall...\n");

  sendbuf = gk_imalloc(n*npes, "sendbuf");
  recvbuf = gk_imalloc(n*npes, "recvbuf");

  for (i=0; i<n*npes; i++)
    sendbuf[i] = i+myrank;

  BDMPI_Alltoall(sendbuf, n, BDMPI_INT, recvbuf, n, BDMPI_INT, BDMPI_COMM_WORLD);

  for (k=0; k<npes; k++) {
    for (i=0; i<n; i++) {
      if (recvbuf[k*n+i] != myrank*n + i+k)
        printf("[%3d]Error: recvbuf[%d]: got %d instead of %d\n", myrank, k*n+i, recvbuf[k*n+i], myrank*n+i+k);
    }
  }

  gk_free((void **)&sendbuf, &recvbuf, LTERM);


  BDMPI_Barrier(BDMPI_COMM_WORLD);
  if (myrank == 0)
    printf("Testing alltoallv...\n");

  sendbuf    = gk_imalloc(n*npes, "sendbuf");
  recvbuf    = gk_imalloc(n*npes, "recvbuf");
  sendcounts = (size_t *)gk_malloc(npes*sizeof(size_t), "sendcounts");
  recvcounts = (size_t *)gk_malloc(npes*sizeof(size_t), "recvcounts");
  sdispls    = (size_t *)gk_malloc(npes*sizeof(size_t), "sdispls");
  rdispls    = (size_t *)gk_malloc(npes*sizeof(size_t), "rdispls");

  for (iter=0; iter<500; iter++) {
    for (i=0; i<n*npes; i++)
      sendbuf[i] = i+myrank;

    for (i=0; i<npes; i++) {
      sendcounts[i] = recvcounts[i] = n;
      sdispls[i] = rdispls[i] = i*n;
    }

    BDMPI_Alltoallv(sendbuf, sendcounts, sdispls, BDMPI_INT, 
                   recvbuf, recvcounts, rdispls, BDMPI_INT, 
                   BDMPI_COMM_WORLD);

    for (k=0; k<npes; k++) {
      for (i=0; i<n; i++) {
        if (recvbuf[k*n+i] != myrank*n + i+k)
          printf("[%3d]Error: recvbuf[%d]: got %d instead of %d\n", myrank, k*n+i, recvbuf[k*n+i], myrank*n+i+k);
      }
    }

    if (iter%9==0)
      BDMPI_Barrier(BDMPI_COMM_WORLD);

    if (myrank == 0) printf("Done with iter: %d\n", (int)iter);
  }

  gk_free((void **)&sendbuf, &recvbuf, &sendcounts, &recvcounts, &sdispls, &rdispls, LTERM);


  BDMPI_Barrier(BDMPI_COMM_WORLD);
  if (myrank == 0)
    printf("Done.\n");

  BDMPI_Finalize();

  return EXIT_SUCCESS;
}
