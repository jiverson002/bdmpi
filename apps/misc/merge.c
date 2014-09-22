/*!
\file
\brief A simple reduce program
\date Started 4/16/2013
\author George
*/

#include <GKlib.h>
#include <bdmpi.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>


#define N 100

/*************************************************************************/
/*! Let the game begin! */
/*************************************************************************/
int main(int argc, char *argv[])
{
  int i, j, myrank, npes;
  pid_t pid, ppid;
  int sids[N], rids[N];
  double svals[N], rvals[N];
  size_t nrecv;

  setbuf(stdout, NULL);
  setbuf(stderr, NULL);

  pid  = getpid();
  ppid = getppid();

  BDMPI_Init(&argc, &argv);

  BDMPI_Comm_size(BDMPI_COMM_WORLD, &npes);
  BDMPI_Comm_rank(BDMPI_COMM_WORLD, &myrank);

  printf("[%d:%d] npes: %d, myrank: %d\n", (int)ppid, (int)pid, npes, myrank);

  srand(10*myrank);

  for (j=0, i=0; i<N; i+=RandomInRange(15), j++) {
    sids[j]  = i;
    svals[j] = 1;
  }
  BDMPI_Barrier(BDMPI_COMM_WORLD);
  BDMPI_Barrier(BDMPI_COMM_WORLD);
  printf("[%2d]IDs: ", myrank);
  for (i=0; i<j; i++)
    printf("%d ", sids[i]);
  printf("\n");
  BDMPI_Barrier(BDMPI_COMM_WORLD);
  BDMPI_Barrier(BDMPI_COMM_WORLD);

  nrecv = N;
  if (myrank == 0)
    BDMPI_Merge(NULL, NULL, 0, rvals, rids, &nrecv, BDMPI_DOUBLE, BDMPI_SUM, 0, BDMPI_COMM_WORLD);
  else
    BDMPI_Merge(svals, sids, j, rvals, rids, &nrecv, BDMPI_DOUBLE, BDMPI_SUM, 0, BDMPI_COMM_WORLD);

  BDMPI_Barrier(BDMPI_COMM_WORLD);
  BDMPI_Barrier(BDMPI_COMM_WORLD);
  if (myrank == 0) {
    printf("Received: %zu elements.\n", nrecv);
    for (i=0; i<nrecv; i++) 
      printf("%4d %le\n", rids[i], rvals[i]);
  }

  BDMPI_Finalize();
}
