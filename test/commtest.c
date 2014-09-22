/*!
\file
\brief A program to test communicator manipulation routines
\date Started 4/20/2013
\author George
*/

#include <GKlib.h>
#include <bdmpi.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>


#define N 20000

/*************************************************************************/
/*! Let the game begin! */
/*************************************************************************/
int main(int argc, char *argv[])
{
  int i, j, myrank, npes, mynrank, nnpes, mycwrank, cwnpes, mysrank, snpes;
  pid_t pid, ppid;
  MPI_Comm gcomm, rowcomm, colcomm;
  int gmyrank, gnpes, gmynrank, gnnpes, rmyrank, rnpes, cmyrank, cnpes;

  setbuf(stdout, NULL);
  setbuf(stderr, NULL);

  pid  = getpid();
  ppid = getppid();

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  MPI_Comm_size(BDMPI_COMM_NODE, &nnpes);
  MPI_Comm_rank(BDMPI_COMM_NODE, &mynrank);

  MPI_Comm_size(BDMPI_COMM_CWORLD, &cwnpes);
  MPI_Comm_rank(BDMPI_COMM_CWORLD, &mycwrank);

  MPI_Comm_size(MPI_COMM_SELF, &snpes);
  MPI_Comm_rank(MPI_COMM_SELF, &mysrank);


  printf("DEFCOMMS: [%d:%d] sizes:[w,s,cw,n]=[%2d,%2d,%2d,%2d]; ranks:[w,s,cw,n]=[%2d,%2d,%2d,%2d]\n",
      (int)ppid, (int)pid, npes, snpes, cwnpes, nnpes,
      myrank, mysrank, mycwrank, mynrank);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    printf("Duplicating the Node communicator...\n");
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Comm_dup(BDMPI_COMM_NODE, &gcomm);
  MPI_Comm_size(gcomm, &gnnpes);
  MPI_Comm_rank(gcomm, &gmynrank);
  printf("COMM_NODE[%d:%d] npes/nnpes: %d/%d, myrank/mynrank: %d/%d [comm: %p]\n", 
      (int)ppid, (int)pid, npes, gnnpes, myrank, gmynrank, (void *)BDMPI_COMM_NODE);
  MPI_Comm_free(&gcomm);

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    printf("Duplicating the World communicator...\n");
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Comm_dup(MPI_COMM_WORLD, &gcomm);
  MPI_Comm_size(gcomm, &gnpes);
  MPI_Comm_rank(gcomm, &gmyrank);
  printf("COMM_WORLD[%d:%d] gnpes: %d, gmyrank: %d [comm: %p]\n", 
      (int)ppid, (int)pid, gnpes, gmyrank, (void *)gcomm);

  MPI_Barrier(gcomm);

  if (gnpes > 3) {
    MPI_Barrier(gcomm);
    if (myrank == 0)
      printf("Splitting gcomm into 3 groups in two ways.\n");

    MPI_Comm_split(gcomm, gmyrank%3, gmyrank%2, &rowcomm);
    MPI_Comm_size(rowcomm, &rnpes);
    MPI_Comm_rank(rowcomm, &rmyrank);
    printf("[%d:%d] rnpes: %d, rmyrank: %d [comm: %p]\n", (int)ppid, (int)pid, rnpes, rmyrank, (void *)rowcomm);
    MPI_Barrier(rowcomm);

    MPI_Comm_split(gcomm, gmyrank/((gnpes+2)/3), gmyrank%2, &colcomm);
    MPI_Comm_size(colcomm, &cnpes);
    MPI_Comm_rank(colcomm, &cmyrank);
    printf("[%d:%d] cnpes: %d, cmyrank: %d [comm: %p]\n", (int)ppid, (int)pid, cnpes, cmyrank, (void *)colcomm);
    MPI_Barrier(colcomm);

    printf("[%d:%d] gnpes: %d, rnpes: %d, cnpes: %d, gmyrank: %d, rmyrank: %d, cmyrank: %d\n",
        (int)ppid, (int)pid, gnpes, rnpes, cnpes, gmyrank, rmyrank, cmyrank);

    MPI_Comm_free(&rowcomm);
    MPI_Comm_free(&colcomm);
  }


  if (gnpes%3 == 0) {
    MPI_Barrier(gcomm);
    if (myrank == 0)
      printf("Splitting gcomm into a %dx%d grid.\n", 3, gnpes/3);

    MPI_Comm_split(gcomm, gmyrank%3, 1, &rowcomm);
    MPI_Comm_size(rowcomm, &rnpes);
    MPI_Comm_rank(rowcomm, &rmyrank);
    printf("[%d:%d] rnpes: %d, rmyrank: %d [comm: %p]\n", (int)ppid, (int)pid, rnpes, rmyrank, (void *)rowcomm);
    MPI_Barrier(rowcomm);

    MPI_Comm_split(gcomm, gmyrank/3, 1, &colcomm);
    MPI_Comm_size(colcomm, &cnpes);
    MPI_Comm_rank(colcomm, &cmyrank);
    printf("[%d:%d] cnpes: %d, cmyrank: %d [comm: %p]\n", (int)ppid, (int)pid, cnpes, cmyrank, (void *)colcomm);
    MPI_Barrier(colcomm);

    printf("[%d:%d] [%d %d] gnpes: %d, rnpes: %d, cnpes: %d, gmyrank: %d, rmyrank: %d, cmyrank: %d\n",
        (int)ppid, (int)pid, gmyrank%3, gmyrank/3, gnpes, rnpes, cnpes, gmyrank, rmyrank, cmyrank);

    MPI_Comm_free(&rowcomm);
    MPI_Comm_free(&colcomm);
  }

  MPI_Barrier(gcomm);
  if (myrank == 0)
    printf("Freeing the gcomm.\n");
  MPI_Comm_free(&gcomm);


  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    printf("Done.\n");

  MPI_Finalize();

  return EXIT_SUCCESS;
}
