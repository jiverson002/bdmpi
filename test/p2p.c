/*!
\file
\brief A program to test point-to-point operations
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
  int myrank, npes, rc, iter;
  pid_t pid, ppid;
  char sendmsg[32768], recvmsg[32768];
  MPI_Status rcv_status;
  struct rlimit rl1, rl2;

  setbuf(stdout, NULL);
  setbuf(stderr, NULL);

  getrlimit(RLIMIT_MSGQUEUE, &rl1);
  getrlimit(RLIMIT_NOFILE, &rl2);
  printf("MSGQUEUE: %d, NOFILE: %d\n", (int)rl1.rlim_max, (int)rl2.rlim_max);

  pid  = getpid();
  ppid = getppid();

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  printf("[%d:%d] npes: %d, myrank: %d\n", (int)ppid, (int)pid, npes, myrank);
  fflush(stdout);

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0) {
    printf("[1]Synchronizing...\n");
    fflush(stdout);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  for (iter=0; iter<1; iter++) {
    if (myrank%2 == 0) {
      printf("[%04d]Sending to %d...\n", myrank, myrank+1);
      fflush(stdout);
  
      sprintf(sendmsg, "Rank %d is sending a message to rank %d.", myrank, (myrank+1));
      rc = MPI_Send(sendmsg, strlen(sendmsg)+1, MPI_CHAR, myrank+1, 1, MPI_COMM_WORLD);
      if (rc != MPI_SUCCESS) 
        printf("[%4d]Send did not succeed. It returned a code of %d\n", myrank, rc);
    }
    else {
      printf("[%04d]Receiving from %d...\n", myrank, myrank-1);
      fflush(stdout);
  
      rc = MPI_Recv(recvmsg, 1024, MPI_CHAR, myrank-1, 1, MPI_COMM_WORLD, &rcv_status);
      if (rc != MPI_SUCCESS) 
        printf("[%4d]Recv did not succeed. It returned a code of %d\n", myrank, rc);
      else 
        printf("[%4d]I got the following message: => %s\n", myrank, recvmsg);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0) {
    printf("[2]Synchronizing...\n");
    fflush(stdout);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if (myrank%2 == 0) {
    printf("[%04d]Sending to %d...\n", myrank, myrank+1);

    sprintf(sendmsg, "Rank %d is sending a message to rank %d.", myrank, (myrank+1));
    sendmsg[1023]  = 'a';
    sendmsg[5023]  = 'b';
    sendmsg[9027]  = 'c';
    sendmsg[11027] = 'd';
    sendmsg[22999] = 'z';
    rc = MPI_Send(sendmsg, 23000, MPI_CHAR, myrank+1, 1, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS) 
      printf("[%4d]Send did not succeed. It returned a code of %d\n", myrank, rc);
  }
  else {
    printf("[%04d]Receiving from %d...\n", myrank, myrank-1);

    rc = MPI_Recv(recvmsg, 23000, MPI_CHAR, myrank-1, 1, MPI_COMM_WORLD, &rcv_status);
    if (rc != MPI_SUCCESS) 
      printf("[%4d]Recv did not succeed. It returned a code of %d\n", myrank, rc);
    else 
      printf("[%4d]I got the following message: => %s [%c%c%c%c%c]\n", myrank, recvmsg,
          recvmsg[1023], recvmsg[5023], recvmsg[9027], recvmsg[11027], recvmsg[22999]);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    printf("[3]Synchronizing...\n");
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();

  return EXIT_SUCCESS;
}
