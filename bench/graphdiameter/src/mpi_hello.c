#include "mpi.h"
#include <stdio.h>
#include <string.h>

#define PING_PONG_BALL 0
#define MSG_MAX 15

int main(int argc, char** argv) {
  int  numtasks, rank, rc;

  rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  /* DEBUG #Tasks */
  // printf ("Number of tasks = %d My rank = %d\n",numtasks,rank);

  /* --- do some circular ping pong --- */

  /* ERROR! inMsg not declared */
//  char *inMsg, *outMsg = "ping_pong_ball";
//  int msgLen = strlen(outMsg);

  /* 1st Solution */
   char inMsg[MSG_MAX] = "ping_pong_ball";
   char *outMsg = inMsg;
   int msgLen = MSG_MAX;

  /* DEBUG Messages */
//  printf("inMsg: %s\n", inMsg);
//  printf("outMsg: %s\n", outMsg);
//  printf("msgLen: %d\n", msgLen);

  MPI_Status status;

  if(rank == 0) {
        int dest = rank + 1,
            src = numtasks - 1;
        MPI_Send(outMsg, msgLen, MPI_CHAR, dest, PING_PONG_BALL, MPI_COMM_WORLD);
        printf("Task %d: Sent %s to task %d\n", rank, outMsg, dest);
        MPI_Recv(inMsg, msgLen, MPI_CHAR, src, PING_PONG_BALL, MPI_COMM_WORLD, &status);
        printf("Task %d: Received %s from task %d\n", rank, inMsg, src);
  }
  else if(rank == (numtasks - 1)) {
      int dest = 0,
          src = numtasks - 2;
      MPI_Recv(inMsg, msgLen, MPI_CHAR, src, PING_PONG_BALL, MPI_COMM_WORLD, &status);
      printf("Task %d: Received %s from task %d\n", rank, inMsg, src);
      MPI_Send(outMsg, msgLen, MPI_CHAR, dest, PING_PONG_BALL, MPI_COMM_WORLD);
      printf("Task %d: Sent %s to task %d\n", rank, outMsg, dest);
  }
  else {
      int dest = rank + 1,
          src = rank - 1;
      MPI_Recv(inMsg, msgLen, MPI_CHAR, src, PING_PONG_BALL, MPI_COMM_WORLD, &status);
      printf("Task %d: Received %s from task %d\n", rank, inMsg, src);
      MPI_Send(outMsg, msgLen, MPI_CHAR, dest, PING_PONG_BALL, MPI_COMM_WORLD);
      printf("Task %d: Sent %s to task %d\n", rank, outMsg, dest);
  }

  int count;
  MPI_Get_count(&status, MPI_CHAR, &count);
  printf("Task %d: Received %d char(s) from task %d with tag %d \n",
         rank, count, status.MPI_SOURCE, status.MPI_TAG);

  MPI_Finalize();
  return 0;
}
