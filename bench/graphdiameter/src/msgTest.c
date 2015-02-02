#include <stdio.h>
#include <string.h>
 
int main(int argc, char *argv[]) 
{
   char *inMsg, *outMsg = "ping_pong_ball"; /* SegFault: inMsg not initialized */
//   int msgLen = strlen(outMsg);

/* CORRECT! */
//   char *inMsg = "ping_pong_ball";
//   char *outMsg = inMsg;

   int inMsgLen = strlen(inMsg);
   int outMsgLen = strlen(outMsg);

   printf("inMsg: %s\n", inMsg);
   printf("outMsg: %s\n", outMsg);
   printf("inMsg has %d chars.\n", inMsgLen);
   printf("outMsg has %d chars.\n", outMsgLen);

   return 0;
}
