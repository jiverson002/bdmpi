/*
 * defs.h
 *
 * This file contains various macro definitions
 *
 * Started 10/5/14
 * Jeremy
 *
 */

#define BDMPI_SLEEP(JOB, MSG)                                   \
do {                                                            \
  for (;;) {                                                    \
    if (-1 == bdmq_recv((JOB)->goMQ, &(MSG), sizeof(bdmsg_t)))  \
      bdprintf("Failed on trying to recv a go message: %s.\n",  \
        strerror(errno));                                       \
    if (BDMPI_MSGTYPE_PROCEED == (MSG).msgtype)                 \
      break;                                                    \
    slv_route(JOB, &(MSG));                                     \
  }                                                             \
} while (0)
