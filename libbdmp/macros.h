/*
 * macros.h
 *
 * This file contains various macro definitions
 *
 * Started 10/5/14
 * Jeremy
 *
 */

/* TODO: I believe that the bug which this hack addresses may actually be
 * related to sb load functionality, because if this hack is place just before
 * BDMPL_SLEEP, but after some other statements which access JOB, then it
 * still fails.  This leads me to believe that the sb load may be broken
 * somehow. */
#define BDMPL_SAVEALL_HACK(JOB)                                         \
do {                                                                    \
  SB_load((JOB), sizeof(sjob_t), SBPAGE_SYNC);                          \
  SB_load((JOB)->goMQ, sizeof(bdmq_t), SBPAGE_SYNC);                    \
  SB_load((JOB)->goMQ->buf, (JOB)->goMQ->msgsize, SBPAGE_DIRTY);        \
} while (0)

#define BDMPL_SLEEP(JOB, MSG)                                           \
do {                                                                    \
  /*bdprintf("sleep beg@%s:%d\n", __FILE__, __LINE__);*/\
  for (;;) {                                                            \
    if (-1 == bdmq_recv((JOB)->goMQ, &(MSG), sizeof(bdmsg_t)))          \
      bdprintf("Failed on trying to recv a go message in sleep: %s.\n", \
        strerror(errno));                                               \
    if (BDMPI_MSGTYPE_PROCEED == (MSG).msgtype)                         \
      break;                                                            \
    slv_route(JOB, &(MSG));                                             \
  }                                                                     \
  /*bdprintf("sleep end@%s:%d\n", __FILE__, __LINE__);*/\
} while (0)
