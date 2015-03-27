/*
 * macros.h
 *
 * This file contains various macro definitions
 *
 * Started 10/5/14
 * Jeremy
 *
 */

#define BDMPL_SLEEP(JOB, MSG)                                               \
do {                                                                        \
  if (BDMPI_SB_SAVEALL == ((JOB)->jdesc->sbopts&(BDMPI_SB_SAVEALL))) {      \
    if ((JOB)->jdesc->nr < (JOB)->jdesc->ns) {                              \
      /* These are the same conditions under which sb_saveall() is called.
       * goMQ->buf must be loaded with write status (SBPAGE_DIRTY) before
       * calling bdmq_recv, else the system call used by bdmq_recv will
       * generate a SEGFAULT causing the reception of the message to fail.
       */                                                                   \
      SB_load((JOB)->goMQ->buf, (JOB)->goMQ->msgsize, SBPAGE_DIRTY);        \
    }                                                                       \
  }                                                                         \
  /*bdprintf("sleep beg@%s:%d\n", basename(__FILE__), __LINE__);*/\
  for (;;) {                                                                \
    if (-1 == bdmq_recv((JOB)->goMQ, &(MSG), sizeof(bdmsg_t)))              \
      bdprintf("Failed on trying to recv a go message in sleep: %s.\n",     \
        strerror(errno));                                                   \
    if (BDMPI_MSGTYPE_PROCEED == (MSG).msgtype)                             \
      break;                                                                \
    slv_route(JOB, &(MSG));                                                 \
  }                                                                         \
  /*bdprintf("sleep end@%s:%d\n", basename(__FILE__), __LINE__);*/\
} while (0)
