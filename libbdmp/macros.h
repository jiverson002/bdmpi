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
      /* these are the same conditions under which sb_saveall() is called.
       */                                                                   \
      SB_load((JOB), sizeof(sjob_t), SBPAGE_SYNC);                          \
      SB_load((JOB)->goMQ->buf, (JOB)->goMQ->msgsize, SBPAGE_DIRTY);        \
    }                                                                       \
  }                                                                         \
  /*bdprintf("sleep beg@%s:%d\n", __FILE__, __LINE__);*/\
  for (;;) {                                                                \
    if (-1 == bdmq_recv((JOB)->goMQ, &(MSG), sizeof(bdmsg_t)))              \
      bdprintf("Failed on trying to recv a go message in sleep: %s.\n",     \
        strerror(errno));                                                   \
    if (BDMPI_MSGTYPE_PROCEED == (MSG).msgtype)                             \
      break;                                                                \
    slv_route(JOB, &(MSG));                                                 \
  }                                                                         \
  /*bdprintf("sleep end@%s:%d\n", __FILE__, __LINE__);*/\
} while (0)
