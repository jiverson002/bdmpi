/*
 * macros.h
 *
 * This file contains various macro definitions
 *
 * Started 10/5/14
 * Jeremy
 *
 */

#if 0
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

#endif

#define P(...)\
(\
  fprintf(stderr, "[%6ld/%6d:%s,%4d]: ", syscall(SYS_gettid),\
    (int)getpid(), basename(__FILE__), __LINE__),\
  fprintf(stderr, __VA_ARGS__)\
)

#define BDMPL_SLEEP(JOB, MSG, IPC)                                          \
do {                                                                        \
  /*bdprintf("sleep beg@%s:%d\n", basename(__FILE__), __LINE__);*/\
  for (;;) {                                                                \
    memset(&(MSG), 0, sizeof(bdmsg_t));                                     \
    if (-1 == bdmq_recv((JOB)->goMQ, &(MSG), sizeof(bdmsg_t))) {            \
      if (EINTR == errno) {                                                 \
        errno = 0;                                                          \
      }                                                                     \
      else {                                                                \
        bdprintf("Failed on trying to recv a go message in sleep: %s.\n",   \
          strerror(errno));                                                 \
      }                                                                     \
    }                                                                       \
    if (BDMPI_MSGTYPE_PROCEED == (MSG).msgtype)                             \
      break;                                                                \
    /*else\
      bdprintf("received other message %d %d,%d\n", (MSG).msgtype, errno,\
        EINTR);*/\
    /*slv_route(JOB, &(MSG));                                               */\
  }                                                                         \
  /*bdprintf("sleep end@%s:%d\n", basename(__FILE__), __LINE__);*/\
} while (0)
