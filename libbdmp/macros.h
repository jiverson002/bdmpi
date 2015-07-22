/*
 * macros.h
 *
 * This file contains various macro definitions
 *
 * Started 10/5/14
 * Jeremy
 *
 */

#define P(...)\
(\
  fprintf(stderr, "[%6ld/%6d:%s,%4d]: ", syscall(SYS_gettid),\
    (int)getpid(), basename(__FILE__), __LINE__),\
  fprintf(stderr, __VA_ARGS__)\
)

#define BDMPL_SLEEP(JOB, MSG, IPC)                                          \
do {                                                                        \
  /*bdprintf("sleep beg@%s:%d\n", basename(__FILE__), __LINE__);*/\
  if (-1 == SBMA_release())                                                 \
    bdprintf("Failed on SBMA release: %s.\n", strerror(errno));             \
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
