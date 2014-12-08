/*
 * defs.h
 *
 * This file contains various macro definitions
 *
 * Started 10/5/14
 * Jeremy
 *
 */

#define BDMPL_WITH_SB_DISCARD   /* enable sb_discard() */
//#define BDMPL_WITH_SB_SAVEALL   /* enable sb_saveall() */
#define BDMPL_WITH_SB_NOTIFY    /* enable memory tracking */
//#define BDMPL_WITH_SB_LAZYREAD  /* enable lazy reading (read on-demand) */

#define BDMPL_SLEEP(JOB, MSG)                                   \
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
