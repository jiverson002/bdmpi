/*
 * macros.h
 *
 * This file contains various macro definitions
 *
 * Started 10/5/14
 * Jeremy
 *
 */

/*
enum sb_opts
{
  BDMPL_WITH_SB_NONE      = 0,
  BDMPL_WITH_SB_DISCARD   = 1,
  BDMPL_WITH_SB_LAZYWRITE = 2,
  BDMPL_WITH_SB_LAZYREAD  = 4
};
*/

/****************************************************************************/
/*!
 *  \details  Enable the use of sb_discard() throughout the BDMPI library.
 */
/****************************************************************************/
#define BDMPL_WITH_SB_DISCARD

/****************************************************************************/
/*!
 *  \details  Enable the use of sb_saveall() throughout the BDMPI library.
 */
/****************************************************************************/
//#define BDMPL_WITH_SB_SAVEALL

/****************************************************************************/
/*!
 *  \details  Enable the "lazy-write" strategy in the sbmalloc library.  This
 *            means that memory allocations controlled by the sbmalloc library
 *            will not be written to disk until there is ``sufficient''
 *            pressure on the total DRAM to warrant such an action.  In this
 *            case, ``sufficient'' is determined by the resident memory
 *            command line parameter `-rm='.
 *
 *  \note     While compatible, it is not recommended to use this option with
 *            the #BDMPL_WITH_SB_SAVEALL option, since the latter will
 *            essentially negate the advantages of this strategy.
 */
/****************************************************************************/
#define BDMPL_WITH_SB_LAZYWRITE

/****************************************************************************/
/*!
 *  \details  Enable the "lazy-read" strategy in the sbmalloc library.  This
 *            means that memory allocations controlled by the sbmalloc library
 *            will not be read from disk and read protected until the
 *            application makes a read / write attempt to the memory location
 *            corresponding to the allocation.  Furthermore, rather than read
 *            the entire allocation chunk, the first time that any system page
 *            within it is accessed, memory is read and protected at a
 *            resolution of an sbpage, which can be any multiple of a system
 *            page.
 */
/****************************************************************************/
#define BDMPL_WITH_SB_LAZYREAD

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
