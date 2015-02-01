#include "bdmplib.h"
#include "xmmalloc.h"


/****************************************************************************/
/*! Charge an allocation to the system */
/****************************************************************************/
static void
_sb_charge(size_t const len)
{
  if (0 == len) {}
#if 0
  bdmsg_t msg, gomsg;

  if (XMOPT_LAZY == xmmallopts[XMOPT_WRITE]) {
    memset(&msg, 0, sizeof(bdmsg_t));
    msg.msgtype = BDMPI_MSGTYPE_MEMLOAD;
    msg.source  = sbinfo->job->rank;
    msg.count   = len;
    bdmq_send(sbinfo->job->reqMQ, &msg, sizeof(bdmsg_t));
    BDMPL_SLEEP(sbinfo->job, gomsg);
  }
#endif
}


/****************************************************************************/
/*! Discharge an allocation to the system */
/****************************************************************************/
static void
_sb_discharge(size_t const len)
{
  if (0 == len) {}
#if 0
  bdmsg_t msg, gomsg;

  if (XMOPT_LAZY == xmmallopts[XMOPT_WRITE]) {
    memset(&msg, 0, sizeof(bdmsg_t));
    msg.msgtype = BDMPI_MSGTYPE_MEMSAVE;
    msg.source  = sbinfo->job->rank;
    msg.count   = len;
    bdmq_send(sbinfo->job->reqMQ, &msg, sizeof(bdmsg_t));
    BDMPL_SLEEP(sbinfo->job, gomsg);
  }
#endif
}


/*************************************************************************/
/*! API: sb_init */
/*************************************************************************/
extern int
sb_init(char const * const fstem, sjob_t * const const job)
{
  xmfstem(fstem);
  if (-1 == xmmallopt(XMOPT_DEBUG, XMDBG_INFO))
    fprintf(stderr, "sb_init: could not set XMOPT_DEBUG\n");
  if (BDMPI_SB_LAZYREAD == (job->jdesc->sbopts&BDMPI_SB_LAZYREAD)) {
    if (-1 == xmmallopt(XMOPT_READ, XMREAD_LAZY))
      fprintf(stderr, "sb_init: could not set XMOPT_READ\n");
  }

  return 1;
}


/*************************************************************************/
/*! API: sb_finalize */
/*************************************************************************/
extern int
sb_finalize(void)
{
  return 1;
}


/*************************************************************************/
/*! API: sb_malloc */
/*************************************************************************/
extern void *
sb_malloc(size_t const len)
{
  // charge
  return xmmalloc(len);
}


/*************************************************************************/
/*! API: sb_calloc */
/*************************************************************************/
extern void *
sb_calloc(size_t const num, size_t const size)
{
  // charge
  return xmcalloc(num, size);
}


/*************************************************************************/
/*! API: sb_realloc */
/*************************************************************************/
extern void *
sb_realloc(void * const oldptr, size_t const len)
{
  // discharge
  //   or
  // charge
  return xmrealloc(oldptr, len);
}


/****************************************************************************/
/*! API: sb_free */
/****************************************************************************/
extern void
sb_free(void * const addr)
{
  // discharge
  xmfree(addr);
}


/*************************************************************************/
/*! API: sb_exists */
/*************************************************************************/
extern int
sb_exists(void * const ptr)
{
  return xmexists(ptr);
}


/*************************************************************************/
/*! API: sb_save */
/*************************************************************************/
extern void
sb_save(void * const buf)
{
  // discharge
  (void)xmsync(buf, SIZE_MAX);
}


/*************************************************************************/
/*! API: sb_saveall */
/*************************************************************************/
extern void
sb_saveall(void)
{
  // discharge
  (void)xmsyncall();
}


/*************************************************************************/
/*! API: sb_saveall_internal */
/*************************************************************************/
extern size_t
sb_saveall_internal(void)
{
  return xmsyncall();
}


/*************************************************************************/
/*! API: sb_load */
/*************************************************************************/
extern void
sb_load(void const * const buf)
{
  // charge
  (void)xmload(buf, SIZE_MAX, XMPAGE_SYNC);
}


/*************************************************************************/
/*! API: sb_loadall */
/*************************************************************************/
extern void
sb_loadall(void)
{
  // charge
  (void)xmloadall(XMPAGE_SYNC);
}


/*************************************************************************/
/* API: sb_discard */
/*************************************************************************/
extern void
sb_discard(void * const ptr, size_t const len)
{
  // discharge
  (void)xmdump(ptr, len);
}
