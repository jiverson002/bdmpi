#include "bdmplib.h"
#include "xmmalloc.h"


static sjob_t * _job=NULL;


/****************************************************************************/
/*! Charge an allocation to the system */
/****************************************************************************/
static void
_sb_charge(size_t const len)
{
#if 0
  if (0 == len) {}
#else
  bdmsg_t msg, gomsg;

  if (BDMPI_SB_LAZYWRITE == (_job->jdesc->sbopts&BDMPI_SB_LAZYWRITE)) {
    memset(&msg, 0, sizeof(bdmsg_t));
    msg.msgtype = BDMPI_MSGTYPE_MEMLOAD;
    msg.source  = _job->rank;
    msg.count   = len*sysconf(_SC_PAGESIZE);
    bdmq_send(_job->reqMQ, &msg, sizeof(bdmsg_t));
    BDMPL_SLEEP(_job, gomsg);
  }
#endif
}


/****************************************************************************/
/*! Discharge an allocation to the system */
/****************************************************************************/
static void
_sb_discharge(size_t const len)
{
#if 0
  if (0 == len) {}
#else
  bdmsg_t msg, gomsg;

  if (BDMPI_SB_LAZYWRITE == (_job->jdesc->sbopts&BDMPI_SB_LAZYWRITE)) {
    memset(&msg, 0, sizeof(bdmsg_t));
    msg.msgtype = BDMPI_MSGTYPE_MEMSAVE;
    msg.source  = _job->rank;
    msg.count   = len*sysconf(_SC_PAGESIZE);
    bdmq_send(_job->reqMQ, &msg, sizeof(bdmsg_t));
    BDMPL_SLEEP(_job, gomsg);
  }
#endif
}


/*************************************************************************/
/*! API: sb_init */
/*************************************************************************/
extern int
sb_init(sjob_t * const const job)
{
  _job = job;

  xmfstem(job->jdesc->wdir);

  if (-1 == xmmallopt(XMOPT_DEBUG, XMDBG_INFO))
    fprintf(stderr, "sb_init: could not set XMOPT_DEBUG\n");

  if (-1 == xmmallopt(XMOPT_MINPAGES, job->jdesc->sbsize))
    fprintf(stderr, "sb_init: could not set XMOPT_PAGESIZE\n");

  if (-1 == xmmallopt(XMOPT_NUMPAGES, job->jdesc->pgsize))
    fprintf(stderr, "sb_init: could not set XMOPT_PAGESIZE\n");

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
  void * newptr=xmmalloc(len);
  _sb_charge(xmprobe(newptr, XMPAGE_SYNC));
  return newptr;
}


/*************************************************************************/
/*! API: sb_calloc */
/*************************************************************************/
extern void *
sb_calloc(size_t const num, size_t const size)
{
  void * newptr=xmcalloc(num, size);
  _sb_charge(xmprobe(newptr, XMPAGE_SYNC));
  return newptr;
}


/*************************************************************************/
/*! API: sb_realloc */
/*************************************************************************/
extern void *
sb_realloc(void * const oldptr, size_t const len)
{
  size_t oldnum, newnum;
  void * newptr;

  oldnum = xmprobe(oldptr, XMPAGE_SYNC|XMPAGE_DIRTY);
  newptr = xmrealloc(oldptr, len);
  newnum = xmprobe(newptr, XMPAGE_SYNC|XMPAGE_DIRTY);

  if (oldnum > newnum)
    _sb_discharge(oldnum-newnum);
  else if (oldnum < newnum)
    _sb_charge(newnum-oldnum);

  return newptr;
}


/****************************************************************************/
/*! API: sb_free */
/****************************************************************************/
extern void
sb_free(void * const addr)
{
  _sb_discharge(xmprobe(addr, XMPAGE_SYNC|XMPAGE_DIRTY));

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
  _sb_discharge(xmsync(buf, SIZE_MAX));
}


/*************************************************************************/
/*! API: sb_saveall */
/*************************************************************************/
extern void
sb_saveall(void)
{
  _sb_discharge(xmsyncall());
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
  _sb_charge(xmload(buf, SIZE_MAX, XMPAGE_SYNC));
}


/*************************************************************************/
/*! API: sb_loadall */
/*************************************************************************/
extern void
sb_loadall(void)
{
  _sb_charge(xmloadall(XMPAGE_SYNC));
}


/*************************************************************************/
/* API: sb_discard */
/*************************************************************************/
extern void
sb_discard(void * const ptr, size_t const len)
{
  _sb_discharge(xmdump(ptr, len));
}
