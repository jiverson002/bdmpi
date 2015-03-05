#include "bdmplib.h"
#include "xmmalloc.h"


static int is_internal=0;
static sjob_t * _job=NULL;


/****************************************************************************/
/*! Charge an allocation to the system */
/****************************************************************************/
static void
_sb_charge(size_t const syspages)
{
  bdmsg_t msg, gomsg;

  if (NULL == _job)
    return;

  if (BDMPI_SB_LAZYWRITE == (_job->jdesc->sbopts&BDMPI_SB_LAZYWRITE)) {
    memset(&msg, 0, sizeof(bdmsg_t));
    msg.msgtype = BDMPI_MSGTYPE_MEMLOAD;
    msg.source  = _job->rank;
    msg.count   = syspages*sysconf(_SC_PAGESIZE);

    if (0 != msg.count) {
      bdmq_send(_job->reqMQ, &msg, sizeof(bdmsg_t));
      BDMPL_SLEEP(_job, gomsg);
    }
  }
}


/****************************************************************************/
/*! Discharge an allocation to the system */
/****************************************************************************/
static void
_sb_discharge(size_t const syspages)
{
  bdmsg_t msg, gomsg;

  if (NULL == _job)
    return;

  if (BDMPI_SB_LAZYWRITE == (_job->jdesc->sbopts&BDMPI_SB_LAZYWRITE)) {
    memset(&msg, 0, sizeof(bdmsg_t));
    msg.msgtype = BDMPI_MSGTYPE_MEMSAVE;
    msg.source  = _job->rank;
    msg.count   = syspages*sysconf(_SC_PAGESIZE);

    if (0 == is_internal && 0 != msg.count) {
      bdmq_send(_job->reqMQ, &msg, sizeof(bdmsg_t));
      BDMPL_SLEEP(_job, gomsg);
    }
  }
}


/*************************************************************************/
/*! API: sb_init */
/*************************************************************************/
extern int
sb_init(sjob_t * const const job)
{
  _job = job;

  (void)xmfstem(job->jdesc->wdir);

  if (-1 == xmmallopt(XMOPT_DEBUG, XMDBG_INFO))
    fprintf(stderr, "sb_init: could not set XMOPT_DEBUG\n");

  if (-1 == xmmallopt(XMOPT_MINPAGES, job->jdesc->sbsize))
    fprintf(stderr, "sb_init: could not set XMOPT_PAGESIZE\n");

  if (-1 == xmmallopt(XMOPT_NUMPAGES, job->jdesc->pgsize))
    fprintf(stderr, "sb_init: could not set XMOPT_PAGESIZE\n");

  if (BDMPI_SB_LAZYWRITE == (job->jdesc->sbopts&BDMPI_SB_LAZYWRITE)) {
    if (-1 == xmbk(&(_sb_charge), &(_sb_discharge)))
      fprintf(stderr, "sb_init: could not set book-keeping functions\n");
  }

  if (BDMPI_SB_LAZYREAD == (job->jdesc->sbopts&BDMPI_SB_LAZYREAD)) {
    if (-1 == xmmallopt(XMOPT_LAZYREAD, 1))
      fprintf(stderr, "sb_init: could not set XMOPT_LAZYREAD\n");
  }

  if (-1 == xmmallopt(XMOPT_MULTITHREAD, job->jdesc->sbnt))
    fprintf(stderr, "sb_init: could not set XMOPT_MULTITHREAD\n");

  if (BDMPI_SB_DLMALLOC == (job->jdesc->sbopts&BDMPI_SB_DLMALLOC)) {
    if (-1 == xmmallopt(XMOPT_DLMALLOC, 1))
      fprintf(stderr, "sb_init: could not set XMOPT_DLMALLOC\n");
  }

  return 1;
}


/*************************************************************************/
/*! API: sb_finalize */
/*************************************************************************/
extern int
sb_finalize(void)
{
  _job = NULL;
  return 1;
}


/*************************************************************************/
/*! API: sb_malloc */
/*************************************************************************/
extern void *
sb_malloc(size_t const len)
{
  return xmmalloc(len);
}


/*************************************************************************/
/*! API: sb_calloc */
/*************************************************************************/
extern void *
sb_calloc(size_t const num, size_t const size)
{
  return xmcalloc(num, size);
}


/*************************************************************************/
/*! API: sb_realloc */
/*************************************************************************/
extern void *
sb_realloc(void * const oldptr, size_t const len)
{
  return xmrealloc(oldptr, len);
}


/****************************************************************************/
/*! API: sb_free */
/****************************************************************************/
extern void
sb_free(void * const addr)
{
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
  (void)xmsync(buf, SIZE_MAX);
}


/*************************************************************************/
/*! API: sb_saveall */
/*************************************************************************/
extern void
sb_saveall(void)
{
  (void)xmsyncall();
}


/*************************************************************************/
/*! API: sb_saveall_internal */
/*************************************************************************/
extern size_t
sb_saveall_internal(void)
{
  size_t num;
  is_internal = 1; /* disable the _sb_discharge function */
  num = xmsyncall();
  is_internal = 0; /* re-enable the _sb_discharge function */
  return num*sysconf(_SC_PAGESIZE);
}


/*************************************************************************/
/*! API: sb_load */
/*************************************************************************/
extern void
sb_load(void const * const buf)
{
  (void)xmload(buf, SIZE_MAX, XMPAGE_SYNC);
}


/*************************************************************************/
/*! API: sb_loadall */
/*************************************************************************/
extern void
sb_loadall(void)
{
  (void)xmloadall(XMPAGE_SYNC);
}


/*************************************************************************/
/* API: sb_discard */
/*************************************************************************/
extern void
sb_discard(void * const ptr, size_t const len)
{
  (void)xmdump(ptr, len);
}
