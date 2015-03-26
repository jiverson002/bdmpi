#include "bdmplib.h"
#include "sbmalloc.h"
#include "klmalloc.h"


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

  (void)SB_fstem(job->jdesc->wdir);

  if (-1 == SB_mallopt(SBOPT_DEBUG, SBDBG_INFO))
    fprintf(stderr, "sb_init: could not set SBOPT_DEBUG\n");

  if (-1 == SB_mallopt(SBOPT_MINPAGES, job->jdesc->sbsize))
    fprintf(stderr, "sb_init: could not set SBOPT_PAGESIZE\n");

  if (-1 == SB_mallopt(SBOPT_NUMPAGES, job->jdesc->pgsize))
    fprintf(stderr, "sb_init: could not set SBOPT_PAGESIZE\n");

  if (BDMPI_SB_LAZYWRITE == (job->jdesc->sbopts&BDMPI_SB_LAZYWRITE)) {
    if (-1 == SB_acct(&(_sb_charge), &(_sb_discharge)))
      fprintf(stderr, "sb_init: could not set book-keeping functions\n");
  }

  if (BDMPI_SB_LAZYREAD == (job->jdesc->sbopts&BDMPI_SB_LAZYREAD)) {
    if (-1 == SB_mallopt(SBOPT_LAZYREAD, 1))
      fprintf(stderr, "sb_init: could not set SBOPT_LAZYREAD\n");
  }

  if (-1 == SB_mallopt(SBOPT_MULTITHREAD, job->jdesc->sbnt))
    fprintf(stderr, "sb_init: could not set SBOPT_MULTITHREAD\n");

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
  return KL_malloc(len);
}


/*************************************************************************/
/*! API: sb_calloc */
/*************************************************************************/
extern void *
sb_calloc(size_t const num, size_t const size)
{
  return KL_calloc(num, size);
}


/*************************************************************************/
/*! API: sb_realloc */
/*************************************************************************/
extern void *
sb_realloc(void * const oldptr, size_t const len)
{
  return KL_realloc(oldptr, len);
}


/****************************************************************************/
/*! API: sb_free */
/****************************************************************************/
extern void
sb_free(void * const addr)
{
  KL_free(addr);
}


/*************************************************************************/
/*! API: sb_exists */
/*************************************************************************/
extern int
sb_exists(void * const ptr)
{
  return SB_exists(ptr);
}


/*************************************************************************/
/*! API: sb_save */
/*************************************************************************/
extern void
sb_save(void * const buf)
{
  (void)SB_sync(buf, SIZE_MAX);
}


/*************************************************************************/
/*! API: sb_saveall */
/*************************************************************************/
extern void
sb_saveall(void)
{
  (void)SB_syncall();
}


/*************************************************************************/
/*! API: sb_saveall_internal */
/*************************************************************************/
extern size_t
sb_saveall_internal(void)
{
  size_t num;
  is_internal = 1; /* disable the _sb_discharge function */
  num = SB_syncall();
  is_internal = 0; /* re-enable the _sb_discharge function */
  return num*sysconf(_SC_PAGESIZE);
}


/*************************************************************************/
/*! API: sb_load */
/*************************************************************************/
extern void
sb_load(void const * const buf)
{
  (void)SB_load(buf, SIZE_MAX, SBPAGE_SYNC);
}


/*************************************************************************/
/*! API: sb_loadall */
/*************************************************************************/
extern void
sb_loadall(void)
{
  (void)SB_loadall(SBPAGE_SYNC);
}


/*************************************************************************/
/* API: sb_discard */
/*************************************************************************/
extern void
sb_discard(void * const ptr, size_t const len)
{
  (void)SB_dump(ptr, len);
}
