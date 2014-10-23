/*!
\file
\brief Implements the storage-backed malloc interface
\date Started 4/1/2014
\author George
*/


#include "bdmplib.h"


/* Stores information associated with a storage-backed memory chunk */
typedef struct sbchunk {
  size_t saddr, eaddr;   /* starting/ending address of the anonymous mapping */
  size_t npages;         /* number of pages allocated */
  size_t nbytes;         /* number of bytes requested */
  uint8_t flags;         /* global chunk-level flag */
  uint8_t *pflags;       /* per-page flag vector */
  char *fname;           /* the file that will store the data */
  pthread_mutex_t mtx;   /* the mutex guarding it */
  struct sbchunk *next;  /* pointer to the next chunk */
} sbchunk_t;

/* Stores global information associated with storage-backed memory */
typedef struct {
  size_t pagesize;       /* the size of a memory page */
  size_t minsize;        /* the minimum allocation in pages handled by sbmalloc */
  size_t npages;         /* the number of pages that defines the unit block */
  char *fstem;           /* the file stem where the data is stored */
  sjob_t *job;           /* the slave job information */
  sbchunk_t *head;       /* the first sbchunk */
  pthread_mutex_t mtx;   /* the mutex guarding it */
  pthread_mutexattr_t mtx_attr;  /* attributes for all mutexes */

  struct sigaction act, oldact;  /* for the SIGSEGV signal handler */
} sbinfo_t;


/* static global variables */
static sbinfo_t *sbinfo=NULL;
__thread size_t last_addr=0;


/* constants */
#define SBCHUNK_NONE    1
#define SBCHUNK_READ    2
#define SBCHUNK_WRITE   4
#define SBCHUNK_ONDISK  8


/* hooks to build-in function */
static void *(*libc_malloc)(size_t) = NULL;
static void *(*libc_calloc)(size_t, size_t) = NULL;
static void *(*libc_realloc)(void *, size_t) = NULL;
static void (*libc_free)(void *) = NULL;
static ssize_t (*libc_read)(int, void *, size_t) = NULL;
static ssize_t (*libc_write)(int, const void *, size_t) = NULL;
static size_t (*libc_fread)(void *, size_t, size_t, FILE *) = NULL;
static size_t (*libc_fwrite)(const void *, size_t, size_t, FILE *) = NULL;
static int (*libc_mlock)(const void *, size_t) = NULL;
static int (*libc_munlock)(const void *, size_t) = NULL;
static int (*libc_mlockall)(int) = NULL;
static int (*libc_munlockall)(void) = NULL;


/* private prototypes */
sbchunk_t *_sb_find(void *ptr);
static void _sb_handler(int sig, siginfo_t *si, void *unused);
void _sb_chunkload(sbchunk_t *sbchunk);
void _sb_chunksave(sbchunk_t *sbchunk);
size_t _sb_chunksave_internal(sbchunk_t *sbchunk);
void _sb_chunkfree(sbchunk_t *sbchunk);
void _sb_chunkfreeall();

#ifdef XXX
/* macros */
#define BD_GET_LOCK(lock) GKASSERT(pthread_mutex_lock(lock) == 0)
#define BD_LET_LOCK(lock) GKASSERT(pthread_mutex_unlock(lock) == 0)
#endif



/*************************************************************************/
/*! Hook: malloc */
/*************************************************************************/
void *malloc(size_t nbytes)
{
  if (libc_malloc == NULL)
    *((void **) &libc_malloc) = dlsym(RTLD_NEXT, "malloc");

  if (sbinfo == NULL)
    return libc_malloc(nbytes);
  else {
    if (sbinfo->minsize == 0 || nbytes <= sbinfo->minsize)
      return libc_malloc(nbytes);
    else
      return sb_malloc(nbytes);
  }

  return NULL;
}


/*************************************************************************/
/*! Hook: calloc */
/*************************************************************************/
void *calloc_no(size_t nmemb, size_t size)
{
  printf("Hi there2\n");
  if (libc_calloc == NULL)
    *((void **) &libc_calloc) = dlsym(RTLD_NEXT, "calloc");

  return libc_calloc(nmemb, size);

  if (sbinfo == NULL)
    return libc_calloc(nmemb, size);
  else {
    if (sbinfo->minsize == 0 || nmemb*size <= sbinfo->minsize)
      return libc_calloc(nmemb, size);
    else
      return sb_malloc(nmemb*size);
  }

  return NULL;
}


/*************************************************************************/
/*! Hook: realloc */
/*************************************************************************/
void *realloc(void *ptr, size_t size)
{
  if (libc_realloc == NULL)
    *((void **) &libc_realloc) = dlsym(RTLD_NEXT, "realloc");

  if (sbinfo == NULL)
    return libc_realloc(ptr, size);
  else {
    if (sb_exists(ptr))
      return sb_realloc(ptr, size);
    else
      return libc_realloc(ptr, size);
  }

  return NULL;

}


/*************************************************************************/
/*! Hook: free */
/*************************************************************************/
void free(void *ptr)
{
  if (libc_free == NULL)
    *((void **) &libc_free) = dlsym(RTLD_NEXT, "free");

  if (sbinfo == NULL) {
    libc_free(ptr);
  }
  else {
    if (sb_exists(ptr))
      sb_free(ptr);
    else
      libc_free(ptr);
  }

  return;
}


/*************************************************************************/
/*! Hook: read */
/*************************************************************************/
ssize_t read(int fd, void *buf, size_t count)
{
  if (libc_read == NULL)
    *((void **) &libc_read) = dlsym(RTLD_NEXT, "read");

  if (sbinfo == NULL) {
    return libc_read(fd, buf, count);
  }
  else {
    if (sb_exists(buf))
      memset(buf, 0, count);

    return libc_read(fd, buf, count);
  }
}


/*************************************************************************/
/*! Hook: write */
/*************************************************************************/
ssize_t write(int fd, const void *buf, size_t count)
{
  if (libc_write == NULL)
    *((void **) &libc_write) = dlsym(RTLD_NEXT, "write");

  if (count > 0)
    sb_load((void *)buf);

  return libc_write(fd, buf, count);
}


/*************************************************************************/
/*! Hook: fread */
/*************************************************************************/
size_t fread(void *buf, size_t size, size_t nmemb, FILE *stream)
{
  if (libc_fread == NULL)
    *((void **) &libc_fread) = dlsym(RTLD_NEXT, "fread");

  if (sbinfo == NULL) {
    return libc_fread(buf, size, nmemb, stream);
  }
  else {
    if (sb_exists(buf))
      memset(buf, 0, size*nmemb);

    return libc_fread(buf, size, nmemb, stream);
  }
}


/*************************************************************************/
/*! Hook: fwrite */
/*************************************************************************/
size_t fwrite(const void *buf, size_t size, size_t nmemb, FILE *stream)
{
  if (libc_fwrite == NULL)
    *((void **) &libc_fwrite) = dlsym(RTLD_NEXT, "fwrite");

  if (nmemb > 0)
    sb_load((void *)buf);

  return libc_fwrite(buf, size, nmemb, stream);
}


/*************************************************************************/
/*! Hook: mlock */
/*************************************************************************/
int mlock(const void *addr, size_t len)
{
  if (libc_mlock == NULL)
    *((void **) &libc_mlock) = dlsym(RTLD_NEXT, "mlock");

  sb_load((void *)addr);

  return libc_mlock(addr, len);
}


/*************************************************************************/
/*! Hook: munlock */
/*************************************************************************/
int munlock(const void *addr, size_t len)
{
  if (libc_munlock == NULL)
    *((void **) &libc_munlock) = dlsym(RTLD_NEXT, "munlock");

  return libc_munlock(addr, len);
}


/*************************************************************************/
/*! Hook: mlockall */
/*************************************************************************/
int mlockall(int flags)
{
  if (libc_mlockall == NULL)
    *((void **) &libc_mlockall) = dlsym(RTLD_NEXT, "mlockall");

  sb_loadall();

  return libc_mlockall(flags);
}


/*************************************************************************/
/*! Hook: munlockall */
/*************************************************************************/
int munlockall(void)
{
  if (libc_munlockall == NULL)
    *((void **) &libc_munlockall) = dlsym(RTLD_NEXT, "munlockall");

  return libc_munlockall();
}


/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

/*************************************************************************/
/*! Initializes the sbmalloc subsystem */
/*************************************************************************/
int sb_init(char *fstem, sjob_t * const job, size_t minsize, size_t npages)
{
  if (sbinfo != NULL) {
    perror("sb_init: sbinfo != NULL");
    exit(EXIT_FAILURE);
  }

  /* hook to the libc version of the key functions */
  *((void **) &libc_malloc)     = dlsym(RTLD_NEXT, "malloc");
  *((void **) &libc_calloc)     = dlsym(RTLD_NEXT, "calloc");
  *((void **) &libc_realloc)    = dlsym(RTLD_NEXT, "realloc");
  *((void **) &libc_free)       = dlsym(RTLD_NEXT, "free");
  *((void **) &libc_read)       = dlsym(RTLD_NEXT, "read");
  *((void **) &libc_write)      = dlsym(RTLD_NEXT, "write");
  *((void **) &libc_fread)      = dlsym(RTLD_NEXT, "fread");
  *((void **) &libc_fwrite)     = dlsym(RTLD_NEXT, "fwrite");
  *((void **) &libc_mlock)      = dlsym(RTLD_NEXT, "mlock");
  *((void **) &libc_munlock)    = dlsym(RTLD_NEXT, "munlock");
  *((void **) &libc_mlockall)   = dlsym(RTLD_NEXT, "mlockall");
  *((void **) &libc_munlockall) = dlsym(RTLD_NEXT, "munlockall");

  if ((sbinfo = libc_malloc(sizeof(sbinfo_t))) == NULL) {
    perror("Failed to allocate sbinfo\n");
    return 0;
  }
  memset(sbinfo, 0, sizeof(sbinfo_t));

  sbinfo->minsize  = minsize*sysconf(_SC_PAGESIZE);
  sbinfo->npages   = npages;
  sbinfo->pagesize = npages*sysconf(_SC_PAGESIZE);
  sbinfo->head     = NULL;

  if ((sbinfo->fstem = (char *)libc_malloc(strlen(fstem)+1)) == NULL) {
    perror("sbinit: failed to allocate memory for fstem\n");
    goto ERROR_EXIT;
  }
  strcpy(sbinfo->fstem, fstem);

  /* setup the signal handler */
  sbinfo->act.sa_flags = SA_SIGINFO;
  sigemptyset(&(sbinfo->act.sa_mask));
  sbinfo->act.sa_sigaction = _sb_handler;
  if (sigaction(SIGSEGV, &(sbinfo->act), &(sbinfo->oldact)) == -1) {
    perror("sbinit: failed to install the signal handler\n");
    goto ERROR_EXIT;
  }

  /* setup the mutex and global attributes */
  GKASSERT(pthread_mutexattr_init(&(sbinfo->mtx_attr)) == 0);
  GKASSERT(pthread_mutexattr_settype(&(sbinfo->mtx_attr), PTHREAD_MUTEX_RECURSIVE) == 0);
  GKASSERT(pthread_mutex_init(&(sbinfo->mtx), &(sbinfo->mtx_attr)) == 0);

  sbinfo->job = job;

  return 1;

ERROR_EXIT:
  if (sbinfo->fstem)
    libc_free(sbinfo->fstem);

  libc_free(sbinfo);
  sbinfo = NULL;

  return 0;
}


/*************************************************************************/
/*! Shutsdown the sbmalloc subsystem */
/*************************************************************************/
int sb_finalize()
{
  if (sbinfo == NULL) {
    perror("sbfinalize: sbinfo -= NULL");
    exit(EXIT_FAILURE);
  }

  if (sigaction(SIGSEGV, &(sbinfo->oldact), NULL) == -1) {
    perror("Failed to install the old signal handler\n");
    return 0;
  }

  GKASSERT(pthread_mutex_destroy(&(sbinfo->mtx)) == 0);
  GKASSERT(pthread_mutexattr_destroy(&(sbinfo->mtx_attr)) == 0);

  _sb_chunkfreeall();
  libc_free(sbinfo->fstem);
  libc_free(sbinfo);
  sbinfo = NULL;

  return 1;
}


/*************************************************************************/
/*! Allocate memory via anonymous mmap */
/*************************************************************************/
void *sb_malloc(size_t nbytes)
{
  ssize_t ip;
  sbchunk_t *sbchunk;

  //printf("sb_malloc: allocating %zu bytes\n", nbytes);

  if (sbinfo == NULL) {
    perror("sb_malloc: sbinfo == NULL");
    exit(EXIT_FAILURE);
  }

  sbchunk = libc_malloc(sizeof(sbchunk_t));
  if (sbchunk == NULL)
    return NULL;

  sbchunk->nbytes = nbytes;

  /* determine the allocation size in terms of pagesize */
  sbchunk->npages = (nbytes+sbinfo->pagesize-1)/sbinfo->pagesize;

  /* allocate the flag array for the pages */
  sbchunk->pflags = (uint8_t *)libc_calloc(sbchunk->npages+1, sizeof(uint8_t));
  if (sbchunk->pflags == NULL) {
    perror("sbmalloc.1 failure\n");
    libc_free(sbchunk);
    return NULL;
  }

  //printf("Allocating %zu bytes\n", sbchunk->npages*sbinfo->pagesize);
  sbchunk->saddr = (size_t) mmap(NULL, sbchunk->npages*sbinfo->pagesize, PROT_NONE,
                                MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);

  if ((void *)sbchunk->saddr == MAP_FAILED) {
    perror("sbmalloc.2 failure");
    libc_free(sbchunk->pflags);
    libc_free(sbchunk);
    return NULL;
  }
  sbchunk->eaddr = sbchunk->saddr+nbytes;


  /* create the filename for storage purposes */
  if ((sbchunk->fname = (char *)libc_malloc(100+strlen(sbinfo->fstem))) == NULL) {
    munmap((void *)sbchunk->saddr, sbchunk->npages*sbinfo->pagesize);
    libc_free(sbchunk->pflags);
    libc_free(sbchunk);
    return NULL;
  }

  if (sprintf(sbchunk->fname, "%s%zx", sbinfo->fstem, sbchunk->saddr) == -1) {
    munmap((void *)sbchunk->saddr, sbchunk->npages*sbinfo->pagesize);
    libc_free(sbchunk->pflags);
    libc_free(sbchunk);
    return NULL;
  }

#if 1
  sbchunk->flags = SBCHUNK_NONE;
#else
  /*----------------------------------------------------------------------*/
#ifdef BDMPL_WITH_SB_NOTFIY
{
  bdmsg_t msg, gomsg;

  /* notify the master that you want to load memory */
  msg.msgtype = BDMPI_MSGTYPE_MEMLOAD;
  msg.source  = sbinfo->job->rank;
  msg.count   = sbchunk->npages*sbinfo->pagesize;
  bdmq_send(sbinfo->job->reqMQ, &msg, sizeof(bdmsg_t));
  BDMPL_SLEEP(sbinfo->job, gomsg);
}
#endif
  /*----------------------------------------------------------------------*/


  if (-1 == mprotect((void *)sbchunk->saddr, sbchunk->npages*sbinfo->pagesize, PROT_READ)) {
    perror("sb_malloc: failed to PROT_READ");
    exit(EXIT_FAILURE);
  }

  sbchunk->flags = SBCHUNK_READ;
  for (ip=0; ip<sbchunk->npages; ip++)
    sbchunk->pflags[ip] = SBCHUNK_READ;
#endif

  /* initialize the mutex */
  GKASSERT(pthread_mutex_init(&(sbchunk->mtx), &(sbinfo->mtx_attr)) == 0);

  BD_GET_LOCK(&(sbinfo->mtx));
  sbchunk->next = sbinfo->head;
  sbinfo->head  = sbchunk;
  BD_LET_LOCK(&(sbinfo->mtx));

  return (void *)sbchunk->saddr;
}


/*************************************************************************/
/*! Reallocates memory via anonymous mmap */
/*************************************************************************/
void *sb_realloc(void *oldptr, size_t nbytes)
{
  size_t ip, new_saddr, new_npages;
  uint8_t *new_pflags=NULL;
  sbchunk_t *sbchunk;

  if (sbinfo == NULL) {
    perror("sb_realloc: sbinfo == NULL");
    exit(EXIT_FAILURE);
  }

  new_npages = (nbytes+sbinfo->pagesize-1)/sbinfo->pagesize;

  BD_GET_LOCK(&(sbinfo->mtx));
  if ((sbchunk = _sb_find(oldptr)) == NULL) {
    perror("sb_realloc: failed to find the sbchunk");
    exit(EXIT_FAILURE);
  }
  BD_LET_LOCK(&(sbinfo->mtx));

  BD_GET_LOCK(&(sbchunk->mtx));

  /* see if we are shrinking */
  if (nbytes <= sbchunk->nbytes) { /* easy case */
    /*----------------------------------------------------------------------*/
#ifdef BDMPL_WITH_SB_NOTIFY
    if (!(sbchunk->flags&SBCHUNK_NONE)) {
      bdmsg_t msg, gomsg;

      /* notify the master that you want to load memory */
      msg.msgtype = BDMPI_MSGTYPE_MEMSAVE;
      msg.source  = sbinfo->job->rank;
      msg.count   = (sbchunk->npages-new_npages)*sbinfo->pagesize;
      bdmq_send(sbinfo->job->reqMQ, &msg, sizeof(bdmsg_t));
      BDMPL_SLEEP(sbinfo->job, gomsg);
    }
#endif
    /*----------------------------------------------------------------------*/


    /* set the now unused pages as DONTNEED */
    if (madvise((void *)(sbchunk->saddr+new_npages*sbinfo->pagesize),
          (sbchunk->npages-new_npages)*sbinfo->pagesize, MADV_DONTNEED) == -1) {
      perror("sb_realloc: failed to MADV_DONTNEED");
      goto ERROR_EXIT;
    }

    sbchunk->nbytes = nbytes;
    sbchunk->npages = new_npages;
    sbchunk->eaddr  = sbchunk->saddr+nbytes;
    sbchunk->pflags[new_npages] = 0;  /* this is for the +1 reset */
  }
  else {
    /*----------------------------------------------------------------------*/
#ifdef BDMPL_WITH_SB_NOTIFY
{
    bdmsg_t msg, gomsg;

    /* notify the master that you want to load memory */
    msg.msgtype = BDMPI_MSGTYPE_MEMLOAD;
    msg.source  = sbinfo->job->rank;
    msg.count   = (new_npages-sbchunk->npages)*sbinfo->pagesize;
    bdmq_send(sbinfo->job->reqMQ, &msg, sizeof(bdmsg_t));
    BDMPL_SLEEP(sbinfo->job, gomsg);
}
#endif
    /*----------------------------------------------------------------------*/


    /* allocate the new pflags */
    if ((new_pflags = (uint8_t *)libc_calloc(new_npages+1, sizeof(uint8_t))) == NULL) {
      perror("sb_realloc: failed to malloc new_pflags");
      goto ERROR_EXIT;
    }

    /* get the new anonymous map */
    new_saddr = (size_t) mmap(NULL, new_npages*sbinfo->pagesize, PROT_READ|PROT_WRITE,
                             MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);

    if ((void *)new_saddr == MAP_FAILED) {
      perror("sb_realloc: failed to get a new mmap");
      goto ERROR_EXIT;
    }

    /* load the sbchunk if not accessible */
    if (sbchunk->flags&SBCHUNK_NONE)
      _sb_chunkload(sbchunk);

    /* copy the data */
    memcpy((void *)new_saddr, (void *)sbchunk->saddr, sbchunk->nbytes);

    /* unmap the old */
    if (munmap((void *)sbchunk->saddr, sbchunk->npages*sbinfo->pagesize) == -1) {
      perror("sb_realloc: failed to unmap the old region");
      exit(EXIT_FAILURE);
    }

    sbchunk->nbytes = nbytes;
    sbchunk->npages = new_npages;
    sbchunk->saddr  = new_saddr;
    sbchunk->eaddr  = sbchunk->saddr+nbytes;

    /* TODO: does the whole chunk need to be re-written */
    /* set remapped sbchunk for writing and remove its ondisk flag */
    sbchunk->flags |= SBCHUNK_WRITE;
    if (sbchunk->flags&SBCHUNK_ONDISK)
      sbchunk->flags ^= SBCHUNK_ONDISK;

    /* set per-page flags */
    for (ip=0; ip<sbchunk->npages; ip++) {
      new_pflags[ip] = sbchunk->pflags[ip];
      if (sbchunk->pflags[ip]&SBCHUNK_ONDISK) {
        new_pflags[ip] ^= SBCHUNK_ONDISK;
        new_pflags[ip] |= SBCHUNK_WRITE;
      }
    }

    libc_free(sbchunk->pflags);
    sbchunk->pflags = new_pflags;
  }

  BD_LET_LOCK(&(sbchunk->mtx));

  return (void *)sbchunk->saddr;

ERROR_EXIT:
  BD_LET_LOCK(&(sbchunk->mtx));
  if (new_pflags != NULL)
    libc_free(new_pflags);

  return NULL;
}


/*************************************************************************/
/*! Frees the allocated sbmemory */
/*************************************************************************/
void sb_free(void *buf)
{
  sbchunk_t *ptr, *pptr=NULL;
  size_t addr=(size_t)buf;

  if (sbinfo == NULL) {
    perror("sb_free: sbinfo == NULL");
    exit(EXIT_FAILURE);
  }

  BD_GET_LOCK(&(sbinfo->mtx));
  /* find the sbchunk */
  for (ptr=sbinfo->head; ptr!=NULL; ptr=ptr->next) {
    if (ptr->saddr == addr)
      break;
    pptr = ptr;
  }
  if (ptr == NULL) {
    printf("Got a free for an unhandled memory location: %zx\n", addr);
    exit(EXIT_FAILURE);
  }

  /* update the link-list */
  if (pptr == NULL)
    sbinfo->head = ptr->next;
  else
    pptr->next = ptr->next;

  BD_LET_LOCK(&(sbinfo->mtx));

  /* destroy the mutex */
  GKASSERT(pthread_mutex_destroy(&(ptr->mtx)) == 0);

  /* free the chunk */
  _sb_chunkfree(ptr);

  return;
}


/*************************************************************************/
/*! Returns true if the supplied pointer is been handled by sbmalloc */
/*************************************************************************/
int sb_exists(void *ptr)
{
  int status;

  BD_GET_LOCK(&(sbinfo->mtx));
  status = (_sb_find(ptr) == NULL ? 0 : 1);
  BD_LET_LOCK(&(sbinfo->mtx));

  return status;
}


/*************************************************************************/
/*! Saves to disk the sbchunk associated with the supplied pointer */
/*************************************************************************/
void sb_save(void *buf)
{
  sbchunk_t *sbchunk;

  if (sbinfo == NULL)
    return;

  BD_GET_LOCK(&(sbinfo->mtx));
  sbchunk = _sb_find(buf);
  BD_LET_LOCK(&(sbinfo->mtx));

  if (sbchunk == NULL)
    return;

  BD_GET_LOCK(&(sbchunk->mtx));
  _sb_chunksave(sbchunk);
  BD_LET_LOCK(&(sbchunk->mtx));
}


/*************************************************************************/
/*! Saves to disk all sbchunks */
/*************************************************************************/
void sb_saveall()
{
  size_t count=0;
  sbchunk_t *sbchunk;

  if (sbinfo == NULL)
    return;

  BD_GET_LOCK(&(sbinfo->mtx));
  for (sbchunk=sbinfo->head; sbchunk!=NULL; sbchunk=sbchunk->next) {
    BD_GET_LOCK(&(sbchunk->mtx));
    if (!(sbchunk->flags&SBCHUNK_NONE))
      count += sbchunk->npages*sbinfo->pagesize;
    _sb_chunksave(sbchunk);
    BD_LET_LOCK(&(sbchunk->mtx));
  }
  BD_LET_LOCK(&(sbinfo->mtx));

  //bdprintf("[%3d] sb_saveall(%zu)\n", sbinfo->job->rank, count);
}

size_t sb_saveall_internal()
{
  size_t count=0;
  sbchunk_t *sbchunk;

  if (sbinfo == NULL)
    return 0;

  BD_GET_LOCK(&(sbinfo->mtx));
  for (sbchunk=sbinfo->head; sbchunk!=NULL; sbchunk=sbchunk->next) {
    BD_GET_LOCK(&(sbchunk->mtx));
    count += _sb_chunksave_internal(sbchunk);
    BD_LET_LOCK(&(sbchunk->mtx));
  }
  BD_LET_LOCK(&(sbinfo->mtx));

  //bdprintf("[%3d] sb_saveall(%zu)\n", sbinfo->job->rank, count);

  return count;
}


/*************************************************************************/
/*! Loads the sbchunk associated with the supplied pointer */
/*************************************************************************/
void sb_load(void *buf)
{
  sbchunk_t *sbchunk;

  if (sbinfo == NULL)
    return;

  BD_GET_LOCK(&(sbinfo->mtx));
  sbchunk = _sb_find(buf);
  BD_LET_LOCK(&(sbinfo->mtx));

  if (sbchunk == NULL)
    return;

  BD_GET_LOCK(&(sbchunk->mtx));
  if (sbchunk->flags&SBCHUNK_NONE)
    _sb_chunkload(sbchunk);
  BD_LET_LOCK(&(sbchunk->mtx));
}


/*************************************************************************/
/*! Loads all unloaded chunks from disk */
/*************************************************************************/
void sb_loadall()
{
  sbchunk_t *sbchunk;

  if (sbinfo == NULL)
    return;

  BD_GET_LOCK(&(sbinfo->mtx));
  for (sbchunk=sbinfo->head; sbchunk!=NULL; sbchunk=sbchunk->next) {
    BD_GET_LOCK(&(sbchunk->mtx));
    if (sbchunk->flags&SBCHUNK_NONE)
      _sb_chunkload(sbchunk);
    BD_LET_LOCK(&(sbchunk->mtx));
  }
  BD_LET_LOCK(&(sbinfo->mtx));
}


/*************************************************************************/
/*! The specified region of memory should be treated as being newly allocated.
    Any disk-backed storage and/or write modifications are discarded.
    For now, this should be called right before an sb_saveall() call as
    it does not properly handles the permission changes. In principle,
    it should be given a READ permission and removed its WRITE permissions.

    If size is -1, then the entire allocation associated with ptr is
    discarded.

*/
/*************************************************************************/
void sb_discard(void *ptr, ssize_t size)
{
  size_t addr, ip, ifirst, iend;
  sbchunk_t *sbchunk;

  if (sbinfo == NULL)
    return;

  /* find the sbchunk */
  BD_GET_LOCK(&(sbinfo->mtx));
  sbchunk = _sb_find(ptr);
  BD_LET_LOCK(&(sbinfo->mtx));

  if (sbchunk == NULL)
    return; /* return if we are not handling it; TODO: Maybe a madvise(MADV_DONTNEED) */

  addr = (size_t)ptr;

  BD_GET_LOCK(&(sbchunk->mtx));
  if (size == -1) { /* discard entire allocation */
    ifirst = 0;
    iend   = sbchunk->npages;

    if (sbchunk->flags&SBCHUNK_WRITE)
      sbchunk->flags ^= SBCHUNK_WRITE;
    if (sbchunk->flags&SBCHUNK_ONDISK)
      sbchunk->flags ^= SBCHUNK_ONDISK;
  }
  else { /* discard supplied range */
    ifirst = (addr - sbchunk->saddr + sbinfo->pagesize-1)/sbinfo->pagesize;
    iend   = (addr+size - sbchunk->saddr)/sbinfo->pagesize;
  }

  /* provide only read permissions, assuming that it had read permissions; this
     will remove the write permissions so in case we do not block, subsequent
     writes will be intercepted correctly */
  if (sbchunk->flags&SBCHUNK_READ) {
    if (mprotect((void *)(sbchunk->saddr+ifirst*sbinfo->pagesize),
                 (iend-ifirst)*sbinfo->pagesize, PROT_READ) == -1) {
      perror("sb_discard: failed to mprotect to PROT_READ");
      exit(EXIT_FAILURE);
    }
  }

  /* update the corresponding pflags[] entries */
  for (ip=ifirst; ip<iend; ip++) {
    if (sbchunk->pflags[ip]&SBCHUNK_WRITE)
      sbchunk->pflags[ip] ^= SBCHUNK_WRITE;
    if (sbchunk->pflags[ip]&SBCHUNK_ONDISK)
      sbchunk->pflags[ip] ^= SBCHUNK_ONDISK;
  }
  BD_LET_LOCK(&(sbchunk->mtx));
}


/*************************************************************************/
/*************************************************************************/
/*! Private APIs */
/*************************************************************************/
/*************************************************************************/


/*************************************************************************/
/*! The SIGSEGV handler */
/*************************************************************************/
static void _sb_handler(int sig, siginfo_t *si, void *unused)
{
  sbchunk_t *sbchunk;
  size_t ip, addr;

  addr = (size_t)si->si_addr;

  //printf("Inside sbhandle for %zx, %zx, %x\n", addr, last_addr, (int)pthread_self());

  /* find the sbchunk */
  BD_GET_LOCK(&(sbinfo->mtx));
  if ((sbchunk = _sb_find((void *)addr)) == NULL) {
    printf("_sb_handler: got a SIGSEGV on an unhandled memory location: %zx\n", addr);
    for (sbchunk=sbinfo->head; sbchunk!=NULL; sbchunk=sbchunk->next)
      printf("%d %zu %zx - %zx\n", sbchunk->flags, sbchunk->nbytes, sbchunk->saddr, sbchunk->eaddr);
    abort();
    exit(EXIT_FAILURE);
  }
  BD_LET_LOCK(&(sbinfo->mtx));

  //printf("  Protection: %x\n", (unsigned int)ptr->flags);

  BD_GET_LOCK(&(sbchunk->mtx));
  /* update protection information */
  if (sbchunk->flags&SBCHUNK_NONE) {
    /* on first exception, load the data in memory */
    _sb_chunkload(sbchunk);
  }
  else if (addr == last_addr) {
    /* this happen due to a write, change the protection of the corresponding page */
    ip = (addr - sbchunk->saddr)/sbinfo->pagesize;
    sbchunk->pflags[ip] |= SBCHUNK_WRITE;
    if (mprotect((void *)(sbchunk->saddr+ip*sbinfo->pagesize), sbinfo->pagesize, PROT_READ|PROT_WRITE) == -1) {
      perror("_sb_handler: failed to mprotect to PROT_READ|PROT_WRITE");
      exit(EXIT_FAILURE);
    }

    sbchunk->flags |= SBCHUNK_WRITE;

    last_addr = 0;
  }
  else {
    last_addr = addr;
  }
  BD_LET_LOCK(&(sbchunk->mtx));

}


/*************************************************************************/
/*! Returns a pointer to an sbchunk that contains the specified address */
/*************************************************************************/
sbchunk_t *_sb_find(void *ptr)
{
  sbchunk_t *sbchunk;
  size_t addr=(size_t)ptr;

  if (sbinfo == NULL)
    return NULL;

  /* find the sbchunk */
  for (sbchunk=sbinfo->head; sbchunk!=NULL; sbchunk=sbchunk->next) {
    if (addr>=sbchunk->saddr && addr<sbchunk->eaddr)
      break;
  }

  return sbchunk;
}


/*************************************************************************/
/*! Loads the supplied sbchunk and makes it readable */
/*************************************************************************/
void _sb_chunkload(sbchunk_t *sbchunk)
{
  ssize_t ip, ifirst, size, tsize, trsize, npages;
  char *buf;
  uint8_t *pflags;
  int fd;


  /*----------------------------------------------------------------------*/
#ifdef BDMPL_WITH_SB_NOTIFY
  if (sbchunk->flags&SBCHUNK_NONE) {
    bdmsg_t msg, gomsg;

    /* notify the master that you want to load memory */
    msg.msgtype = BDMPI_MSGTYPE_MEMLOAD;
    msg.source  = sbinfo->job->rank;
    msg.count   = sbchunk->npages*sbinfo->pagesize;
    bdmq_send(sbinfo->job->reqMQ, &msg, sizeof(bdmsg_t));
    BDMPL_SLEEP(sbinfo->job, gomsg);
  }
#endif
  /*----------------------------------------------------------------------*/


  npages = sbchunk->npages;
  pflags = sbchunk->pflags;

  //printf("_sb_chunkload: loading %s\n", sbchunk->fname);

  /* if required, read the data from the file */
  if (sbchunk->flags&SBCHUNK_ONDISK) {
    trsize = 0;
    /* make it writeable for load/clear */
    if (mprotect((void *)sbchunk->saddr, npages*sbinfo->pagesize, PROT_WRITE) == -1) {
      perror("_sb_chunkload: failed to PROT_WRITE");
      exit(EXIT_FAILURE);
    }
    if ((fd = open(sbchunk->fname, O_RDONLY)) == -1) {
      perror("_sb_chunkload: failed to open file");
      exit(EXIT_FAILURE);
    }

    for (ifirst=-1, ip=0; ip<=npages; ip++) {
      if (pflags[ip]&SBCHUNK_ONDISK) {
        if (ifirst == -1)
          ifirst = ip;
      }
      else if (ifirst != -1) {
        //printf("_sb_chunkload: reading from [%zu...%zu) out of %zu.\n", ifirst, ip, npages);
        if (lseek(fd, ifirst*sbinfo->pagesize, SEEK_SET) == -1) {
          perror("_sb_chunkload: failed on lseek");
          exit(EXIT_FAILURE);
        }

        tsize = (ip-ifirst)*sbinfo->pagesize;
        buf = (char *)(sbchunk->saddr + ifirst*sbinfo->pagesize);
        trsize += tsize;
        do {
          if ((size = libc_read(fd, buf, tsize)) == -1) {
            perror("_sb_chunkload: failed to read the required data");
            exit(EXIT_FAILURE);
          }
          buf   += size;
          tsize -= size;
        } while (tsize > 0);

        ifirst = -1;
      }
    }

    if (close(fd) == -1) {
      perror("sb_realloc: failed to close the fd");
      exit(EXIT_FAILURE);
    }
    //printf("_sb_chunkload: %s, size: %zu\n", sbchunk->fname, trsize);
  }

  if (mprotect((void *)sbchunk->saddr, npages*sbinfo->pagesize, PROT_READ) == -1) {
    perror("sb_realloc: failed to PROT_READ");
    exit(EXIT_FAILURE);
  }

  sbchunk->flags ^= SBCHUNK_NONE;
  sbchunk->flags |= SBCHUNK_READ;

  for (ip=0; ip<npages; ip++) {
    pflags[ip] |= SBCHUNK_READ;
    if (pflags[ip]&SBCHUNK_NONE)
      pflags[ip] ^= SBCHUNK_NONE;
  }
}


/*************************************************************************/
/*! Saves the supplied sbchunk to disk; this is an internal routine */
/*************************************************************************/
void _sb_chunksave(sbchunk_t *sbchunk)
{
  size_t ip, ifirst, npages, size, tsize, twsize=0, count=0;
  char *buf;
  uint8_t *pflags;
  int fd;

  /*----------------------------------------------------------------------*/
#ifdef BDMPL_WITH_SB_NOTIFY
  if (!(sbchunk->flags&SBCHUNK_NONE)) {
    bdmsg_t msg, gomsg;

    /* notify the master that you have saved memory */
    msg.msgtype = BDMPI_MSGTYPE_MEMSAVE;
    msg.source  = sbinfo->job->rank;
    msg.count   = sbchunk->npages*sbinfo->pagesize;
    bdmq_send(sbinfo->job->reqMQ, &msg, sizeof(bdmsg_t));
    BDMPL_SLEEP(sbinfo->job, gomsg);
  }
#endif
  /*----------------------------------------------------------------------*/

  _sb_chunksave_internal(sbchunk);
}


/*************************************************************************/
/*! Saves the supplied sbchunk to disk; this is an internal routine */
/*************************************************************************/
size_t _sb_chunksave_internal(sbchunk_t *sbchunk)
{
  size_t ip, ifirst, npages, size, tsize, twsize=0, count=0;
  char *buf;
  uint8_t *pflags;
  int fd;

  if (!(sbchunk->flags&SBCHUNK_NONE))
    count = sbchunk->npages*sbinfo->pagesize;

  npages = sbchunk->npages;
  pflags = sbchunk->pflags;

  if (sbchunk->flags&SBCHUNK_WRITE) { /* some data have been modified */
    /* open the file for writing, and create it if it does not exist */
    if ((fd = open(sbchunk->fname, O_RDWR|O_CREAT, S_IRUSR|S_IWUSR)) == -1) {
      perror("_sb_chunksave: failed to open file for sbchunk");
      exit(EXIT_FAILURE);
    }

    /* go over the pages and write the ones that have changed.
       perform the writes in contigous chunks of changed pages. */
    for (ifirst=-1, ip=0; ip<=npages; ip++) {
      if (pflags[ip]&SBCHUNK_WRITE) {
        pflags[ip] |= SBCHUNK_ONDISK;
        if (ifirst == -1)
          ifirst = ip;
      }
      else if (ifirst != -1) {
        /* write from [ifirst...ip) */
        if (lseek(fd, ifirst*sbinfo->pagesize, SEEK_SET) == -1) {
          perror("_sb_chunksave: failed on lseek");
          exit(EXIT_FAILURE);
        }

        /* write the data */
        tsize = (ip-ifirst)*sbinfo->pagesize;
        buf = (char *)(sbchunk->saddr + ifirst*sbinfo->pagesize);
        twsize += tsize;
        do {
          if ((size = libc_write(fd, buf, tsize)) == -1) {
            perror("_sb_chunksave: failed to write the required data");
            exit(EXIT_FAILURE);
          }
          buf   += size;
          tsize -= size;
        } while (tsize > 0);

        ifirst = -1;
      }
    }

    if (close(fd) == -1) {
      perror("_sb_chunksave: failed to close the fd data for write");
      exit(EXIT_FAILURE);
    }

    sbchunk->flags |= SBCHUNK_ONDISK;

    //printf("_sb_chunksave: %s, size: %zu\n", sbchunk->fname, twsize);
  }

  /* reset flags */
  if (sbchunk->flags&SBCHUNK_WRITE)
    sbchunk->flags ^= SBCHUNK_WRITE;
  if (sbchunk->flags&SBCHUNK_READ)
    sbchunk->flags ^= SBCHUNK_READ;
  sbchunk->flags |= SBCHUNK_NONE;

  for (ip=0; ip<npages; ip++) {
    pflags[ip] |= SBCHUNK_NONE;
    if (pflags[ip]&SBCHUNK_READ)
      pflags[ip] ^= SBCHUNK_READ;
    if (pflags[ip]&SBCHUNK_WRITE)
      pflags[ip] ^= SBCHUNK_WRITE;
  }

  /* set madvise */
  if (madvise((void *)sbchunk->saddr, npages*sbinfo->pagesize, MADV_DONTNEED) == -1) {
    perror("_sb_chunksave: failed to MADV_DONTNEED");
    exit(EXIT_FAILURE);
  }

  /* unlock */
  /*if (0 != munlock((void *)sbchunk->saddr, npages*sbinfo->pagesize)) {
    perror("_sb_chunksave: failed to munlock");
    exit(EXIT_FAILURE);
  }*/

  /* change protection to PROT_NONE for next time */
  if (mprotect((void *)sbchunk->saddr, npages*sbinfo->pagesize, PROT_NONE) == -1) {
    perror("_sb_chunksave: failed to PROT_NONE");
    exit(EXIT_FAILURE);
  }

  return count;
}


/*************************************************************************/
/*! Frees an sbchunk and its anonymous mmap */
/*************************************************************************/
void _sb_chunkfree(sbchunk_t *sbchunk)
{
  /*----------------------------------------------------------------------*/
#ifdef BDMPL_WITH_SB_NOTIFY
{
  if (!(sbchunk->flags&SBCHUNK_NONE)) {
    bdmsg_t msg, gomsg;

    /* notify the master that you want to load memory */
    msg.msgtype = BDMPI_MSGTYPE_MEMSAVE;
    msg.source  = sbinfo->job->rank;
    msg.count   = sbchunk->npages*sbinfo->pagesize;
    bdmq_send(sbinfo->job->reqMQ, &msg, sizeof(bdmsg_t));
    BDMPL_SLEEP(sbinfo->job, gomsg);
  }
}
#endif
  /*----------------------------------------------------------------------*/


  /* delete the file, if on disk */
  if (sbchunk->flags&SBCHUNK_ONDISK) {
    if (unlink(sbchunk->fname) == -1) {
      perror("Failed to remove the file for sbchunk");
      exit(EXIT_FAILURE);
    }
  }

  /* unmap and remove sbchunk related malloced space */
  if (munmap((void *)sbchunk->saddr, sbchunk->npages*sbinfo->pagesize) == -1) {
    perror("Failed to unmap for for _sb_chunkfree");
    exit(EXIT_FAILURE);
  }

  libc_free(sbchunk->pflags);
  libc_free(sbchunk->fname);
  libc_free(sbchunk);

  return;
}


/*************************************************************************/
/*! Frees all the sbchunks */
/*************************************************************************/
void _sb_chunkfreeall()
{
  sbchunk_t *ptr, *nptr;

  ptr = sbinfo->head;
  while (ptr != NULL) {
    nptr = ptr->next;
    _sb_chunkfree(ptr);
    ptr = nptr;
  }
  sbinfo->head = NULL;

  return;
}
