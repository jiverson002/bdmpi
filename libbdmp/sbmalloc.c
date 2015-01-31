/*!
\file
\brief Implements the storage-backed malloc interface
\date Started 4/1/2014
\author George
*/


#include "bdmplib.h"
#include "asio.h"


/* Stores information associated with a storage-backed memory chunk */
typedef struct sbchunk {
  int sig;               /* async signal */
  size_t saddr, eaddr;   /* starting/ending address of the anonymous mapping */
  size_t npages;         /* number of pages allocated */
  size_t ldpages;        /* number of loaded pages */
  size_t nbytes;         /* number of bytes requested */
  uint8_t flags;         /* global chunk-level flag */
  uint8_t *pflags;       /* per-page flag vector */
  char *fname;           /* the file that will store the data */
  sem_t sem;             /* the async semaphore */
  pthread_mutex_t mtx;   /* the mutex guarding it */
  struct sbchunk *next;  /* pointer to the next chunk */
} sbchunk_t;

/* Stores global information associated with storage-backed memory */
typedef struct {
  int opts;              /* the sb library options */
  size_t pagesize;       /* the size of a memory page */
  size_t minsize;        /* the minimum allocation in pages handled by sbmalloc */
  char *fstem;           /* the file stem where the data is stored */
  sjob_t *job;           /* the slave job information */
  sbchunk_t *head;       /* the first sbchunk */
  pthread_mutex_t mtx;   /* the mutex guarding it */
  pthread_mutexattr_t mtx_attr;  /* attributes for all mutexes */

  struct sigaction act, oldact;  /* for the SIGSEGV signal handler */
} sbinfo_t;


/* static global variables */
static sbinfo_t *sbinfo=NULL;
static struct asio_env *sbasio=NULL;


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


/* async I/O prototypes */
static int _sb_init_sbchunk(sbchunk_t * const sbchunk);
static int _sb_quit_sbchunk(sbchunk_t * const sbchunk);
static int _sb_free_sbchunk(sbchunk_t * const sbchunk);
static void _sb_asio_cancel(void * const arg);
static void _sb_asio_cb(void * const arg);
void _sb_chunkload_mt(sbchunk_t * const sbchunk, size_t const ip);

/* private prototypes */
sbchunk_t *_sb_find(void *ptr);
static void _sb_handler(int sig, siginfo_t *si, void *unused);
void _sb_chunkload(sbchunk_t *sbchunk);
size_t _sb_chunksave(sbchunk_t *sbchunk, int const flag);
void _sb_chunkfree(sbchunk_t *sbchunk);
void _sb_chunkfreeall();

void _sb_pageload(sbchunk_t * const sbchunk, size_t const ip);
void _sb_pagesave(sbchunk_t * const sbchunk, size_t const ip);


#ifdef XXX
/* macros */
#define BD_GET_LOCK(lock) GKASSERT(pthread_mutex_lock(lock) == 0)
#define BD_LET_LOCK(lock) GKASSERT(pthread_mutex_unlock(lock) == 0)
#endif

#define MLOCK(SADDR, SIZE)                                                  \
do {                                                                        \
  if (-1 == mlock((void*)(SADDR), SIZE)) {                                  \
    fprintf(stderr, "%s:%d: failed to mlock: %s\n", __func__, __LINE__,     \
      strerror(errno));                                                     \
    exit(EXIT_FAILURE);                                                     \
  }                                                                         \
} while (0)

#define MUNLOCK(SADDR, SIZE)                                                \
do {                                                                        \
  if (-1 == munlock((void*)(SADDR), SIZE)) {                                \
    fprintf(stderr, "%s:%d: failed to mlock: %s\n", __func__, __LINE__,     \
      strerror(errno));                                                     \
    exit(EXIT_FAILURE);                                                     \
  }                                                                         \
} while (0)

#define MPROTECT(SADDR, PGSIZE, PROT)                                       \
do {                                                                        \
  /*MPROTECT((sbchunk->saddr+ip*sbinfo->pagesize), sbinfo->pagesize,  */    \
  /*  PROT_READ|PROT_WRITE);                                          */    \
  if (-1 == mprotect((void*)(SADDR), PGSIZE, PROT)) {                       \
    fprintf(stderr, "%s:%d: failed to mprotect: %s\n", __func__, __LINE__,  \
      strerror(errno));                                                     \
    /*fprintf(stderr, "%zu %zu %zu %zu\n", (size_t)sbchunk->saddr,    */    \
    /*  (size_t)sbchunk->eaddr, (size_t)(SADDR), PGSIZE);             */    \
    exit(EXIT_FAILURE);                                                     \
  }                                                                         \
} while (0)


//#define __cplusplus
#ifdef __cplusplus
void *dlmalloc(size_t nbytes);
void *dlrealloc(void *ptr, size_t nbytes);
void *dlfree(void *ptr);
#endif


/*************************************************************************/
/*! Hook: malloc */
/*************************************************************************/
void *malloc(size_t nbytes)
{
  if (libc_malloc == NULL)
    *((void **) &libc_malloc) = dlsym(RTLD_NEXT, "malloc");
  // should do error check to make sure symbol was found?
  //if (libc_malloc == NULL)
  //  return NULL;

  //printf("%zu\n", nbytes);

  if (sbinfo == NULL)
    return libc_malloc(nbytes);
  else {
#ifndef __cplusplus
    if (sbinfo->minsize == 0 || nbytes <= sbinfo->minsize)
      return libc_malloc(nbytes);
    else
      return sb_malloc(nbytes);
#else
      return dlmalloc(nbytes);
#endif
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
#ifndef __cplusplus
      return sb_realloc(ptr, size);
#else
      return dlrealloc(ptr, size);
#endif
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
#ifndef __cplusplus
      sb_free(ptr);
#else
      dlfree(ptr);
#endif
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
int sb_init(char *fstem, sjob_t * const job)
{
  int i, j;

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

  sbinfo->job      = job;
  sbinfo->opts     = job->jdesc->sbopts;
  sbinfo->minsize  = job->jdesc->sbsize*sysconf(_SC_PAGESIZE);
  sbinfo->pagesize = job->jdesc->pgsize*sysconf(_SC_PAGESIZE);
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


  /***************************************/
  /* setup the multi-threaded I/O struct */
  /***************************************/
  SB_SB_IFSET(BDMPI_SB_ASIO) {
    if (sbasio != NULL) {
      perror("sb_init: sbasio != NULL");
      exit(EXIT_FAILURE);
    }

    if (NULL == (sbasio=libc_malloc(sizeof(struct asio_env)))) {
      perror("Failed to allocate sbasio\n");
      goto ERROR_EXIT;
    }
    memset(sbasio, 0, sizeof(struct asio_env));

    if (0 != asio_init(sbasio, 2, 64, &_sb_asio_cb)) {
      fprintf(stderr, "Failed to initialize asio.\n");
      goto ERROR_EXIT;
    }
  }

  return 1;

ERROR_EXIT:
  if (sbinfo->fstem)
    libc_free(sbinfo->fstem);

  libc_free(sbinfo);
  sbinfo = NULL;

  SB_SB_IFSET(BDMPI_SB_ASIO) {
    if (NULL != sbasio) {
      (void)asio_free(sbasio);
      libc_free(sbasio);
      sbasio = NULL;
    }
  }

  return 0;
}


/*************************************************************************/
/*! Shutsdown the sbmalloc subsystem */
/*************************************************************************/
int sb_finalize()
{
  int i;

  if (sbinfo == NULL) {
    perror("sbfinalize: sbinfo -= NULL");
    exit(EXIT_FAILURE);
  }

  if (sigaction(SIGSEGV, &(sbinfo->oldact), NULL) == -1) {
    perror("Failed to install the old signal handler\n");
    return 0;
  }

  _sb_chunkfreeall();

  SB_SB_IFSET(BDMPI_SB_ASIO) {
    if (NULL == sbasio) {
      perror("sbfinalize: sbasio -= NULL");
      exit(EXIT_FAILURE);
    }

    if (0 != asio_free(sbasio)) {
      fprintf(stderr, "Could not destroy sbasio\n");
      return 0;
    }
    libc_free(sbasio);
    sbasio = NULL;
  }

  GKASSERT(pthread_mutex_destroy(&(sbinfo->mtx)) == 0);
  GKASSERT(pthread_mutexattr_destroy(&(sbinfo->mtx_attr)) == 0);

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

  if (sbinfo == NULL) {
    perror("sb_malloc: sbinfo == NULL");
    exit(EXIT_FAILURE);
  }
  if (NULL == libc_malloc) {
    perror("sb_malloc: libc_malloc == NULL");
    exit(EXIT_FAILURE);
  }

  sbchunk = libc_malloc(sizeof(sbchunk_t));
  if (sbchunk == NULL)
    return NULL;

  sbchunk->sig     = 0;
  sbchunk->ldpages = 0;
  sbchunk->nbytes  = nbytes;

  /* determine the allocation size in terms of pagesize */
  sbchunk->npages = (nbytes+sbinfo->pagesize-1)/sbinfo->pagesize;

  /* allocate the flag array for the pages */
  sbchunk->pflags = (uint8_t *)libc_malloc((sbchunk->npages+1)*sizeof(uint8_t));
  if (sbchunk->pflags == NULL) {
    fprintf(stderr, "[%zu, %zu] ", nbytes, sbchunk->npages+1);
    perror("sbmalloc.1 failure");
    libc_free(sbchunk);
    return NULL;
  }

  sbchunk->saddr = (size_t) mmap(NULL, sbchunk->npages*sbinfo->pagesize,  \
    PROT_NONE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);

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

  /* updated sb permissions */
  sbchunk->flags = SBCHUNK_NONE;

  for (ip=0; ip<sbchunk->npages; ++ip)
    sbchunk->pflags[ip] = SBCHUNK_NONE;
  sbchunk->pflags[sbchunk->npages] = 0;

  /* Initialize sbchunk for asio */
  GKASSERT(0 == _sb_init_sbchunk(sbchunk));

  BD_GET_LOCK(&(sbinfo->mtx));
  sbchunk->next = sbinfo->head;
  sbinfo->head  = sbchunk;
  BD_LET_LOCK(&(sbinfo->mtx));

  //printf("MALLOC (%p) %zu\n", (void*)sbchunk, sbchunk->npages);

  return (void *)sbchunk->saddr;
}


/*************************************************************************/
/*! Reallocates memory via anonymous mmap */
/*************************************************************************/
void *sb_realloc(void *oldptr, size_t nbytes)
{
  size_t ip, npages, new_saddr, new_npages, count;
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

  npages = sbchunk->npages;

  /* see if we are shrinking */
  if (nbytes <= sbchunk->nbytes) { /* easy case */
    for (count=0,ip=new_npages; ip<sbchunk->npages; ++ip) {
      if (0 == (sbchunk->pflags[ip]&SBCHUNK_NONE))
        count++;
    }
    GKASSERT(sbchunk->ldpages >= count);
    sbchunk->ldpages -= count;

    /* set the now unused pages as DONTNEED */
    if (madvise((void *)(sbchunk->saddr+new_npages*sbinfo->pagesize),
          (sbchunk->npages-new_npages)*sbinfo->pagesize, MADV_DONTNEED) == -1)
    {
      perror("sb_realloc: failed to MADV_DONTNEED");
      goto ERROR_EXIT;
    }

    sbchunk->nbytes = nbytes;
    sbchunk->npages = new_npages;
    sbchunk->eaddr  = sbchunk->saddr+nbytes;
    sbchunk->pflags[new_npages] = 0;  /* this is for the +1 reset */

    /*----------------------------------------------------------------------*/
    SB_SB_IFSET(BDMPI_SB_LAZYWRITE) {
      if (0 == (sbchunk->flags&SBCHUNK_NONE)) {
        bdmsg_t msg, gomsg;

        /* notify the master that you want to save memory */
        memset(&msg, 0, sizeof(bdmsg_t));
        msg.msgtype = BDMPI_MSGTYPE_MEMSAVE;
        msg.source  = sbinfo->job->rank;
        msg.count   = (npages-new_npages)*sbinfo->pagesize;
        bdmq_send(sbinfo->job->reqMQ, &msg, sizeof(bdmsg_t));
        BDMPL_SLEEP(sbinfo->job, gomsg);
        //printf("[%d, %d] (%p: -%zu)\n", getpid(), __LINE__, (void*)sbchunk,
        //  msg.count);
        //fflush(stdout);
      }
    }
    /*----------------------------------------------------------------------*/
  }
  else {
    /*----------------------------------------------------------------------*/
    SB_SB_IFSET(BDMPI_SB_LAZYWRITE) {
      bdmsg_t msg, gomsg;

      /* notify the master that you want to load memory */
      memset(&msg, 0, sizeof(bdmsg_t));
      msg.msgtype = BDMPI_MSGTYPE_MEMLOAD;
      msg.source  = sbinfo->job->rank;
      msg.count   = (new_npages-sbchunk->npages)*sbinfo->pagesize;
      bdmq_send(sbinfo->job->reqMQ, &msg, sizeof(bdmsg_t));
      BDMPL_SLEEP(sbinfo->job, gomsg);
      //printf("[%d, %d] (%p: +%zu)\n", getpid(), __LINE__, (void*)sbchunk,
      //  msg.count);
      //fflush(stdout);
    }
    /*----------------------------------------------------------------------*/

    /* allocate the new pflags */
    if ((new_pflags = (uint8_t *)libc_malloc((new_npages+1)*sizeof(uint8_t))) == NULL) {
      perror("sb_realloc: failed to malloc new_pflags");
      goto ERROR_EXIT;
    }

    /* get the new anonymous map */
    new_saddr = (size_t) mmap(NULL, new_npages*sbinfo->pagesize, PROT_WRITE,
                             MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);

    if ((void *)new_saddr == MAP_FAILED) {
      perror("sb_realloc: failed to get a new mmap");
      goto ERROR_EXIT;
    }

    /* Load the sbchunk if not accessible. */
    if (sbchunk->ldpages < sbchunk->npages)
      _sb_chunkload(sbchunk);

    /* error check */
    GKASSERT(sbchunk->flags&SBCHUNK_READ);
    for (ip=0; ip<sbchunk->npages; ip++)
      GKASSERT(sbchunk->pflags[ip]&SBCHUNK_READ);

    /* copy the data */
    memcpy((void *)new_saddr, (void *)sbchunk->saddr, sbchunk->nbytes);

    /* unmap the old */
    if (munmap((void *)sbchunk->saddr, sbchunk->npages*sbinfo->pagesize) == -1) {
      perror("sb_realloc: failed to unmap the old region");
      exit(EXIT_FAILURE);
    }

    /* set permissions to read only */
    MPROTECT(new_saddr, new_npages*sbinfo->pagesize, PROT_READ);

    /* set remapped sbchunk for writing and remove its ondisk flag */
    if (sbchunk->flags&SBCHUNK_ONDISK) {
      sbchunk->flags ^= SBCHUNK_ONDISK;
      sbchunk->flags |= SBCHUNK_WRITE;
    }

    /* set per-page flags */
    for (ip=0; ip<sbchunk->npages; ip++) {
      new_pflags[ip] = sbchunk->pflags[ip];
      if (SBCHUNK_WRITE == (new_pflags[ip]&SBCHUNK_WRITE) ||
          SBCHUNK_ONDISK == (new_pflags[ip]&SBCHUNK_ONDISK))
      {
        MPROTECT(new_saddr+ip*sbinfo->pagesize, sbinfo->pagesize, \
          PROT_READ|PROT_WRITE);

        if (SBCHUNK_ONDISK == (new_pflags[ip]&SBCHUNK_ONDISK))
          new_pflags[ip] ^= SBCHUNK_ONDISK;
        new_pflags[ip] |= SBCHUNK_WRITE;
      }
    }
    for (; ip<new_npages; ip++)
      new_pflags[ip] = SBCHUNK_READ;
    new_pflags[new_npages] = 0;

    libc_free(sbchunk->pflags);

    sbchunk->nbytes  = nbytes;
    sbchunk->npages  = new_npages;
    sbchunk->ldpages = new_npages;
    sbchunk->saddr   = new_saddr;
    sbchunk->eaddr   = sbchunk->saddr+nbytes;
    sbchunk->pflags  = new_pflags;
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
  if (NULL == libc_free) {
    perror("sb_free: libc_free == NULL");
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

  /* Initialize sbchunk for asio */
  GKASSERT(0 == _sb_free_sbchunk(ptr));

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
  _sb_chunksave(sbchunk, 1);
  BD_LET_LOCK(&(sbchunk->mtx));
}


/*************************************************************************/
/*! Saves to disk all sbchunks */
/*************************************************************************/
void sb_saveall()
{
  sbchunk_t *sbchunk;

  if (sbinfo == NULL)
    return;

  BD_GET_LOCK(&(sbinfo->mtx));
  for (sbchunk=sbinfo->head; sbchunk!=NULL; sbchunk=sbchunk->next) {
    BD_GET_LOCK(&(sbchunk->mtx));
    (void)_sb_chunksave(sbchunk, 1);
    BD_LET_LOCK(&(sbchunk->mtx));
  }
  BD_LET_LOCK(&(sbinfo->mtx));
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
    count += _sb_chunksave(sbchunk, 0);
    BD_LET_LOCK(&(sbchunk->mtx));
  }
  BD_LET_LOCK(&(sbinfo->mtx));

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
  if (sbchunk->ldpages < sbchunk->npages)
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
    if (sbchunk->ldpages < sbchunk->npages)
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
#if 0
  size_t addr, ip, ifirst, iend, count;
  sbchunk_t *sbchunk;

  if (NULL == sbinfo)
    return;
  if (0 == size)
    return;
  if (size < -1)
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
    /* can only discard pages fully within range, thus ifirst is a ceil
     * operation and iend is a floor operation. */
    ifirst = (addr == sbchunk->saddr) ? 0 : 1+((addr-sbchunk->saddr-1)/sbinfo->pagesize); /* ceiling */
    iend   = (addr+size == sbchunk->saddr+sbchunk->nbytes) ? sbchunk->npages :
      (addr+size-sbchunk->saddr)/sbinfo->pagesize;  /* floor */
  }

  if (ifirst < iend) {
    for (count=0,ip=ifirst; ip<iend; ++ip) {
      if (0 == (sbchunk->pflags[ip]&SBCHUNK_NONE))
        count++;
      else {
        GKASSERT(0 == (sbchunk->pflags[ip]&SBCHUNK_READ));
        GKASSERT(0 == (sbchunk->pflags[ip]&SBCHUNK_WRITE));
      }
    }

    if (0 != count) {
      /*------------------------------------------------------------------*/
      SB_SB_IFSET(BDMPI_SB_LAZYWRITE) {
        //printf("DISCARD (%p) %zu / %zu / %zu\n", (void*)sbchunk, count,
        //  sbchunk->ldpages, sbchunk->npages);
#if 1
        bdmsg_t msg, gomsg;

        /* notify the master that you want to save memory */
        memset(&msg, 0, sizeof(bdmsg_t));
        msg.msgtype = BDMPI_MSGTYPE_MEMSAVE;
        msg.source  = sbinfo->job->rank;
        msg.count   = count*sbinfo->pagesize;
        bdmq_send(sbinfo->job->reqMQ, &msg, sizeof(bdmsg_t));
        BD_LET_LOCK(&(sbchunk->mtx));
        BDMPL_SLEEP(sbinfo->job, gomsg);
        BD_GET_LOCK(&(sbchunk->mtx));
#endif
      }
      /*------------------------------------------------------------------*/

#if 1
      /* provide no permissions this will remove the write permissions so in
       * case we do not block, subsequent writes will be intercepted correctly
       * */
      MPROTECT(sbchunk->saddr+ifirst*sbinfo->pagesize,  \
        (iend-ifirst)*sbinfo->pagesize, PROT_NONE);

      /* update the corresponding pflags[] entries */
      for (ip=ifirst; ip<iend; ip++) {
        if (sbchunk->pflags[ip]&SBCHUNK_READ)
          sbchunk->pflags[ip] ^= SBCHUNK_READ;
        if (sbchunk->pflags[ip]&SBCHUNK_WRITE)
          sbchunk->pflags[ip] ^= SBCHUNK_WRITE;
        if (sbchunk->pflags[ip]&SBCHUNK_ONDISK)
          sbchunk->pflags[ip] ^= SBCHUNK_ONDISK;
        sbchunk->pflags[ip] |= SBCHUNK_NONE;
      }
      sbchunk->ldpages -= count;

      if (0 == sbchunk->ldpages) {
        if (sbchunk->flags&SBCHUNK_READ)
          sbchunk->flags ^= SBCHUNK_READ;
        if (sbchunk->flags&SBCHUNK_WRITE)
          sbchunk->flags ^= SBCHUNK_WRITE;
        sbchunk->flags |= SBCHUNK_NONE;
      }
#endif
    }
  }
  BD_LET_LOCK(&(sbchunk->mtx));
#endif
}


/*************************************************************************/
/*************************************************************************/
/*! Private APIs */
/*************************************************************************/
/*************************************************************************/


/*************************************************************************/
/*! Initialize a sbchunk for asio */
/*************************************************************************/
static int _sb_init_sbchunk(sbchunk_t * const sbchunk)
{
  sbchunk->sig = 0;
  if (0 != sem_init(&(sbchunk->sem), 0, 1))
    goto CLEANUP;
  if (0 != pthread_mutex_init(&(sbchunk->mtx), NULL))
    goto CLEANUP;

  return 0;

CLEANUP:
  return -1;
}


/*************************************************************************/
/*! Tell async handlers to release an sbchunk */
/*************************************************************************/
static int _sb_quit_sbchunk(sbchunk_t * const sbchunk)
{
#ifndef NDEBUG
  int sval;
#endif
  int hassem=0;

  BD_LET_LOCK(&(sbchunk->mtx));

  /* Check if chunk is async I/O, 0==hassem:yes, 1==hassem:no.  If chunk is
   * async I/O, signal any async threads to quit. */
  BD_TRY_SEM(&(sbchunk->sem), hassem);
  if (0 == hassem) {
    /* Signal that the async handler of this piece of work, if any, should
     * break. */
    BD_GET_LOCK(&(sbchunk->mtx));
    assert(0 == sbchunk->sig);
    sbchunk->sig = SIGQUIT;
    BD_LET_LOCK(&(sbchunk->mtx));

    /* Wait until async handler has finished */
    BD_GET_SEM(&(sbchunk->sem));
  }

  BD_GET_LOCK(&(sbchunk->mtx));
  BD_LET_SEM(&(sbchunk->sem));

  return 0;
}


/*************************************************************************/
/*! Destroy a sbchunk for asio */
/*************************************************************************/
static int _sb_free_sbchunk(sbchunk_t * const sbchunk)
{
#ifndef NDEBUG
  int sval;
#endif
  int hassem=0;

  /* Check if chunk is async I/O, 0==hassem:yes, 1==hassem:no.  If chunk is
   * async I/O, signal any async threads to quit. */
  BD_TRY_SEM(&(sbchunk->sem), hassem);
  if (0 == hassem) {
    /* Signal that the async handler of this piece of work, if any, should
     * break. */
    BD_GET_LOCK(&(sbchunk->mtx));
    assert(0 == sbchunk->sig);
    sbchunk->sig = SIGQUIT;
    BD_LET_LOCK(&(sbchunk->mtx));

    /* Wait until async handler has finished */
    BD_GET_SEM(&(sbchunk->sem));
  }

  /* Free work -- wait until lock is re-acquired to be sure that async thread
   * is finished. */
  BD_GET_LOCK(&(sbchunk->mtx));
  assert(0 == sem_getvalue(&(sbchunk->sem), &sval));
  assert(0 == sval);
  BD_LET_LOCK(&(sbchunk->mtx));

  if (0 != sem_destroy(&(sbchunk->sem)))
    goto CLEANUP;
  if (0 != pthread_mutex_destroy(&(sbchunk->mtx)))
    goto CLEANUP;

  return 0;

CLEANUP:
  return -1;
}


/*************************************************************************/
/*! The pthread cancel handler for _sb_asio_cb() */
/*************************************************************************/
static void _sb_asio_cancel(void * const arg)
{
  sbchunk_t * const sbchunk=(sbchunk_t *)arg;

  BD_LET_SEM(&(sbchunk->sem));
  BD_LET_LOCK(&(sbchunk->mtx));
}


/*************************************************************************/
/*! The asio callback function */
/*************************************************************************/
static void _sb_asio_cb(void * const arg)
{
#ifndef NDEBUG
  int sval;
#endif
  int fd=-1;
  size_t saddr, pgsize;
  ssize_t i, ip, iend, ifirst, size, tsize, npages;
  struct timespec ts;
  char * buf;
  uint8_t * pflags;
  sbchunk_t * sbchunk;

  if (NULL == arg) {
    fprintf(stderr, "Error: Callback received invalid arg at line %d of " \
      "file %s.\n", __LINE__, __FILE__);
    abort();
  }
  sbchunk = (sbchunk_t *)arg;

  pthread_cleanup_push(_sb_asio_cancel, sbchunk);
  BD_GET_LOCK(&(sbchunk->mtx));

  assert(0 == sem_getvalue(&(sbchunk->sem), &sval));
  assert(0 == sval);
  assert(SBCHUNK_NONE != (sbchunk->flags&SBCHUNK_NONE));

  if (SBCHUNK_ONDISK == (sbchunk->flags&SBCHUNK_ONDISK)) {
    if (-1 == (fd=open(sbchunk->fname, O_RDONLY))) {
      perror("_sb_asio_cb: failed to open file");
      exit(EXIT_FAILURE);
    }
  }

  for (i=0; i<sbchunk->npages; i+=32) {
    /***************************************************/
    /* Handle any received ``signals'' */
    /***************************************************/
    BD_LET_LOCK(&(sbchunk->mtx));

    /* Time-out so that other threads can lock mutex and update sbchunk->sig
     * if need be. */
    ts.tv_sec  = 0;
    ts.tv_nsec = 5000000;
    nanosleep(&ts, NULL);

    BD_GET_LOCK(&(sbchunk->mtx));

    switch (sbchunk->sig) {
    case 0:
      break;
    case SIGQUIT:
      goto QUIT;
    default:
      fprintf(stderr, "Error: Callback received invalid signal (%d) at "  \
        "line %d of file %s.\n", sbchunk->sig, __LINE__, __FILE__);
      abort();
    }

    /***************************************************/
    /* Do some background work */
    /***************************************************/
    iend   = i+32 < sbchunk->npages ? i+32 : sbchunk->npages;
    saddr  = sbchunk->saddr;
    pflags = sbchunk->pflags;
    pgsize = sbinfo->pagesize;

    /* Incase the chunk was realloc'd while the mutex was released */
    if (i >= sbchunk->npages)
      break;

#if 0
    /* if required, read the data from the file */
    if (SBCHUNK_ONDISK == (sbchunk->flags&SBCHUNK_ONDISK)) {
      /***************************************************/
      /* read any required pages */
      /***************************************************/
      MPROTECT(saddr+(i*pgsize), (iend-i)*pgsize, PROT_WRITE);

      for (ifirst=-1, ip=i; ip<=iend; ip++) {
        if (0 == (pflags[ip]&SBCHUNK_READ)  &&
            0 == (pflags[ip]&SBCHUNK_WRITE) && /* uneeded? */
            SBCHUNK_ONDISK == (pflags[ip]&SBCHUNK_ONDISK))
        {
          if (-1 == ifirst)
            ifirst = ip;
        }
        else if (-1 != ifirst) {
          if (-1 == lseek(fd, ifirst*pgsize, SEEK_SET)) {
            perror("_sb_asio_cb: failed on lseek");
            exit(EXIT_FAILURE);
          }

          tsize = (ip-ifirst)*pgsize;
          buf = (char *)(saddr+(ifirst*pgsize));
          do {
            if (-1 == (size=libc_read(fd, buf, tsize))) {
              perror("_sb_asio_cb: failed to read the required data");
              exit(EXIT_FAILURE);
            }
            buf   += size;
            tsize -= size;
          } while (tsize > 0);

          ifirst = -1;
        }
      }
    }
#endif

#if 1
    if (SBCHUNK_ONDISK != (sbchunk->flags&SBCHUNK_ONDISK)) {
    /***************************************************/
    /* give final protection and set appropriate flags */
    /***************************************************/
    MPROTECT(saddr+(i*pgsize), (iend-i)*pgsize, PROT_READ);

    for (ip=i; ip<iend; ip++) {
      if (SBCHUNK_NONE == (pflags[ip]&SBCHUNK_NONE)) {
        pflags[ip] ^= SBCHUNK_NONE;
        sbchunk->ldpages++;
      }
      else if (SBCHUNK_WRITE == (pflags[ip]&SBCHUNK_WRITE)) {
        MPROTECT(saddr+(ip*pgsize), pgsize, PROT_READ|PROT_WRITE);
      }
      pflags[ip] |= SBCHUNK_READ;
    }
    }
#endif
  }

QUIT:
  if (-1 != fd && -1 == close(fd)) {
    perror("_sb_asio_cb: failed to close the fd");
    exit(EXIT_FAILURE);
  }

  pthread_cleanup_pop(1);
}


/*************************************************************************/
/*! The SIGSEGV handler */
/*************************************************************************/
static void _sb_handler(int sig, siginfo_t *si, void *unused)
{
  size_t ip=0;
  size_t const addr=(size_t)si->si_addr;
  sbchunk_t * sbchunk=NULL;

  /* find the sbchunk */
  BD_GET_LOCK(&(sbinfo->mtx));
  if (NULL == (sbchunk=_sb_find((void*)addr))) {
    printf("_sb_handler: got a SIGSEGV on an unhandled memory location: " \
      "%zx\n", addr);
    abort();
  }
  BD_LET_LOCK(&(sbinfo->mtx));

  /* update protection information */
  BD_GET_LOCK(&(sbchunk->mtx));
  ip = (addr-sbchunk->saddr)/sbinfo->pagesize;

  /* Only look at SBCHUNK_NONE, SBCHUNK_READ, and SBCHUNK_WRITE bits */
  switch (sbchunk->pflags[ip]&0x7) {
  case SBCHUNK_NONE:
    SB_SB_IFSET(BDMPI_SB_LAZYREAD) {
      /* on first exception, load the data page into memory */
      _sb_pageload(sbchunk, ip);
      //_sb_chunkload(sbchunk);
    }
    else SB_SB_IFSET(BDMPI_SB_ASIO) {
      /* on first exception, load the data page and setup chunk for async
       * I/O */
      _sb_chunkload_mt(sbchunk, ip);
    }
    else {
      /* on first exception, load the data chunk into memory */
      _sb_chunkload(sbchunk);
    }
    break;

  case SBCHUNK_READ:
    /* on second exception, change the protection of the page to writeable */
    MPROTECT(sbchunk->saddr+ip*sbinfo->pagesize, sbinfo->pagesize,        \
      PROT_READ|PROT_WRITE);

    sbchunk->flags      |= SBCHUNK_WRITE;
    sbchunk->pflags[ip] |= SBCHUNK_WRITE;
    break;

  case SBCHUNK_READ|SBCHUNK_WRITE:
    break;

  default:
    fprintf(stderr, "_sb_handler: got a SIGSEGV at memory location with " \
      "unknown sb-permissions: %zx (%d)\n", addr, sbchunk->flags);
    abort();
  }

#if 0
  if (SBCHUNK_NONE == (sbchunk->pflags[ip]&SBCHUNK_NONE)) {
    SB_SB_IFSET(BDMPI_SB_LAZYREAD) {
      /* on first exception, load the data page into memory */
      _sb_pageload(sbchunk, ip);
    }
    else SB_SB_IFSET(BDMPI_SB_ASIO) {
      /* on first exception, load the data page and setup chunk for async
       * I/O */
      _sb_chunkload_mt(sbchunk, ip);
    }
    else {
      /* on first exception, load the data chunk into memory */
      _sb_chunkload(sbchunk);
    }
    GKASSERT(SBCHUNK_WRITE != (sbchunk->pflags[ip]&SBCHUNK_WRITE));
  }
  else if (SBCHUNK_READ == (sbchunk->pflags[ip]&SBCHUNK_READ)) {
    if (SBCHUNK_WRITE == (sbchunk->pflags[ip]&SBCHUNK_WRITE)) {
      fprintf(stderr, "_sb_handler: got a SIGSEGV at a sb-readable and "  \
        "sb-writeable location: %zx (%d)\n", addr, sbchunk->pflags[ip]);
      abort();
    }

    /* on second exception, change the protection of the page to writeable */
    MPROTECT(sbchunk->saddr+ip*sbinfo->pagesize, sbinfo->pagesize,        \
      PROT_READ|PROT_WRITE);

    sbchunk->flags      |= SBCHUNK_WRITE;
    sbchunk->pflags[ip] |= SBCHUNK_WRITE;
  }
  else {
    fprintf(stderr, "_sb_handler: got a SIGSEGV at memory location with " \
      "unknown sb-permissions: %zx (%d)\n", addr, sbchunk->flags);
    abort();
  }
#endif

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
  size_t count;
  ssize_t ip, ifirst, size, tsize, npages;
  char *buf;
  uint8_t *pflags;
  int fd;

  SB_SB_IFSET(BDMPI_SB_LAZYWRITE) {
    if (SBCHUNK_NONE == (sbchunk->flags&SBCHUNK_NONE)) {
      bdmsg_t msg, gomsg;

      /* notify the master that you want to load memory */
      memset(&msg, 0, sizeof(bdmsg_t));
      msg.msgtype = BDMPI_MSGTYPE_MEMLOAD;
      msg.source  = sbinfo->job->rank;
      msg.count   = sbchunk->npages*sbinfo->pagesize;
      bdmq_send(sbinfo->job->reqMQ, &msg, sizeof(bdmsg_t));
      BDMPL_SLEEP(sbinfo->job, gomsg);
      //printf("[%d, %d] (%p: +%zu)\n", getpid(), __LINE__, (void*)sbchunk,
      //  msg.count);
      //fflush(stdout);
    }
  }

  for (count=0,ip=0; ip<sbchunk->npages; ++ip) {
    if (SBCHUNK_NONE == (sbchunk->pflags[ip]&SBCHUNK_NONE))
      count++;
  }
  GKASSERT(sbchunk->ldpages+count == sbchunk->npages);

  npages = sbchunk->npages;
  pflags = sbchunk->pflags;

  /* if required, read the data from the file */
  if (sbchunk->flags&SBCHUNK_ONDISK) {
    /* make it writeable for load/clear */
    MPROTECT(sbchunk->saddr, npages*sbinfo->pagesize, PROT_WRITE);

    if ((fd = open(sbchunk->fname, O_RDONLY)) == -1) {
      perror("_sb_chunkload: failed to open file");
      exit(EXIT_FAILURE);
    }

    for (ifirst=-1, ip=0; ip<=npages; ip++) {
      if (0 == (pflags[ip]&SBCHUNK_READ)  &&
          0 == (pflags[ip]&SBCHUNK_WRITE) &&
          SBCHUNK_ONDISK == (pflags[ip]&SBCHUNK_ONDISK))
      {
        if (ifirst == -1)
          ifirst = ip;
      }
      else if (ifirst != -1) {
        if (lseek(fd, ifirst*sbinfo->pagesize, SEEK_SET) == -1) {
          perror("_sb_chunkload: failed on lseek");
          exit(EXIT_FAILURE);
        }

        tsize = (ip-ifirst)*sbinfo->pagesize;
        buf = (char *)(sbchunk->saddr + ifirst*sbinfo->pagesize);
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
  }

  MPROTECT(sbchunk->saddr, npages*sbinfo->pagesize, PROT_READ);
  //MLOCK(sbchunk->saddr, npages*sbinfo->pagesize);

  if (sbchunk->flags&SBCHUNK_NONE)
    sbchunk->flags ^= SBCHUNK_NONE;
  sbchunk->flags |= SBCHUNK_READ;

  for (ip=0; ip<npages; ip++) {
    if (SBCHUNK_NONE == (pflags[ip]&SBCHUNK_NONE))
      pflags[ip] ^= SBCHUNK_NONE;
    else if (SBCHUNK_WRITE == (pflags[ip]&SBCHUNK_WRITE))
      MPROTECT(sbchunk->saddr+ip*sbinfo->pagesize, sbinfo->pagesize,  \
        PROT_READ|PROT_WRITE);
    pflags[ip] |= SBCHUNK_READ;
  }

  sbchunk->ldpages = sbchunk->npages;

  //printf("LOAD (%p) %zu / %zu\n", (void*)sbchunk, sbchunk->ldpages,
  //  sbchunk->npages);
}


/*************************************************************************/
/*! Loads the supplied sbchunk and makes it readable via multi-threading */
/*************************************************************************/
void _sb_chunkload_mt(sbchunk_t * const sbchunk, size_t const ip)
{
  int hassem=0;

  /* If page is still not read, read it. */
  if (SBCHUNK_NONE == (sbchunk->pflags[ip]&SBCHUNK_NONE))
    _sb_pageload(sbchunk, ip);

  /* Check if chunk is async I/O, 0==hassem:yes, 1==hassem:no. */
  BD_TRY_SEM(&(sbchunk->sem), hassem);

  /* Release sbchunk lock. */
  BD_LET_LOCK(&(sbchunk->mtx));

  /* If sbchunk is not already async I/O, make it so. */
  if (1 == hassem) {
    GKASSERT(0 == asio_addw(sbasio, sbchunk));
  }

  /* Re-acquire sbchunk lock. */
  BD_GET_LOCK(&(sbchunk->mtx));
}


/*************************************************************************/
/*! Saves the supplied sbchunk to disk; this is an internal routine */
/*************************************************************************/
size_t _sb_chunksave(sbchunk_t *sbchunk, int const flag)
{
  size_t ip, ifirst, npages, size, tsize;
  char *buf;
  uint8_t *pflags;
  int fd;

  if (SBCHUNK_NONE == (sbchunk->flags&SBCHUNK_NONE))
    return 0;

  SB_SB_IFSET(BDMPI_SB_ASIO) {
    _sb_quit_sbchunk(sbchunk);
  }

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
  }

  /* reset flags */
  sbchunk->flags |= SBCHUNK_NONE;
  if (sbchunk->flags&SBCHUNK_WRITE)
    sbchunk->flags ^= SBCHUNK_WRITE;
  if (sbchunk->flags&SBCHUNK_READ)
    sbchunk->flags ^= SBCHUNK_READ;

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

  /* change protection to PROT_NONE for next time */
  //MUNLOCK(sbchunk->saddr, npages*sbinfo->pagesize);
  MPROTECT(sbchunk->saddr, npages*sbinfo->pagesize, PROT_NONE);

  sbchunk->ldpages = 0;

  /*----------------------------------------------------------------------*/
  SB_SB_IFSET(BDMPI_SB_LAZYWRITE) {
    if (1 == flag) {
      bdmsg_t msg, gomsg;

      /* notify the master that you have saved memory */
      memset(&msg, 0, sizeof(bdmsg_t));
      msg.msgtype = BDMPI_MSGTYPE_MEMSAVE;
      msg.source  = sbinfo->job->rank;
      msg.count   = sbchunk->npages*sbinfo->pagesize;
      bdmq_send(sbinfo->job->reqMQ, &msg, sizeof(bdmsg_t));
      BDMPL_SLEEP(sbinfo->job, gomsg);
    }
  }
  /*----------------------------------------------------------------------*/

  return sbchunk->npages*sbinfo->pagesize;
}


/*************************************************************************/
/*! Loads the supplied ip from sbchunk and makes it readable */
/*************************************************************************/
void _sb_pageload(sbchunk_t * const sbchunk, size_t const ip)
{
  ssize_t size, tsize;
  char *ptr, *buf;
  int fd;

  SB_SB_IFSET(BDMPI_SB_LAZYWRITE) {
#if 1
    /* notify master for each chunk */
    if (SBCHUNK_NONE == (sbchunk->flags&SBCHUNK_NONE)) {
      bdmsg_t msg, gomsg;

      /* notify the master that you want to load memory */
      memset(&msg, 0, sizeof(bdmsg_t));
      msg.msgtype = BDMPI_MSGTYPE_MEMLOAD;
      msg.source  = sbinfo->job->rank;
      msg.count   = sbchunk->npages*sbinfo->pagesize;
      bdmq_send(sbinfo->job->reqMQ, &msg, sizeof(bdmsg_t));
      BDMPL_SLEEP(sbinfo->job, gomsg);
    }
#else
    /* notify master for each page */
    if (SBCHUNK_NONE == (sbchunk->pflags[ip]&SBCHUNK_NONE)) {
      bdmsg_t msg, gomsg;

      /* notify the master that you want to load memory */
      memset(&msg, 0, sizeof(bdmsg_t));
      msg.msgtype = BDMPI_MSGTYPE_MEMLOAD;
      msg.source  = sbinfo->job->rank;
      msg.count   = sbinfo->pagesize;
      bdmq_send(sbinfo->job->reqMQ, &msg, sizeof(bdmsg_t));
      BDMPL_SLEEP(sbinfo->job, gomsg);
    }
#endif
  }

  /* TODO: maybe instead of having to open and close the file each page load,
           the memory could just be file-back mmaped and read from the mmap
           into the memory buffer each page load. */

  ptr = (char *)(sbchunk->saddr + ip*sbinfo->pagesize);

  if (sbchunk->pflags[ip]&SBCHUNK_ONDISK) {
    GKASSERT(0 == (sbchunk->pflags[ip]&SBCHUNK_READ));

    /* make page writeable for load */
    /* have reason to believe that for big enough files and small enough sb
     * parameters `-sb=' and `-pg=', mprotect'ing each page causes the kernel
     * to return a resource shortage error on the call to mprotect. */
    MPROTECT(ptr, sbinfo->pagesize, PROT_WRITE);

    /* load page data from file */
    if (-1 == (fd=open(sbchunk->fname, O_RDONLY))) {
      perror("_sb_pageload: failed to open file");
      exit(EXIT_FAILURE);
    }

    /* need to lseek */
    if (-1 == lseek(fd, ip*sbinfo->pagesize, SEEK_SET)) {
      perror("_sb_pageload: failed on lseek");
      exit(EXIT_FAILURE);
    }

    tsize = sbinfo->pagesize;
    buf = ptr;
    do {
      if (-1 == (size=libc_read(fd, buf, tsize))) {
        perror("_sb_pageload: failed to read the required data");
        exit(EXIT_FAILURE);
      }
      buf   += size;
      tsize -= size;
    } while (tsize > 0);

    if (-1 == close(fd)) {
      perror("sb_pageload: failed to close the fd");
      exit(EXIT_FAILURE);
    }
  }

  /* make page readable */
  MPROTECT(ptr, sbinfo->pagesize, PROT_READ);

  /* update flags to reflect readable */
  if (sbchunk->flags&SBCHUNK_NONE)
    sbchunk->flags ^= SBCHUNK_NONE;
  sbchunk->flags |= SBCHUNK_READ;

  if (sbchunk->pflags[ip]&SBCHUNK_NONE)
    sbchunk->pflags[ip] ^= SBCHUNK_NONE;
  sbchunk->pflags[ip] |= SBCHUNK_READ;

  GKASSERT(SBCHUNK_WRITE != (sbchunk->pflags[ip]&SBCHUNK_WRITE));

  sbchunk->ldpages++;
  GKASSERT(sbchunk->ldpages <= sbchunk->npages);
}


/*************************************************************************/
/*! Saves the supplied ip from sbchunk to disk                           */
/*************************************************************************/
void _sb_pagesave(sbchunk_t * const sbchunk, size_t const ip)
{
  size_t size, tsize;
  char * ptr, * buf;
  int fd;

  ptr = (char *)(sbchunk->saddr + ip*sbinfo->pagesize);

  if (sbchunk->pflags[ip]&SBCHUNK_WRITE) { /* some data have been modified */
    /* open the file for writing, and create it if it does not exist */
    if (-1 == (fd=open(sbchunk->fname, O_RDWR|O_CREAT, S_IRUSR|S_IWUSR))) {
      perror("_sb_pagesave: failed to open file for sbchunk");
      exit(EXIT_FAILURE);
    }

    sbchunk->pflags[ip] |= SBCHUNK_ONDISK;

    if (-1 == lseek(fd, ip*sbinfo->pagesize, SEEK_SET)) {
      perror("_sb_pagesave: failed on lseek");
      exit(EXIT_FAILURE);
    }

    /* write the data */
    tsize = sbinfo->pagesize;
    buf = ptr;
    do {
      if (-1 == (size=libc_write(fd, buf, tsize))) {
        perror("_sb_chunksave: failed to write the required data");
        exit(EXIT_FAILURE);
      }
      buf   += size;
      tsize -= size;
    } while (tsize > 0);

    if (-1 == close(fd)) {
      perror("_sb_pagesave: failed to close the fd data for write");
      exit(EXIT_FAILURE);
    }
  }

  sbchunk->pflags[ip] |= SBCHUNK_NONE;
  if (sbchunk->pflags[ip]&SBCHUNK_READ)
    sbchunk->pflags[ip] ^= SBCHUNK_READ;
  if (sbchunk->pflags[ip]&SBCHUNK_WRITE)
    sbchunk->pflags[ip] ^= SBCHUNK_WRITE;

  /* set madvise */
  if (-1 == madvise(ptr, sbinfo->pagesize, MADV_DONTNEED)) {
    perror("_sb_pagesave: failed to MADV_DONTNEED");
    exit(EXIT_FAILURE);
  }

  /* change protection to PROT_NONE for next time */
  MPROTECT(ptr, sbinfo->pagesize, PROT_NONE);
}


/*************************************************************************/
/*! Frees an sbchunk and its anonymous mmap */
/*************************************************************************/
void _sb_chunkfree(sbchunk_t *sbchunk)
{
  uint8_t flags=sbchunk->flags;
  size_t npages=sbchunk->npages;

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

  /*----------------------------------------------------------------------*/
  SB_SB_IFSET(BDMPI_SB_LAZYWRITE) {
    if (0 == (flags&SBCHUNK_NONE)) {
      bdmsg_t msg, gomsg;

      /* notify the master that you want to save memory */
      memset(&msg, 0, sizeof(bdmsg_t));
      msg.msgtype = BDMPI_MSGTYPE_MEMSAVE;
      msg.source  = sbinfo->job->rank;
      msg.count   = npages*sbinfo->pagesize;
      bdmq_send(sbinfo->job->reqMQ, &msg, sizeof(bdmsg_t));
      BDMPL_SLEEP(sbinfo->job, gomsg);
      //printf("[%d, %d] (%p: -%zu)\n", getpid(), __LINE__, (void*)sbchunk,
      //  msg.count);
      //fflush(stdout);
    }
  }
  /*----------------------------------------------------------------------*/
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
