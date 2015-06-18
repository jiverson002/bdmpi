#ifndef _GNU_SOURCE
# define _GNU_SOURCE
#endif


#ifdef NDEBUG
# undef NDEBUG
#endif


#include <fcntl.h>     /* O_RDWR, O_CREAT, O_EXCL */
#include <stdint.h>    /* uint8_t, uintptr_t */
#include <stddef.h>    /* NULL, size_t */
#include <stdio.h>     /* FILENAME_MAX */
#include <string.h>    /* memcpy */
#include <sys/mman.h>  /* mmap, mremap, munmap, madvise, mprotect */
#include <sys/stat.h>  /* S_IRUSR, S_IWUSR */
#include <sys/types.h> /* truncate */
#include <unistd.h>    /* truncate */
#include "config.h"
#include "ipc.h"
#include "mmu.h"
#include "vmm.h"


/****************************************************************************/
/*! mtouch function prototype. */
/****************************************************************************/
extern ssize_t __ooc_mtouch__(void * const __addr, size_t const __len);


/****************************************************************************/
/*! Allocate memory via anonymous mmap. */
/****************************************************************************/
extern void *
__ooc_malloc__(size_t const __size)
{
  int ret, fd;
  size_t page_size, s_pages, n_pages, f_pages;
  uintptr_t addr;
  struct ate * ate;
  char fname[FILENAME_MAX];

  /* shortcut */
  if (0 == __size)
    return NULL;

  assert(__size > 0);

  /* compute allocation sizes */
  page_size = vmm.page_size;
  s_pages   = 1+((sizeof(struct ate)-1)/page_size);
  n_pages   = 1+((__size-1)/page_size);
  f_pages   = 1+((n_pages*sizeof(uint8_t)-1)/page_size);

  /* check memory file to see if there is enough free memory to complete this
   * allocation. */
  if (VMM_LZYWR == (vmm.opts&VMM_LZYWR)) {
    ret = __ipc_madmit__(&(vmm.ipc), __vmm_to_sys__(s_pages+n_pages+f_pages));
    if (-1 == ret)
      return NULL;
  }

  /* allocate memory */
  addr = (uintptr_t)mmap(NULL, (s_pages+n_pages+f_pages)*page_size,
    PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS|MAP_NORESERVE|MAP_LOCKED,\
    -1, 0);
  if ((uintptr_t)MAP_FAILED == addr)
    return NULL;

  /* read-only protect application pages -- this will avoid the double SIGSEGV
   * for new allocations */
  ret = mprotect((void*)(addr+(s_pages*page_size)), n_pages*page_size,\
    PROT_READ);
  if (-1 == ret)
    return NULL;

  /* create and truncate the file to size */
  if (0 > snprintf(fname, FILENAME_MAX, "%s%d-%zx", vmm.fstem, (int)getpid(),\
    addr))
  {
    return NULL;
  }
  if (-1 == (fd=libc_open(fname, O_WRONLY|O_CREAT|O_EXCL, S_IRUSR|S_IWUSR)))
    return NULL;
  /*if (-1 == ftruncate(fd, n_pages*page_size))
    return NULL;*/
  if (-1 == close(fd))
    return NULL;

  /* set pointer for the allocation table entry */
  ate = (struct ate*)addr;

  /* populate ate structure */
  ate->n_pages = n_pages;
  ate->l_pages = n_pages;
  ate->base    = addr+(s_pages*page_size);
  ate->flags   = (uint8_t*)(addr+((s_pages+n_pages)*page_size));

  /* initialize ate lock */
  ret = LOCK_INIT(&(ate->lock));
  if (-1 == ret)
    return NULL;

  /* insert ate into mmu */
  ret = __mmu_insert_ate__(&(vmm.mmu), ate);
  if (-1 == ret)
    return NULL;

  /* track number of syspages currently loaded, currently allocated, and high
   * water mark number of syspages */
  __vmm_track__(curpages, __vmm_to_sys__(s_pages+n_pages+f_pages));
  __vmm_track__(numpages, __vmm_to_sys__(s_pages+n_pages+f_pages));
  __vmm_track__(maxpages, vmm.curpages>vmm.maxpages?vmm.curpages-vmm.maxpages:0);

  return (void*)ate->base;
}


/****************************************************************************/
/*! Function for API consistency, adds no additional functionality. */
/****************************************************************************/
extern void *
__ooc_calloc__(size_t const __num, size_t const __size)
{
  return __ooc_malloc__(__num*__size);
}


/****************************************************************************/
/*! Frees the memory associated with an anonymous mmap. */
/****************************************************************************/
extern int
__ooc_free__(void * const __ptr)
{
  int ret;
  size_t page_size, s_pages, n_pages, f_pages, l_pages;
  struct ate * ate;
  char fname[FILENAME_MAX];

  page_size = vmm.page_size;
  s_pages   = 1+((sizeof(struct ate)-1)/page_size);
  ate       = (struct ate*)((uintptr_t)__ptr-(s_pages*page_size));
  n_pages   = ate->n_pages;
  f_pages   = 1+((n_pages*sizeof(uint8_t)-1)/page_size);
  l_pages   = ate->l_pages;

  /* remove the file */
  if (0 > snprintf(fname, FILENAME_MAX, "%s%d-%zx", vmm.fstem, (int)getpid(),\
    (uintptr_t)ate))
  {
    return -1;
  }
  ret = unlink(fname);
  if (-1 == ret && ENOENT != errno)
    return -1;

  /* invalidate ate */
  ret = __mmu_invalidate_ate__(&(vmm.mmu), ate);
  if (-1 == ret)
    return -1;

  /* destory ate lock */
  ret = LOCK_FREE(&(ate->lock));
  if (-1 == ret)
    return -1;

  /* free resources */
  ret = munmap((void*)ate, (s_pages+n_pages+f_pages)*page_size);
  if (-1 == ret)
    return -1;

  /* update memory file */
  if (VMM_LZYWR == (vmm.opts&VMM_LZYWR)) {
    ret = __ipc_mevict__(&(vmm.ipc),\
      -__vmm_to_sys__(s_pages+l_pages+f_pages));
    if (-1 == ret)
      return -1;
  }

  /* track number of syspages currently loaded and allocated */
  __vmm_track__(curpages, -__vmm_to_sys__(s_pages+l_pages+f_pages));
  __vmm_track__(numpages, -__vmm_to_sys__(s_pages+n_pages+f_pages));

  return 0;
}


/****************************************************************************/
/*! Re-allocate memory via anonymous mmap. */
/****************************************************************************/
extern void *
__ooc_realloc__(void * const __ptr, size_t const __size)
{
  return NULL;
}


/****************************************************************************/
/*! Remap an address range to a new address. */
/****************************************************************************/
extern int
__ooc_remap__(void * const __nbase, void * const __obase, size_t const __size,
              size_t const __off)
{
  int ret;
  size_t i, page_size, s_pages, beg, end;
  uintptr_t oaddr, naddr;
  void * optr, * nptr;
  uint8_t * oflags, * nflags;
  struct ate * oate, * nate;
  char ofname[FILENAME_MAX], nfname[FILENAME_MAX];

  page_size = vmm.page_size;
  s_pages   = 1+((sizeof(struct ate)-1)/page_size);
  oate      = (struct ate*)((uintptr_t)__obase-(s_pages*page_size));
  oaddr     = (uintptr_t)oate;
  oflags    = oate->flags;
  optr      = (void*)((uintptr_t)__obase+__off);
  nate      = (struct ate*)((uintptr_t)__nbase-(s_pages*page_size));
  naddr     = (uintptr_t)nate;
  nflags    = nate->flags;
  nptr      = (void*)((uintptr_t)__nbase+__off);

  assert((uintptr_t)__obase == oate->base);
  assert((uintptr_t)__nbase == nate->base);
  assert(nate->l_pages == nate->n_pages);

  /* need to make sure that all bytes are captured, thus beg is a floor
   * operation and end is a ceil operation. */
  beg = ((uintptr_t)nptr-nate->base)/page_size;
  end = 1+(((uintptr_t)nptr+__size-nate->base-1)/page_size);

  /* load old memory */
  ret = __ooc_mtouch__(optr, __size);
  if (-1 == ret)
    return -1;

  /* grant read-write permission to new memory */
  ret = mprotect((void*)((uintptr_t)nptr-__off), __size+__off,\
    PROT_READ|PROT_WRITE);
  if (-1 == ret)
    return -1;

  /* copy memory */
  memcpy(nptr, optr, __size);

  /* grant read-only permission to new memory */
  ret = mprotect((void*)((uintptr_t)nptr-__off), __size+__off, PROT_READ);
  if (-1 == ret)
    return -1;

  for (i=beg; i<end; ++i) {
    assert(MMU_RSDNT != (nflags[i]&MMU_RSDNT));

    /* copy zfill and dirty bit from old flag for clean pages */
    if (MMU_DIRTY != (nflags[i]&MMU_DIRTY))
      nflags[i] |= (oflags[i]&(MMU_ZFILL|MMU_DIRTY));

    /* grant read-write permission to dirty pages of new memory */
    if (MMU_DIRTY == (nflags[i]&MMU_DIRTY)) {
      ret = mprotect((void*)((uintptr_t)__nbase+(i*page_size)), page_size,\
        PROT_READ|PROT_WRITE);
      if (-1 == ret)
        return -1;
    }
  }

  /* move old file to new file and truncate to size */
  ret = snprintf(nfname, FILENAME_MAX, "%s%d-%zx", vmm.fstem, (int)getpid(),\
    naddr);
  if (ret < 0)
    return -1;
  ret = snprintf(ofname, FILENAME_MAX, "%s%d-%zx", vmm.fstem, (int)getpid(),\
    oaddr);
  if (ret < 0)
    return -1;
  ret = rename(ofname, nfname);
  if (-1 == ret)
    return -1;
  /*ret = truncate(nfname, nn_pages*page_size);
  if (-1 == ret)
    return -1;*/

  return 0;
}
