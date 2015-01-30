#include "klutil.h"

#if !defined(__INTEL_COMPILER) && !defined(__GNUC__)
  static int kl_builtin_popcountl(size_t v) {
    int c = 0;
    for (; v; c++) {
      v &= v-1;
    }
    return c;
  }
  static int kl_builtin_clzl(size_t v) {
    /* result will be nonsense if v is 0 */
    int i;
    for (i=sizeof(size_t)*CHAR_BIT-1; i>=0; --i) {
      if (v & (1ul << i)) {
        break;
      }
    }
    return sizeof(size_t)*CHAR_BIT-i-1;
  }
  static int kl_builtin_ctzl(size_t v) {
    /* result will be nonsense if v is 0 */
    int i;
    for (i=0; i<sizeof(size_t)*CHAR_BIT; ++i) {
      if (v & (1ul << i)) {
        return i;
      }
    }
    return i;
  }
  #define kl_popcount(V)  kl_builtin_popcountl(V)
  #define kl_clz(V)       kl_builtin_clzl(V)
  #define kl_ctz(V)       kl_builtin_ctzl(V)
#else
  #define kl_popcount(V)  __builtin_popcountl(V)
  #define kl_clz(V)       __builtin_clzl(V)
  #define kl_ctz(V)       __builtin_ctzl(V)
#endif

void * sb_malloc(size_t sz);
void sb_free(void * ptr);

size_t
kl_ilog2(
  size_t v
)
{
  return sizeof(size_t)*CHAR_BIT-1-__builtin_clzl(v);
}

size_t
kl_pow2up(
  size_t v
)
{
  /* assumes that sizeof(size_t) == __WORDSIZE/8 */
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
#if defined(__WORDSIZE) && __WORDSIZE == 64
  v |= v >> 32;
#endif
  return v+1;
}

size_t
kl_multup(
  size_t v,
  size_t m
)
{
  if (m == 0) {
    return 0;
  }

  return (v+m-1)/m*m;
}

void *
kl_mmap_aligned(
  size_t const sz
)
{
  return sb_malloc(sz);
}

void
kl_munmap_aligned(
  void * const ptr,
  size_t const sz
)
{
  sb_free(ptr);
}

size_t
kl_vmem()
{
  int                pid, ppid, pgrp, session, tty_nr, tpgid;
  char               state;
  unsigned           flags;
  long int           cutime, cstime, priority, nice, num_threads, itrealvalue;
  long unsigned      minflt, cminflt, majflt, cmajflt, utime, stime, vsize;
  /*long */long unsigned starttime;
  char               comm[1024];

  FILE * f = fopen("/proc/self/stat", "r");
  
  if(fscanf(f, "%d %s %c %d %d %d %d %d %u %lu %lu %lu %lu %lu %lu %ld %ld "
               "%ld %ld %ld %ld %lu %lu", &pid, comm, &state, &ppid, &pgrp,
                &session, &tty_nr, &tpgid, &flags, &minflt, &cminflt, &majflt,
                &cmajflt, &utime, &stime, &cutime, &cstime, &priority, &nice,
                &num_threads, &itrealvalue, &starttime, &vsize) != 23) {}
  fclose(f);

  return vsize;
}
