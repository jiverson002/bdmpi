#include "hll_counter.h"
#include <assert.h>
#include "check_ptr.h"
#include "debug.h"

unsigned long long int HLL_SEED = 0x9fba92b0fa0ce343L;

/**
 * Function to compute the hash function from node IDs.
 *
 * Taken from the WebGraph framework, specifically the class
 * IntHyperLogLogCounterArray.
 *
 * Note that the `x` parameter is a `Long`, but the function will also work
 * with `Int` values.
 *
 * @param x the element to hash, i.e. the node ID
 * @param seed the seed to set up internal state.
 * @return the hashed value of `x`
 */
unsigned long long int jenkins(unsigned long long int x, unsigned long long int seed ) {
  /* Set up the internal state */
  unsigned long long int a = seed + x;
  unsigned long long int b = seed;
  unsigned long long int c = 0x9e3779b97f4a7c13L; /* the golden ratio; an arbitrary value */

  // on unsigned values, the right shift is the logical shift
  a -= b; a -= c; a ^= c >> 43;
  b -= c; b -= a; b ^= a << 9;
  c -= a; c -= b; c ^= b >> 8;
  a -= b; a -= c; a ^= c >> 38;
  b -= c; b -= a; b ^= a << 23;
  c -= a; c -= b; c ^= b >> 5;
  a -= b; a -= c; a ^= c >> 35;
  b -= c; b -= a; b ^= a << 49;
  c -= a; c -= b; c ^= b >> 11;
  a -= b; a -= c; a ^= c >> 12;
  b -= c; b -= a; b ^= a << 18;
  c -= a; c -= b; c ^= b >> 22;

  return c;
}

inline void hll_cnt_init (hll_counter_t *counter, size_t bits) {
  assert( bits < sizeof(hll_hash_t)*8
         && "You cannot use more bits than there are in the hash value");
  counter->b = bits;
  counter->m = 1 << bits; // 2^bits
  if (bits != 0) {
    counter->mask = ~(~0 << (sizeof(hll_hash_t)*8 - bits));
  } else {
    counter->mask = ~0;
  }

  size_t mem = counter->m*sizeof(hll_reg_t);

  counter->registers = (hll_reg_t*) malloc(mem);
  check_ptr(counter->registers);
  memset(counter->registers, 0, mem);
}

inline void hll_cnt_free (hll_counter_t * counter) {
  assert(counter->registers != NULL);
  free(counter->registers);
}

inline hll_counter_t * hll_cnt_new(size_t bits) {
  hll_counter_t *c = malloc(sizeof(hll_counter_t));
  check_ptr(c);
  hll_cnt_init(c, bits);
  return c;
}

inline void hll_cnt_delete(hll_counter_t * counter) {
  hll_cnt_free(counter);
  free(counter);
}

inline hll_counter_t * hll_cnt_copy(hll_counter_t * counter) {
  assert(counter->registers != NULL);
  hll_counter_t *copy = malloc(sizeof(hll_counter_t));
  check_ptr(copy);
  memcpy(copy, counter, sizeof(hll_counter_t));
  copy->registers = malloc(counter->m*sizeof(hll_reg_t));
  check_ptr(copy->registers);
  hll_cnt_copy_to(counter, copy);
  return copy;
}

inline void hll_cnt_copy_to (hll_counter_t *from, hll_counter_t *to) {
  assert(from->b == to->b);
  assert(from->registers != NULL);
  assert(to->registers != NULL);
  memcpy(to->registers, from->registers, from->m*sizeof(hll_reg_t));
}

inline hll_reg_t hll_cnt_rho(hll_hash_t x, hll_hash_t mask) {
  if( (x & mask) == 0 ) {
    return sizeof(hll_hash_t)*8;
  }
  hll_reg_t i = 1;
  hll_hash_t k;
  for(k = x & mask; (k & 1) != 1 ; k = k >> 1) {
    ++i;
  }
  return i;
}

inline void hll_cnt_add(hll_hash_t elem, hll_counter_t * cnt) {
  assert(cnt->registers != NULL);
  hll_hash_t x = jenkins(elem, HLL_SEED);
  // Computes the register index
  size_t j = (cnt->b == 0)? 0 : x >> (sizeof(hll_hash_t)*8 - cnt->b);
  hll_reg_t r = hll_cnt_rho(x, cnt->mask);
  // assigns the maximum between the current value and rho to the register
  cnt->registers[j] = (cnt->registers[j] > r)? cnt->registers[j] : r;
}

inline hll_cardinality_t hll_cnt_size(hll_counter_t * counter) {
  assert(counter->registers != NULL);
  double denominator = 0;
  int i;
  for (i = 0; i<counter->m; ++i) {
    // this `unsigned long long` 1literal is to avoid overflows. Since
    // tipically registers are uint8_t values,
    // using a literal `1` would result in a shift
    // of 8 position of an 8 bit number, resulting in a `0` operand at
    // the denominator.
    denominator += 1.0 / (1ull << counter->registers[i]); // 1 / 2^register
  }
  return counter->m*counter->m*HLL_ALPHA/denominator;
}

#define MAX(a, b) (a >= b)? a : b

inline hll_counter_t * hll_cnt_union(hll_counter_t * a, hll_counter_t * b) {
  assert(a->b == b->b);
  assert(a->registers != NULL);
  assert(b->registers != NULL);
  hll_counter_t * c = hll_cnt_new(a->b);
  int i = 0;
  for(; i<a->m; ++i) {
    c->registers[i] = MAX(a->registers[i], b->registers[i]);
  }
  return c;
}

inline void hll_cnt_union_i(hll_counter_t *a, hll_counter_t *b) {
  assert(a->b == b->b);
  assert(a->registers != NULL);
  assert(b->registers != NULL);
  int i=0;
  for(; i<a->m; ++i) {
    a->registers[i] = MAX(a->registers[i], b->registers[i]);
  }
}

#undef MAX

inline int hll_cnt_equals(hll_counter_t * a, hll_counter_t * b) {
  assert(a->registers != NULL);
  assert(b->registers != NULL);
  if(a->b != b->b) {
    return 0;
  }
  if(memcmp(a->registers, b->registers, a->m*sizeof(hll_reg_t)) != 0) {
    return 0;
  }
  return 1;
}

inline void hll_cnt_replace_registers (hll_counter_t *cnt, hll_reg_t *reg) {
  free(cnt->registers);
  cnt->registers = reg;
}
