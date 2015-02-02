/**
 * @file   hll_counter.h
 * @brief  Header for HyperLogLog counters.
 *
 * This header defines data structures and operations related to
 * HyperLogLog counters
 * [[FlFuGaMe07](http://algo.inria.fr/flajolet/Publications/FlFuGaMe07.pdf)].
 *
 * Introduction
 * ------------
 *
 * This counters are probabilistic data structures used to keep track
 * of the cardinality of large multisets.
 *
 * Operations allowed on the counter
 * ---------------------------------
 *
 * On a counter we can perform the following operations
 *
 *  - add a new element: hll_cnt_add(hll_hash_t, hll_counter_t)
 *  - query for the estimated size of the underlying multiset:
 *    hll_cnt_size(hll_counter_t)
 *  - get the union of two counters: hll_cnt_union(hll_counter_t, hll_counter_t)
 *
 * Algorithm
 * ---------
 *
 * The HyperLogLog enables to keep track of the cardinalities of very large
 * multisets in very little space.
 *
 * In fact this algorithm only needs \f$O(\log\log n)\f$
 * memory (hence the name) to keep track of cardinalities up to \f$n\f$.
 *
 * In the description that follows we assume that we have an hash function
 * \f$h(x)\f$ that maps the set elements into binary strings.
 * The function \f$\rho(y)\f$ is used to get the position of the leftmost 1-bit
 * in the binary string \f$y\f$.
 * 
 * This algorithm, uses stochastic averaging to control
 * the standard error on the estimated value. The input multiset
 * \f$\mathcal{M}\f$ is partitioned in \f$m = 2^b\f$ multisets
 * \f$\mathcal{M}_1 \dots \mathcal{M}_m\f$ using the first \f$b\f$ bits of the
 * hashed values.
 * 
 * Then for each \f$\mathcal{M}_j\f$ the
 * observable is
 * \f[
 *   M^j = \max_{x \in \mathcal{M}_j} \rho(x)
 * \f]
 * 
 * The intuition underlying the algorithm is the following: each multiset
 * \f$\mathcal{M}_j\f$ is expected to contain \f$n/m\f$ distinct elements at
 * the end of the execution. Hence the parameter \f$M^j\f$ should be close to
 * \f$\log_2(n/m)\f$. The harmonic mean of the quantities \f$2^{M^j}\f$ should
 * then be in the order of \f$n/m\f$.
 * 
 * The algorithm returns the following estimate of the cardinality
 * \f[
 *  N = \frac{\alpha m^2}{\sum_{j=1}^{m} 2^{-M^j}}
 * \f]
 * 
 * The parameter \f$\alpha\f$ is used to correct a multiplicative bias of the
 * algorithm and is defined in macro ::HLL_ALPHA.
 *
 * Compilation notes
 * =================
 *
 * When compiling, the dimension of the has values is conditionally determined
 * based on a couple of macros:
 *
 *  - `HASH_SMALL`: uses 8 bit hashes, that is `uint8_t`
 *  - `HASH_MEDIUM`: use `uint16_t` integers for hash values
 *  - `HASH_BIG`: use `uint32_t` integers, this is the default
 *  - `HASH_HUGE`: use 64 bit hashes, that is `uint64_t` integers
 *
 * The rationale behind using different integer types is that graphs that have
 * a number of nodes significantly smaller than the maximum representable by
 * the datatype will underutilize the registers. The problem with this is that
 * in this way we get only the first couple of registers with nonzero values,
 * heavily biasing the cardinality extimation.
 *
 * To select the right flag at compilation time, specify it on the command line
 * of `cmake`:
 *
 *     mkdir gdem_small
 *     cd gdem_small
 *     cmake -DCMAKE_BUILD_TYPE=Release -DHASH_SIZE=Small ..
 *     make
 *
 * Thus we will get different binaries, each suitable for different graph sizes.
 */

#ifndef _HLL_COUNTER_H_
#define _HLL_COUNTER_H_

#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

/**
 * @brief Numeric type of estimated cardinalities.
 *
 * This typedef is used to hide the actual numeric type used to represent
 * multiset cardinalities.
 */
typedef double hll_cardinality_t;

// ----------------------------------------------------------------------------
//   Counter Definitions
// ----------------------------------------------------------------------------

/**
 * Correction constant for the estimation
 */
#define HLL_ALPHA 0.72134

/**
 * @brief Type of counter registers.
 *
 * This data type represents a single register of the counter.
 */
typedef uint8_t hll_reg_t;

/**
 * @brief Type of hash values.
 *
 * This data type is the input of the add operation.
 */
#if defined( HASH_SMALL )
typedef uint8_t hll_hash_t;
#elif defined( HASH_MEDIUM )
typedef uint16_t hll_hash_t;
#elif defined( HASH_HUGE )
typedef uint64_t hll_hash_t;
#else // HASH_BIG is the default
typedef uint32_t hll_hash_t;
#endif

/**
 * @brief This struct is the actual counter.
 *
 * It contains an array of registers and their number, the number of
 * bits to be used in partitioning operations and the mask to be used during
 * update operations.
 */
struct hll_counter {
  /**
   * The number of registers.
   */
  size_t m;
  /**
   * The number of bits used to select the register to update.
   */
  unsigned int b;
  /**
   * The mask used during updates.
   */
  hll_hash_t mask;
  /**
   * The array of registers of this counter.
   */
  hll_reg_t *registers;
};

/**
 * @brief Type of the counter.
 */
typedef struct hll_counter hll_counter_t;

// ----------------------------------------------------------------------------
//   Memory Management for counters
// ----------------------------------------------------------------------------

/**
 * @brief Creates a new hyperLogLog counter, with the specified number of
 * bits to index the registers.
 *
 * The number of registers is m = 2^bits.
 *
 * **Attention**: the memory for the registers of the counter is dynamically
 * allocated. Hence, once used, the counter _must_ be deallocated with the
 * function hll_cnt_delete(hll_counter_t).
 *
 * @param bits: the bits to index the registers.
 * @return a new hll_counter_t instance
 */
hll_counter_t * hll_cnt_new(size_t bits);

/**
 * @brief Deletes a counter.
 *
 * This functions is responsible of deallocating the right amount of memory.
 *
 * @param counter the counter to be deallocated.
 */
void hll_cnt_delete(hll_counter_t * counter);

/**
 * @brief Copies a counter.
 *
 * A deep copy of the counter is performed: the pointer ot the registers array
 * is not shared.
 *
 * @param counter the counter to copy
 * @return the new copy of the counter.
 */
hll_counter_t * hll_cnt_copy(hll_counter_t * counter);

/**
 * @brief Copies the values of the registers of from in the registers of to
 *
 * **Attention**: the `to` counter _must_ be already initialized.
 *
 * @param from the original counter
 * @param to the counter we are copying to
 */
void hll_cnt_copy_to(hll_counter_t * from, hll_counter_t * to);

/**
 * @brief Initializes an already allocated counter.
 * @param counter the counter to initialize
 * @param bits the number of bits to index the registers
 */
void hll_cnt_init(hll_counter_t * counter, size_t bits);

/**
 * @brief Deallocates the memory of the registers of the counter.
 */
void hll_cnt_free(hll_counter_t * counter);

// ----------------------------------------------------------------------------
//   Counter operations
// ----------------------------------------------------------------------------

/**
 * @brief Returns the position of the rightmost 1-bit in a hash.
 *
 * **Attention**: bit indexes returned by this method
 * start at 1, not at 0!
 *
 * @param elem the hash to be evaluated
 * @param mask the bit mask to hide the first `b` bits of the hash.
 * @return the position of the rightmost 1-bit in the hash, starting from 1.
 */
hll_reg_t hll_cnt_rho(hll_hash_t elem, hll_hash_t mask);

/**
 * @brief Add a new element to the counter.
 *
 * @param elem the element to be added.
 * @param counter the counter to be updated.
 */
void hll_cnt_add(hll_hash_t elem, hll_counter_t * counter);

/**
 * @brief Estimates the cardinality of the underlying multiset.
 *
 * @param counter the counter used to get the estimation.
 * @return the estimated cardinality of the multiset.
 */
hll_cardinality_t hll_cnt_size(hll_counter_t * counter);

/**
 * @brief Performs the union of two counters.
 *
 * The union of two counters is defined as the register by register maximum.
 *
 * @param a the first counter
 * @param b the second counter
 * @return the union of the two counters
 */
hll_counter_t * hll_cnt_union(hll_counter_t * a, hll_counter_t * b);

/**
 * @brief Performs the inplace union of two counters.
 *
 * This is the same as hll_cnt_union(hll_counter_t, hll_counter_t) but instead
 * of returning a newly allocated counter the result is placed in the first
 * parameter.
 *
 * @param a the first counter and the one that will hold the result
 * @param b the second counter
 */
void hll_cnt_union_i(hll_counter_t * firstAndResult, hll_counter_t * second);

/**
 * @brief Checks if two counters are equals.
 *
 * Two counters are equal if they have the same `b` and their registers
 * are equals.
 *
 * @param a the first counter
 * @param b the second counter
 * @return `0` if the counters are not equals.
 */
int hll_cnt_equals(hll_counter_t * a, hll_counter_t * b);

/**
 * @brief replaces the array of registers with the given one
 *
 * Frees the memory allocated to the old registers.
 * Note that it is your responsibility to ensure that the new array has the same
 * size fo the old one
 *
 * @param cnt the counter
 * @param reg the array of registers we want to inject
 */
void hll_cnt_replace_registers(hll_counter_t * cnt, hll_reg_t * reg);

#endif // _HLL_COUNTER_H_
