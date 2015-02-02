#include "hll_counter.h"
#include <check.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <time.h>

START_TEST (hll_counter_rho_test) {
  hll_counter_t *c = hll_cnt_new(3);

  ck_assert_int_eq(1, hll_cnt_rho(0b1010001, c->mask));

  ck_assert_int_eq(5, hll_cnt_rho(0b010000, c->mask));

  ck_assert_int_eq(5, hll_cnt_rho(0b010000, c->mask));

  hll_cnt_delete(c);
}
END_TEST

START_TEST (hll_test_size) {
  /*
   * We will set the register value at 5.
   * The value that the hll_size function should return is
   *
   *     0.72134 / (2^-5) = 23.08 -> 23
   *
   */
  hll_counter_t *c = hll_cnt_new(0);
  c->registers[0] = 5;
  hll_cardinality_t actual = hll_cnt_size(c);
  ck_assert_int_eq(23, actual);

  hll_cnt_delete(c);
}
END_TEST

START_TEST (hll_test_size_2) {
  /*
   * We are going to set the first register at 5 and the second to 3.
   * The value that the hll_size function should return is
   *
   *     0.72134*4 / ( (2^-5) + (2^-3) ) = 18.46 -> 18
   *
   */
  hll_counter_t *c = hll_cnt_new(1);
  c->registers[0] = 5;
  c->registers[1] = 3;
  hll_cardinality_t actual = hll_cnt_size(c);
  ck_assert_int_eq(18, actual);

  hll_cnt_delete(c);
}
END_TEST

int compare (const void * a, const void * b)
{
  return ( *(hll_hash_t*)a - *(hll_hash_t*)b );
}

START_TEST (hll_size_precision) {
#if defined( HASH_SMALL )
  size_t bits = 4;
  int numInsertions = 500;
#else
  size_t bits = 9;
  int numInsertions = 999999;
#endif
  hll_counter_t *c = hll_cnt_new(bits);

  hll_hash_t stream[numInsertions];
  int i=0;
  // Populate stream
  for(; i<numInsertions; ++i) {
    // If we don't do like this, the test fails: the problem is not in the
    // counter but in the random number generator.
    stream[i] = rand() + rand() + rand() + rand() + rand() + rand();
  }
  // Count number of different elements while pushing into counter
  qsort(stream, numInsertions, sizeof(hll_hash_t), compare);
  int n=0;
  int current = stream[0];
  for(i=1; i<numInsertions; ++i) {
    hll_cnt_add(stream[i], c);
    if(stream[i] != current) {
      ++n;
      current = stream[i];
    }
  }

  int estimate = hll_cnt_size(c);
  double error = ((double) estimate / n) * 100 - 100;

  printf("n: %d\nEstimate: %d\nError: %f\n", n, estimate, error);

  ck_assert_int_ne(0, estimate);
  for(i = 0; i < c->m; ++i) {
    ck_assert_int_ne(0, c->registers[i]);
  }
  ck_assert(-10 <= error && error <= 10);

  hll_cnt_delete(c);
}
END_TEST

START_TEST (hll_equals) {
  hll_counter_t
      *c1 = hll_cnt_new(2),
      *c2 = hll_cnt_new(2);

  ck_assert(hll_cnt_equals(c1,c2));

  hll_cnt_add(123, c1);
  hll_cnt_add(123, c2);

  ck_assert(hll_cnt_equals(c1,c2));

  c2->registers[0] = 14; // set the first register to  different value
  ck_assert(!hll_cnt_equals(c1,c2));

  hll_cnt_delete(c1);
  hll_cnt_delete(c2);
}
END_TEST

START_TEST (hll_copy_1) {
  hll_counter_t *c1 = hll_cnt_new(2);
  hll_counter_t *c2 = hll_cnt_copy(c1);

  ck_assert(hll_cnt_equals(c1,c2));

  hll_cnt_delete(c1);
  hll_cnt_delete(c2);
}
END_TEST

START_TEST (hll_copy_2) {
  hll_counter_t *c1 = hll_cnt_new(2);

  c1->registers[0] = 123;
  c1->registers[1] = 43;
  c1->registers[2] = 32;
  c1->registers[3] = 21;

  hll_counter_t *c2 = hll_cnt_copy(c1);

  ck_assert(hll_cnt_equals(c1,c2));

  hll_cnt_delete(c1);
  hll_cnt_delete(c2);
}
END_TEST

START_TEST (hll_copy_to_1) {
  hll_counter_t c1, c2;

  hll_cnt_init(&c1, 2);
  hll_cnt_init(&c2, 2);

  hll_cnt_add(123, &c1);

  hll_cnt_copy_to(&c1, &c2);

  ck_assert(hll_cnt_equals(&c1,&c2));

  hll_cnt_free(&c1);
  hll_cnt_free(&c2);
}
END_TEST

START_TEST (hll_union_i) {
  hll_counter_t *first, *second;

  // first test
  // ==========
  //
  // two counters with 4 registers
  first = hll_cnt_new(2);
  second = hll_cnt_new(2);

  first->registers[0] = 2;
  second->registers[3] = 4;

  hll_cnt_union_i(first, second);

  ck_assert_int_eq(first->registers[0], 2);
  ck_assert_int_eq(first->registers[1], 0);
  ck_assert_int_eq(first->registers[2], 0);
  ck_assert_int_eq(first->registers[3], 4);

  hll_cnt_delete(first);
  hll_cnt_delete(second);

  // second test
  // ===========
  //
  // two counters with 4 registers
  first = hll_cnt_new(2);
  second = hll_cnt_new(2);

  first->registers[0] = 2;
  first->registers[1] = 4;
  first->registers[3] = 6;
  second->registers[1] = 5;
  second->registers[3] = 4;

  hll_cnt_union_i(first, second);

  ck_assert_int_eq(first->registers[0], 2);
  ck_assert_int_eq(first->registers[1], 5);
  ck_assert_int_eq(first->registers[2], 0);
  ck_assert_int_eq(first->registers[3], 6);

  hll_cnt_delete(first);
  hll_cnt_delete(second);
}
END_TEST

START_TEST (hll_union_monotonic) {
  // Each element in the union of two counters should be greater or equals than
  // the corresponding element in the original counters.
  hll_counter_t *first, *second, *u;

  // first test
  // ==========
  //
  // two counters with 4 registers
  first = hll_cnt_new(2);
  second = hll_cnt_new(2);

  first->registers[0] = 2;
  second->registers[3] = 4;

  u = hll_cnt_union(first, second);

  for(int i=0; i<4; ++i) {
    ck_assert(u->registers[i] >= first->registers[i]);
    ck_assert(u->registers[i] >= second->registers[i]);
  }

  hll_cnt_delete(first);
  hll_cnt_delete(second);
  hll_cnt_delete(u);

  // second test
  // ===========
  //
  // two counters with 4 registers
  first = hll_cnt_new(2);
  second = hll_cnt_new(2);

  first->registers[0] = 2;
  first->registers[1] = 4;
  first->registers[3] = 6;
  second->registers[1] = 5;
  second->registers[3] = 4;

  u = hll_cnt_union(first, second);

  for(int i=0; i<4; ++i) {
    ck_assert(u->registers[i] >= first->registers[i]);
    ck_assert(u->registers[i] >= second->registers[i]);
  }

  hll_cnt_delete(first);
  hll_cnt_delete(second);
  hll_cnt_delete(u);
}
END_TEST

START_TEST (hll_union_i_monotonic) {
  // Each element in the union of two counters should be greater or equals than
  // the corresponding element in the original counters.
  hll_counter_t *first, *second, *u;

  // first test
  // ==========
  //
  // two counters with 4 registers
  first = hll_cnt_new(2);
  second = hll_cnt_new(2);

  first->registers[0] = 2;
  second->registers[3] = 4;
  u = hll_cnt_copy(first);

  hll_cnt_union_i(u, second);

  for(int i=0; i<4; ++i) {
    ck_assert(u->registers[i] >= first->registers[i]);
    ck_assert(u->registers[i] >= second->registers[i]);
  }

  hll_cnt_delete(first);
  hll_cnt_delete(second);
  hll_cnt_delete(u);

  // second test
  // ===========
  //
  // two counters with 4 registers
  first = hll_cnt_new(2);
  second = hll_cnt_new(2);

  first->registers[0] = 2;
  first->registers[1] = 4;
  first->registers[3] = 6;
  second->registers[1] = 5;
  second->registers[3] = 4;

  u = hll_cnt_copy(first);

  hll_cnt_union_i(u, second);

  for(int i=0; i<4; ++i) {
    ck_assert(u->registers[i] >= first->registers[i]);
    ck_assert(u->registers[i] >= second->registers[i]);
  }

  hll_cnt_delete(first);
  hll_cnt_delete(second);
  hll_cnt_delete(u);
}
END_TEST

#define RANDOM rand() + rand() + rand() + rand()
START_TEST(hll_size_monotonic) {
  // the estimated size of the union of two counters should be greater or
  // equals than the estimated size of each original counter.
  hll_counter_t a, b;

  hll_cnt_init(&a, 2);
  hll_cnt_init(&b, 2);

  for(int i=0; i<1024; ++i) {
    hll_cnt_add(RANDOM, &a);
    hll_cnt_add(RANDOM, &b);
  }

  hll_cardinality_t
      size_a = hll_cnt_size(&a),
      size_b = hll_cnt_size(&b);

  hll_cnt_union_i(&a, &b);

  hll_cardinality_t size_union = hll_cnt_size(&a);

  printf("size_a: %u\n"
         "size_b: %u\n"
         "size_union: %u\n",
         size_a, size_b, size_union);

  ck_assert(size_union >= size_a);
  ck_assert(size_union >= size_b);

  hll_cnt_free(&a);
  hll_cnt_free(&b);
}
END_TEST
#undef RANDOM

START_TEST(hll_size_monotonic_2) {
  // this test case arised during debugging: using the number of bits of the
  // hash type for the `0` element will cause the whole program to fail, since
  // the neighbourhood functions will be decreasing.
  hll_counter_t a, b;

  hll_cnt_init(&a, 1);
  hll_cnt_init(&b, 1);

  a.registers[0] = 32;
  a.registers[1] = 0;

  b.registers[0] = 2;
  b.registers[1] = 0;

  hll_cardinality_t
      size_a = hll_cnt_size(&a),
      size_b = hll_cnt_size(&b);

  printf("size_a: %u\n"
         "size_b: %u\n",
         size_a, size_b);

  ck_assert(size_a >= size_b);
}
END_TEST

Suite * hll_counter_suite () {
  Suite *s = suite_create("HyperLogLog Counter");
  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, hll_equals);
  tcase_add_test(tc_core, hll_union_i);
  tcase_add_test(tc_core, hll_union_monotonic);
  tcase_add_test(tc_core, hll_union_i_monotonic);
  tcase_add_test(tc_core, hll_copy_1);
  tcase_add_test(tc_core, hll_copy_2);
  tcase_add_test(tc_core, hll_copy_to_1);
  tcase_add_test(tc_core, hll_counter_rho_test);
  tcase_add_test(tc_core, hll_test_size);
  tcase_add_test(tc_core, hll_test_size_2);
  tcase_add_test(tc_core, hll_size_precision);
  tcase_add_test(tc_core, hll_size_monotonic);
  tcase_add_test(tc_core, hll_size_monotonic_2);
  suite_add_tcase(s, tc_core);

  return s;
}

int main() {
  printf(
        "Platform information:\n"
        "sizeof(int): %lu\n"
        "sizeof(unsigned int): %lu\n"
        "sizeof(hll_reg_t): %lu\n"
        "sizeof(hll_hash_t): %lu\n",
        sizeof(int),
        sizeof(unsigned int),
        sizeof(hll_reg_t),
        sizeof(hll_hash_t)
        );

  srand(time(NULL));

  int number_failed = 0;
  Suite *s = hll_counter_suite();
  SRunner *sr = srunner_create(s);
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
