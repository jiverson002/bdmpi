#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "hyper_anf_serial.h"
#include "hll_counter.h"
#include "check_ptr.h"
#include "debug.h"

int diameter (node_t *graph,
              size_t n,
              int bits,
              int max_iteration)
{
  // Declare counter arrays
  hll_counter_t *counters = malloc(n*sizeof(hll_counter_t));
  check_ptr(counters);
  hll_counter_t *counters_prev = malloc(n*sizeof(hll_counter_t));
  check_ptr(counters_prev);

  // Declare the neighbourhood function array
  hll_cardinality_t *neighbourhood_function =
      malloc(max_iteration*sizeof(hll_cardinality_t));
  check_ptr(neighbourhood_function);
  memset(neighbourhood_function, 0, max_iteration*sizeof(hll_cardinality_t));

  // Initialize counter arrays++
  for(int i=0; i<n; ++i) {
    hll_cnt_init(& (counters[i]), bits);
    hll_cnt_init(& (counters_prev[i]), bits);

    // and add the node itself to each counter
    hll_cnt_add(graph[i].id, &(counters[i]));
    hll_cnt_add(graph[i].id, &(counters_prev[i]));
  }

  int iteration = 0;
  int num_changed = n;

  // Main iteration
  while (num_changed != 0 && iteration < max_iteration) {
    // reset the number of changed nodes
    num_changed = 0;
    // Compute the neighbourhood function before updating the counters.
    hll_cardinality_t sum = 0;
    for(int h=0; h<n; ++h) {
      printd("sum = %u\n", sum);
      sum += hll_cnt_size(& (counters[h]));
    }
    neighbourhood_function[iteration] = sum;
    printf("N(%d) = %d\n", iteration, sum);

    // Update the counters:
    // For each node, take its counter
    for (int i=0; i<n; ++i) {
      printd("node: %d\n", i);
      hll_counter_t node_counter = counters[i];
      // Then, for each neighbour
      for (int j=0; j<graph[i].num_out; ++j) {
        node_id_t neighbour = graph[i].out[j];
        hll_counter_t neighbour_counter = counters_prev[neighbour];
        // perform the union of the two counters
        hll_cnt_union_i(&node_counter, &neighbour_counter);
#ifndef NDEBUG // check that the counter increased
        hll_counter_t prev = counters_prev[i];
        printf("| current |    prev |\n");
        for(int r=0; r<node_counter.m; ++r) {
          printf("| %7u | %7u |\n",
                 node_counter.registers[r], prev.registers[r]);
          assert(node_counter.registers[r] >= prev.registers[r]);
        }
        printf("current size: %u previous size: %u\n",
               hll_cnt_size(&node_counter), hll_cnt_size(&prev));
        assert(hll_cnt_size(&node_counter) >= hll_cnt_size(&prev));
#endif
      }
    }

    // Check the counters that changed and then copy the new value in the
    // `prev` version of the counter
    for (int i=0; i<n; ++i) {
      hll_counter_t
          node_counter = counters[i],
          node_counter_prev = counters_prev[i];
      assert(hll_cnt_size(&node_counter) >= hll_cnt_size(&node_counter_prev));
      // If the counter changed with respect to its previous value, increment
      // num_changed variable
      if (!hll_cnt_equals(&node_counter, &node_counter_prev)) {
        ++num_changed;
      }
      // Then copy from `counters` to `counters prev`
      hll_cnt_copy_to(&node_counter, &node_counter_prev);
    }

    ++iteration;
  }

  // Free resources
  for (int i=0; i<n; ++i) {
    hll_cnt_free(&counters[i]);
    hll_cnt_free(&counters_prev[i]);
  }
  free(counters);
  free(counters_prev);
  free(neighbourhood_function);

  return iteration-1;
}
