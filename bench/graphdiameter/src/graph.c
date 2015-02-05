#include "graph.h"
#include "null.h"
#include "hll_counter.h"

static size_t _ne=0;
static node_id_t * edges;

static size_t _nn=0;
static hll_counter_t * neighbourhood;

void init_edges(size_t ne) {
  if (NULL == (edges=malloc(ne*sizeof(node_id_t))))
    abort();
}

void free_edges(void) {
  free(edges);
}

void init_neighbourhood(size_t ne) {
  if (NULL == (neighbourhood=malloc(ne*sizeof(hll_counter_t))))
    abort();
}

void free_neighbourhood(void) {
  free(neighbourhood);
}

void node_init (node_t *node, node_id_t id, size_t out_n, size_t in_n) {
  node->id = id;
  node->num_out = out_n;
  node->num_in = in_n;
  if (out_n == 0) {
    node->out = NULL;
  } else {
    node->out = edges+_ne;
    _ne += out_n;
  }
  if (in_n == 0) {
    node->in = NULL;
  } else {
    node->in = edges+_ne;
    _ne += in_n;
  }
}

void neighbour_init (hll_counter_t ** counters, size_t n) {
  *counters = neighbourhood+_nn;
  _nn += n;
}
