#ifndef _HYPER_ANF_SERIAL_H_
#define _HYPER_ANF_SERIAL_H_

#include "graph.h"

/**
 * @brief Computes the diameter of the given graph.
 *
 * @param graph
 * @param n
 * @param alpha
 * @param bits
 * @param max_iteration
 * @return
 */
int diameter( node_t *graph,
              size_t n,
              int bits,
              int max_iteration);

#endif // _HYPER_ANF_SERIAL_H_
