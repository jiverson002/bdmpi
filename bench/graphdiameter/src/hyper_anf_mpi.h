#ifndef _HYPER_ANF_MPI_H_
#define _HYPER_ANF_MPI_H_

/**
 * @file hyper_anf_mpi.h
 */

#include <mpi.h>
#include "graph.h"
#include "hll_counter.h"

/**
 * Defines an association between an array of counters and its dimension.
 * Useful to keep track of the neighbours' counters.
 */
struct mpi_neighbourhood {
  size_t dimension;
  hll_counter_t *counters;
};
typedef struct mpi_neighbourhood mpi_neighbourhood_t;

/**
 * Struct that carries all the information needed by the diameter computation.
 */
struct context {
  size_t num_nodes; /**< The number of nodes under the responsibility of
                      *  this processor.
                      */
  node_t *nodes; /**< The nodes under the responsibility of this processor */
  hll_counter_t *counters; /**< The array of current counters, one for each
                             *  node. The size is `num_nodes`
                             */
  hll_counter_t *counters_prev; /**< The array of previous counters, one for
                                  *  each node. The size is `num_nodes`
                                  */
  mpi_neighbourhood_t *neighbourhoods; /**< The array of out neighbourhoods.
                                       *  Size is `num_nodes`
                                       */
  MPI_Request *requests; /**< Array that holds all the receive and send
                           * requests. The number of elements of this array is
                           * $$
                           * \\sum_{v \\in nodes} |\\delta^+(v)| +
                           *                      |\\delta^-(v)|
                           * $$
                           */
  size_t num_requests; /**< The number of requests. See requests.
                         */
  int num_registers; /**< The number of registers used in counters */
  int iteration; /**< The current algorithm iteration */
  int num_changed; /**< The number of nodes changed since the last iteration */
  int num_processors; /**< The number of processors */
  int max_iteration; /**< The maximum iteration number we are allowed to do */
  int rank; /**< The rank of the processor */
  double alpha; /**< We compute the effective diameter at `alpha` */
  double *neighbourhood_function; /**< The neighbourhood function as an
                                       array of length max_iteration */
};
typedef struct context context_t;

void init_context(context_t * context,
                  node_t *partial_graph,
                  size_t partial_graph_cardinality,
                  int bits,
                  int max_iteration,
                  double alpha);

/**
 * @brief free_context
 *
 * **Attention**: does not free the contents of the `nodes` array.
 *
 * @param context
 */
void free_context(context_t * context);

/**
 * @brief Compute the diameter of the graph, given the partial view.
 *
 * **Attention**: assumes that MPI_init has already been called.
 *
 * Usage
 * -----
 *
 * A context must be allocated before using this function:
 *
 * ~~~~~c
 * MPI_Init(argc, argv);
 * // ...
 * context_t ctx;
 * init_context(&ctx, partial_graph, num_nodes, bits, max_iteration);
 * mpi_diameter(&ctx);
 * free_context(&ctx);
 * // ...
 * MPI_Finalize();
 * ~~~~~
 *
 * @param partial_graph
 * @param partial_graph_cardinality
 * @param bits
 * @param max_iteration
 * @return
 */
double mpi_diameter(context_t *context);

double effective_diameter(context_t *context);

void compute_neighbourhood_function(context_t * context);

/**
 * @brief Returns the global rank of the processor responsible for the node.
 * @param node
 * @param num_processors
 * @return
 */
int get_processor_rank( node_id_t node, int num_processors );

void mpi_neighbourhood_init(mpi_neighbourhood_t *neigh, int n, int bits);

void mpi_neighbourhood_free(mpi_neighbourhood_t *neigh);

/**
 * @brief loads the partial graph under the responsibility of this node.
 *
 * @param rank The rank of the processor loading the file
 * @param filename the file name of the partial graph.
 * @param nodes a pointer to a `node_id_t` array that will be
 *              initialized with the ids of the nodes under the
 *              responsibility of the current processor.
 * @param n a pointer to the integer that will define the lenght of
 *          the `nodes` array
 */
void load_partial_graph( int rank,
                         char * basename,
                         node_t **nodes,
                         int *n);

#endif // _HYPER_ANF_MPI_H_
