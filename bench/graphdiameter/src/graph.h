#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <inttypes.h>
#include <stdlib.h>

// ----------------------------------------------------------------------------
//    Data structures
// ----------------------------------------------------------------------------

typedef uint32_t node_id_t;

typedef struct _node {
  node_id_t id;     /**< Unique ID of the node */
  size_t num_out;   /**< Number of outgoing neighbours */
  size_t num_in;    /**< Number of outgoing neighbours */
  node_id_t *out;   /**< Array of outgoing neighbours */
  node_id_t *in;    /**< Array of incoming neighbours */
} node_t;


// ----------------------------------------------------------------------------
//    Memory management routines
// ----------------------------------------------------------------------------

/**
 * @brief Dynamically allocates a new node.
 *
 * @see node_init(node_t *, node_id_t, size_t, size_t);
 *
 * @param id the id of the node
 * @param out_n the cardinality of the forward star of the node
 * @param in_n the cardinality of the reverse star
 * @return a pointer to the heap-allocated node.
 */
node_t * node_new(node_id_t id, size_t out_n, size_t in_n);

/**
 * @brief Deletes a heap-allocated node.
 * @param node the node to delete.
 */
void node_delete(node_t * node);

/**
 * @brief Initializes an already allocated node_t struct
 *
 * To be used to initialize already allocated nodes. This may be useful to
 * init nodes that are allcoated in the stack or in an array.
 *
 * Example usage
 *
 * ~~~~{.c}
 * node_t node;
 * node_init (&node, 1, 3, 4);
 * ~~~~
 *
 * **Attention**: if `out_n` or `in_n` are 0, the corresponding array will be
 * set to `NULL`.
 *
 * @param node Pointer to the node to initialize
 * @param id the id of the node
 * @param out_n the out degree of the node
 * @param in_n the in degree of the node
 */
void node_init(node_t * node, node_id_t id, size_t out_n, size_t in_n);

/**
 * @brief Frees the resources associated to a node.
 *
 * To be used with nodes allocated with
 * node_init(node_t * , node_id_t , size_t , size_t )
 *
 * This function will free the arrays associated with the nodes and that live
 * on the heap. It will _not_ free the memory associated with the node itself.
 *
 * @param node the node to free
 */
void node_free(node_t * node);

#endif // _GRAPH_H_
