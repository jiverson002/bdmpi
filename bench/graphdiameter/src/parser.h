/**
 * @file parser.h
 * @brief parser definitions.
 *
 * This header contains functions that parse files and strings to build graph
 * nodes.
 */

#ifndef _PARSER_H_
#define _PARSER_H_

#include "graph.h"

/**
 * @brief Builds a heap-allocated node from a string description
 *
 * This function takes a string that contains the description of a node, parses
 * it, allocates a new node_t struct and populates it with the correct values.
 *
 * **Attention**: it the description os not well formed, the function
 * returns NULL.
 *
 * Description format
 * ------------------
 *
 * The format used in node descriptions _must_ follow this specification
 *
 *     id | <forward star> | <reverse star>
 *
 * that is: the id of the node, a pipe symbol, a list of identifiers for the
 * outgoing neighbours, a pipe symbol, a list of identifiers for the incomping
 * neighbours.
 *
 * ### Example ###
 *
 * For example we have a node, `0`, whose outgoing neighbours are `1`, `3` and
 * `5` and whose incoming neighbour is `2`. The string to describe this node is
 * the following
 *
 *     0 | 1 3 5 | 2
 *
 * A node with no neighbours is represented in this way (suppose its id is `2`)
 *
 *     2 | |
 *
 * That is, the pipes are **mandatory**
 *
 * @param descr the description of the node
 * @return a heap allocated node struct, NULL if there are errors
 */
node_t * parse_node_descr(char * descr);

/**
 * @brief Parses a node description populating the given node.
 *
 * This is the same of parse_node_descr(char *), but instead of initializing
 * a node_t on the heap, populates an already existing one.s
 *
 * @see parse_node_descr(char *)
 *
 * @param descr the description to parse
 * @param node the node to populate
 * @return 0 on success, 1 otherwise
 */
int parse_node_descr_to(char * descr, node_t * node);

/**
 * @brief Loads an array of node_t structs from the given file.
 *
 * The array that will hold the nodes will be dinamically allocated by the
 * function itself;
 *
 * Example usage
 *
 * ~~~~~{.c}
 * int n;
 * node_t *nodes;
 *
 * int rc = parse_graph_file("graph.dat", &nodes, &n);
 *
 * if (rc < 0) {
 *   // ERROR!
 * } else {
 *   // Everything ok, you can use the nodes array that has been
 *   // correctly allocated
 * }
 * ~~~~~
 *
 * @param filename the name of the file containing the nodes
 * @param nodes the (uninitialized) array of node_t that will hold the
 *              loaded nodes
 * @param n pointer to the int that will hold the size of the array.
 * @return the number of parsed nodes on success, -1 otherwise
 */
int parse_graph_file(char *filename, node_t **nodes, int *n);

/**
 * @brief Parses a string description of a graph, possibly loaded from a file.
 *
 * Used as a subroutine by parse_graph_file.
 *
 * @param str the string description of the graph
 * @param nodes the node_t array to populate
 * @param n pointer to the int that will hold the size of the nodes array
 * @return the number of nodes parsed.
 */
int parse_graph_string(char *str, node_t **nodes, int *n);

/*
 * Given a string, counts the number of distinct numbers in it.
 * Returns `-1` on error
 */
int count_numbers(char * str);

/*
 * Reads the number of lines in the given string.
 */
int count_lines(char * str);

/*
 * Populates the `adj` array with `n` integers extracted from the string
 * `adj_str`.
 * Returns -1 on error, the number of inserted elements otherwise
 */
int populate_adjacency(char * adj_str, node_id_t * adj, int n);

#endif // _PARSER_H_
