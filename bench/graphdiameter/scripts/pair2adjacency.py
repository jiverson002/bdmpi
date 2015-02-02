#!/usr/bin/env python

#################################################################################
#
# Script that accepts a graph incidence matrix on standard input and outputs the 
# adjacency list representation of the graph, in the format accepted and 
# described in parser.h
# 
# The documentation can be generated from these comments using
# [`pycco`](http://fitzgen.github.io/pycco/).
#
#################################################################################

# Purpose
# =======
# Most graph datasets are stored as plain text files where each line corresponds
# to an edge written as a pair:
#
#     source_id destination_id
#
# For undirected graphs the same edge is not written twice.
# 
# The `gdem` program accepts as input files in the following format
#
#     node_id | outgoing neighbours | incoming neighbours
#
# The purpose of this script is to convert between the two formats.
#
#
# Usage
# =====
# 
# The script reads from standard input the graph in pair format and writes to
# standard output the adjacency list format.
# 
# To convert a graph stored in a file named `graph.pair` issue the following
# command
#
#     cat graph.pair | python pair2adjacency > graph.adj
#
# for undirected graphs. In the case of directed graphs the command is the
# following
#
#     cat graph.pair | python pair2adjacency -d > graph.adj
#


# Implementation
# ==============


import sys
from optparse import OptionParser


# Adjacencies dictionary 
# ----------------------
# The key is the node ID, the value is
# a pair: the first position is a list of the outgoing neighbours, the second is
# a list of the incoming neighbours
adjacencies = dict()

# Adding edges to the dictionary
# ------------------------------
# 
# Adds an edge between `u` and `v`. If the parameter `directed` is set to `True`
# then the edge is directed.
def add_edge(u, v, directed = False):
    # Add `v` to the outgoing nodes of `u`
    if u in adjacencies:
        adjacencies[u][0].append(v)
    else:
        adjacencies[u] = ([v],[])
    # Add `u` to the incoming nodes of `v`
    if v in adjacencies:
        adjacencies[v][1].append(u)
    else:
        adjacencies[v] = ([],[u])
    # If the graph si undirected
    if not directed:
        # Add `u` to the outgoing nodes of `v`
        if v in adjacencies:
            adjacencies[v][0].append(u)
        else:
            adjacencies[v] = ([u],[])
        # Add `v` to the incoming nodes of `u`
        if u in adjacencies:
            adjacencies[u][1].append(v)
        else:
            adjacencies[u] = ([],[u])
  

# Converting between representations
# ----------------------------------
# 
# Reads each line coming from standard input and updates the `adjacency` dictionary
# using the `add_edge` function.
#
# Once the dictionary is populated, then its contents are printed to standard
# output in the format expected by `gdem`.
def convert(directed = False):
    for line in sys.stdin.readlines():
        edge = line.split()
        if len(edge) != 2:
            raise Exception("edge with a number of nodes different than 2 " + str(edge))
        add_edge(int(edge[0]), int(edge[1]), directed)
    for (node_id, (out_edges, in_edges)) in adjacencies.iteritems():
        sys.stdout.write("%d |" % (node_id))
        for v in out_edges:
            sys.stdout.write(" %d" % (v))
        sys.stdout.write(" |")
        for v in in_edges:
            sys.stdout.write(" %d" % (v))
        sys.stdout.write("\n")


# Main function
# -------------
# 
# The main entry point of the script: sets a command line option (`-d`)
# and invokes `convert`
def main():
    parser = OptionParser()
    parser.add_option("-d", "--directed", 
                      action="store_true",
                      dest="directed", default=False,
                      help="Create a directed graph")

    (options, args) = parser.parse_args()

    sys.stderr.write("Directed graph? " + str(options.directed) + "\n")

    convert(options.directed)


if __name__ == "__main__":
    main()

