#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include "hll_counter.h"
#include "parser.h"
#include "null.h"
#include "debug.h"
#include "check_ptr.h"
#include "hyper_anf_serial.h"
#include <sys/time.h>

#define COUNTER_BITS 1
#define MAX_ITER 300


void timeval_print(struct timeval *tv)
{
    char buffer[30];
    time_t curtime;

    printf("%ld.%06ld", tv->tv_sec, tv->tv_usec);
    curtime = tv->tv_sec;
    strftime(buffer, 30, "%m-%d-%Y  %T", localtime(&curtime));
    printf(" = %s.%06ld\n", buffer, tv->tv_usec);
}

int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1)
{
    long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);
    result->tv_sec = diff / 1000000;
    result->tv_usec = diff % 1000000;

    return (diff<0);
}

int main(int argc, char **argv) {

  if (argc != 2) {
    fprintf(stderr, "Please provide a graph file name\n");
    exit(EXIT_FAILURE);
  }

  int n = 0;
  node_t *nodes;

  int rc = parse_graph_file(argv[1], &nodes, &n);
  if (rc < 0) {
    fprintf(stderr, "ERROR: could not parse graph file correctly.\n");
    exit(EXIT_FAILURE);
  } else if (rc != n) {
    fprintf(stderr, "ERROR: Some nodes were not parsed correctly, exiting.\n");
    exit(EXIT_FAILURE);
  }

  printf("Loaded graph with %d nodes\n", n);

  struct timeval tvBegin, tvEnd, tvDiff;

  gettimeofday(&tvBegin, NULL);

  int diam = diameter(nodes, n, COUNTER_BITS, MAX_ITER);

  gettimeofday(&tvEnd, NULL);
  timeval_subtract(&tvDiff, &tvEnd, &tvBegin);

  printf("Time for the computation: %03ld.%06ld\n",
         tvDiff.tv_sec, tvDiff.tv_usec);

  printf("Diameter is %d\n", diam);

  // free resources
  for(int i=0; i<n; ++i) {
    node_free(&nodes[i]);
  }
  free(nodes);

  return 0;
}
