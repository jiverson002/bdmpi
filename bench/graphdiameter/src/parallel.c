#include "hyper_anf_mpi.h"
#include "graph.h"
#include <stdio.h>
#include "boolean.h"

struct options {
  int bits; /**< The number of bits to use in hll counters */
  int max_iter; /**< The maximum number of iterations */
  double alpha; /**< We compute the effective diameter at alpha */
  char *basename; /**< The base name of the graph under exam. */
};

typedef struct options options_t;

void print_options(options_t opts) {
  printf(
        "Options:\n"
        "  bits: %d\n"
        "  max_iter: %d\n"
        "  alpha: %f\n"
        "  basename: %s\n",
        opts.bits, opts.max_iter, opts.alpha, opts.basename);
}

void print_usage() {

}

options_t parse_options(int argc, char **argv) {
  options_t opts;
  int bits_set = FALSE,
      max_iter_set = FALSE,
      basename_set = FALSE,
      alpha_set = FALSE;
  for (int idx = 1;idx < argc; ++idx) {
    if(argv[idx][0] == '-') {
      switch (argv[idx][1]){
      case 'h':
        print_usage();
        break;
      case 'b': //basename
        opts.basename = argv[++idx];
        basename_set = TRUE;
        break;
      case 'm': // max iteration number
        opts.max_iter = atol(argv[++idx]);
        max_iter_set = TRUE;
        break;
      case 'k': // bits
        opts.bits = atol(argv[++idx]);
        bits_set = TRUE;
        break;
      case 'a': // alpha
        opts.alpha = strtod(argv[++idx], NULL);
        alpha_set = TRUE;
        break;

      default:
        printf("Wrong Argument: %s\n", argv[idx]);
        print_usage();
        exit(1);
      }
    } else {
      printf("Wrong Argument: %s\n", argv[idx]);
      print_usage();
      exit(1);
    }
  }

  if(!bits_set)
    opts.bits = 1;
  if(!max_iter_set)
    opts.max_iter = 10;
  if(!basename_set)
    opts.basename = "graph";
  if(!alpha_set)
    opts.alpha = 1;
  return opts;
}

int main(int argc, char **argv) {

  options_t opts = parse_options(argc, argv);
  print_options(opts);

  int rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  node_t * nodes = 0;
  int n;
  load_partial_graph(rank, opts.basename, &nodes, &n);
  printf("(Process %d) loaded %d nodes\n", rank, n);


  context_t context;
  init_context(&context, nodes, n, opts.bits, opts.max_iter, opts.alpha);

  double diameter = mpi_diameter(&context);

  if(rank == 0) {
    printf("Effective diameter at %f = %f\n", opts.alpha, diameter);
  }

  if(nodes != 0) {
    free(nodes);
  }
  free_context(&context);

  printf("Done!!\n");

  MPI_Finalize();

  return 0;
}
