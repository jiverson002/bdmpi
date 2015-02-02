#include "hyper_anf_mpi.h"
#include "check_ptr.h"
#include "graph.h"
#include "parser.h"
#include "debug.h"
#include <mpi.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>

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

void init_context ( context_t *context,
                    node_t *partial_graph,
                    size_t partial_graph_cardinality,
                    int bits,
                    int max_iteration,
                    double alpha)
{
  context->alpha = alpha;
  context->num_nodes = partial_graph_cardinality;
  assert(context->num_nodes > 0);
  context->iteration = 0;
  context->num_changed = INT_MAX;
  context->max_iteration = max_iteration;
  context->num_registers = 1 << bits; // 2^bits
  // get the number of processors
  MPI_Comm_size(MPI_COMM_WORLD, & context->num_processors);
  // get the rank of this processor
  MPI_Comm_rank(MPI_COMM_WORLD, & context->rank);

  context->nodes = partial_graph;

  // init neighbourhoods
  context->neighbourhoods =
        malloc(context->num_nodes * sizeof(mpi_neighbourhood_t));
  check_ptr(context->neighbourhoods);
  for (int i=0; i<context->num_nodes; ++i) {
    mpi_neighbourhood_init(
          &context->neighbourhoods[i], partial_graph[i].num_out, bits);
  }

  // init neighbourhood function
  context->neighbourhood_function =
      malloc(context->max_iteration * sizeof(double));
  check_ptr(context->neighbourhood_function);
  memset(context->neighbourhood_function, 0, context->max_iteration * sizeof(double));

  // init counters
  context->counters = malloc(context->num_nodes * sizeof(hll_counter_t));
  check_ptr(context->counters);
  context->counters_prev = malloc(context->num_nodes * sizeof(hll_counter_t));
  check_ptr(context->counters_prev);
  for (int i = 0; i < context->num_nodes; ++i) {
    hll_cnt_init(&context->counters[i], bits);
    hll_cnt_init(&context->counters_prev[i], bits);
    // add the node itself to each counter
    hll_cnt_add(partial_graph[i].id, &context->counters[i]);
    hll_cnt_add(partial_graph[i].id, &context->counters_prev[i]);
  }

  // init array of requests.
  size_t sum = 0;
  for (int i = 0; i < context->num_nodes; ++i) {
    sum += context->nodes[i].num_out + context->nodes[i].num_in;
  }
  context->num_requests = sum;
  context->requests = malloc(sum * sizeof(MPI_Request));
  check_ptr(context->requests);
}

void free_context (context_t *context) {
  for (int i = 0; i < context->num_nodes; ++i) {
    hll_cnt_free(& context->counters[i]);
    hll_cnt_free(& context->counters_prev[i]);
    mpi_neighbourhood_free(& context->neighbourhoods[i]);
  }
  free(context->counters);
  free(context->counters_prev);
  free(context->neighbourhoods);
  free(context->requests);
  free(context->neighbourhood_function);
}

void receive_counters (context_t * context, int i, size_t * request_idx) {
  size_t recv_req_idx = 0;
  context_t ctx = *context;
  mpi_neighbourhood_t neighbourhood = ctx.neighbourhoods[i];

  for (int j = 0; j < ctx.nodes[i].num_out; ++j) {
    node_id_t neighbour = ctx.nodes[i].out[j];
    int neighbour_processor = get_processor_rank(
          neighbour, ctx.num_processors);

    MPI_Irecv( neighbourhood.counters[j].registers, // the buffer of data that receivs the result
               ctx.num_registers, // the number of data items being sent
               MPI_UNSIGNED_CHAR, // the type of data being sent
               neighbour_processor, // the source id.
               neighbour, // the tag of message: the neighbour's ID
               MPI_COMM_WORLD, // the communicator
               & ctx.requests[recv_req_idx++] // the requests to store
             );
  } // end receive from neghbours
  *request_idx += recv_req_idx;
}

void send_counters (context_t * context, int i, size_t * request_idx) {
  context_t ctx = *context;
  size_t snd_req_idx = 0;

  hll_counter_t node_counter = ctx.counters[i];
  hll_reg_t *node_registers = node_counter.registers;

  // this ID will be the tag of the message.
  node_id_t node_id = ctx.nodes[i].id;

  for (int j = 0; j < ctx.nodes[i].num_in; ++j) {
    node_id_t neighbour_id = ctx.nodes[i].in[j];
    int neighbour_processor = get_processor_rank(
          neighbour_id, ctx.num_processors);
    MPI_Isend( node_registers, // the buffer being sent
               ctx.num_registers, // the number of elements being sent
               MPI_UNSIGNED_CHAR, // the type of message elements
               neighbour_processor, // the destination of the message
               node_id, // the tag of the message: the node's ID
               MPI_COMM_WORLD, // the communicator
               & ctx.requests[snd_req_idx++] // the request to store
             );
  } // end send to neighbours
  *request_idx += snd_req_idx;
}

void exchange_counters (context_t * context) {

  struct timeval tvBegin, tvEnd, tvDiff;

  size_t request_idx = 0;
  // - for each node in partial graph
  printf("(Process %d) Distributing counters\n", context->rank);

  gettimeofday(&tvBegin, NULL);

  for (int i = 0; i < context->num_nodes; ++i) {
    //   * expect to receive counters from out neighbours
    receive_counters(context, i, &request_idx);
    //   * send counter to in neighbours
    send_counters(context, i, &request_idx);
  }

  gettimeofday(&tvEnd, NULL);
  timeval_subtract(&tvDiff, &tvEnd, &tvBegin);

  // check if we have sent all requests
  printd("Requests %d over %d\n", request_idx, context->num_requests);

  assert(request_idx == context->num_requests);

  printf("(Process %d) Time for the counters exchange: %ld.%06ld\n",
         context->rank, tvDiff.tv_sec, tvDiff.tv_usec);

}

int update_counters (context_t *context) {
  context_t ctx = *context;

  printf("(Process %d) Updating counters\n", ctx.rank);
  // set up the counter of changed nodes
  int changed_counters = 0;
  // - for each node in partial graph
  for (int i = 0; i < ctx.num_nodes; ++i) {
    hll_counter_t
        node_counter = ctx.counters[i],
        node_counter_prev = ctx.counters_prev[i];
    assert(hll_cnt_size(&node_counter) >= hll_cnt_size(&node_counter_prev));
    //   * update counters
    for (int j = 0; j < ctx.neighbourhoods[i].dimension; ++j) {
      hll_cnt_union_i(
            &node_counter, &ctx.neighbourhoods[i].counters[j]);
    }
    if(!hll_cnt_equals(&node_counter, &node_counter_prev)) {
      ++changed_counters;
    }
    // update the prev counter
    hll_cnt_copy_to(&node_counter, &node_counter_prev);
  }
  return changed_counters;
}

void count_changed (context_t * context, int local_changed) {
  struct timeval tvBegin, tvEnd, tvDiff;
  gettimeofday(&tvBegin, NULL);

  int changed_counters = local_changed;
  printf("(Process %d) %d counters changed\n",
         context->rank, changed_counters);
  // - use mpi_reduce to compute the number of changed nodes.
  int total_changed = 0;
  MPI_Reduce( &changed_counters,
              &total_changed,
              1,
              MPI_INT,
              MPI_SUM,
              0, // root of the sum
              MPI_COMM_WORLD);

  if(context->rank == 0) {
    printf("(Process %d) A total of %d counters changed\n",
           context->rank, total_changed);
  }
  MPI_Bcast(&total_changed, 1, MPI_INT, 0, MPI_COMM_WORLD);

  context->num_changed = total_changed;

  gettimeofday(&tvEnd, NULL);
  timeval_subtract(&tvDiff, &tvEnd, &tvBegin);

  printf("(Process %d) Time for the counters counting: %ld.%06ld\n",
         context->rank, tvDiff.tv_sec, tvDiff.tv_usec);

}

// **Attention**: maybe it's not really important to tag messages. When we
// perform the union of counters we don't care from where the counters arrive.
double mpi_diameter( context_t * context )
{
  struct timeval tvBegin, tvEnd, tvDiff;

  gettimeofday(&tvBegin, NULL);

  while ( context->num_changed != 0 &&
          context->iteration < context->max_iteration)
  {
    printd("(Process %d) Start iteration %d\n",
           context->rank, context->iteration);
    // reset the number of changed nodes
    context->num_changed = 0;
    // compute the neighbourhood function before updating counters.
    compute_neighbourhood_function(context);

    exchange_counters(context);

    int changed_counters = update_counters(context);

    count_changed(context, changed_counters);

    printd("(Process %d) finish iteration %d\n",
           context->rank, context->iteration);

    ++context->iteration;
  }

  gettimeofday(&tvEnd, NULL);
  timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
  printf("(Process %d) Elapsed time: %ld.%06ld\n",
         context->rank, tvDiff.tv_sec, tvDiff.tv_usec);

  return effective_diameter(context);
}

double effective_diameter(context_t *context) {
  int last_idx = context->iteration - 1;
  // broadcast neighbourhood function
  MPI_Bcast(context->neighbourhood_function, last_idx+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  double nf_max = context->neighbourhood_function[last_idx];
  printf("(Process %d) max value = %f\n", context->rank, nf_max);
  double nf_elem = 0;
  int i = 0;
  for(; i<=last_idx; ++i) {
    nf_elem = context->neighbourhood_function[i];
//    printf("(Process %d) N(%d) = %f\n", context->rank, i, nf_elem);
    if(nf_elem/nf_max >= context->alpha) {
      break;
    }
  }

  if (i == 0) {
    return i + (context->alpha * nf_max - nf_elem) / nf_elem;
  } else {
    double num = (context->alpha * nf_max - nf_elem);
    double den = (nf_elem - context->neighbourhood_function[i-1]);
//    printf("(Process %d) num = %f, den = %f\n", num, den);
    return i + num / den;
  }
}

void compute_neighbourhood_function (context_t * context) {
  hll_cardinality_t local_sum = 0;
  for (int i = 0; i < context->num_nodes; ++i) {
    local_sum += hll_cnt_size(&context->counters[i]);
  }
  hll_cardinality_t sum;
  MPI_Reduce(&local_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if ( context->rank == 0 ) {
    printf("N(%d) = %f\n", context->iteration, sum);
    context->neighbourhood_function[context->iteration] = sum;
  }
}


inline int get_processor_rank (node_id_t node, int num_processors) {
  return node % num_processors;
}

void mpi_neighbourhood_init (mpi_neighbourhood_t *neigh, int n, int bits) {
  neigh->dimension = n;
  neigh->counters = malloc(n*sizeof(hll_counter_t));
  check_ptr(neigh->counters);
  for (int i = 0; i < n; ++i) {
    hll_cnt_init(&neigh->counters[i], bits);
  }
}

void mpi_neighbourhood_free (mpi_neighbourhood_t *neigh) {
  // FIXME: it seems that this free is invalid
//  for (int i = 0; i < neigh->dimension; ++i) {
//    hll_cnt_free(&neigh->counters[i]);
//  }
  free(neigh->counters);
}

#define MAX_FILENAME_LENGTH 128

void load_partial_graph ( int rank,
                          char *basename,
                          node_t **nodes,
                          int *n) {
  char filename[MAX_FILENAME_LENGTH];
  sprintf(filename, "%s-%d.adj", basename, rank);
  printf("(Process %d) Loading file %s\n", rank, filename);
  parse_graph_file(filename, nodes, n);
}
