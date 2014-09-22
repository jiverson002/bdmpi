#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>

#include <mpi.h>

#define DEBUG   0
#define MAXLINE 256

typedef struct {
  int     m, n, nnz;
  int     rs, re;
  int     * rnnz;
  int     * colptr;
  int     * rowind;
  double  * values;
} csr_t;

typedef struct {
  int     m;
  int     * l2g;
  double  * local;
  double  * nonlocal;
  double  * values;
} vec_t;

typedef struct {
  int n_nbrs;
  int * nbr_ids;
  int * ptr;
  int * ind;
} cinfo_t;

typedef struct {
  int k;
  int v;
} kv_t;

double gk_clockseconds();

double ji_memusage(int sz);

int min(const int a, const int b);

int cmp0 (const void * a, const void * b);

int cmp1 (const void * a, const void * b);

void read_files(const char * f1, const char * f2, csr_t * mat, vec_t * vec);

void write_file(const char * f1, vec_t * vec);

void scatter_data(csr_t * imat, csr_t * mat, vec_t * ivec, vec_t * vec, int np, 
      int rank);

void setup_rinfo(cinfo_t * rinfo, int m, int nnz, int * ineed, int np, 
      int rank);

void setup_sinfo(cinfo_t * sinfo, cinfo_t * rinfo, int np, int rank);

void exchange_nonlocal(int m, double * local_val, double ** recv_val_, 
      cinfo_t * sinfo, cinfo_t * rinfo, int np, int rank);

void cat(int mm, int nn, double * local, double * nonlocal, double ** values_,
      int np, int rank);

void g2l(int m, int nnz, int * ineed, cinfo_t * rinfo, int ** l2g_, int np, 
      int rank);

void l2g(int n, int * l2g, int * ineed, int np, int rank);

void sparsemv(csr_t * mat, vec_t * bvec, vec_t * yvec, int np, int rank);

void gather_data(int m, vec_t * yvec, vec_t * ovec, int np, int rank);

#define clear_timer(tmr) (tmr = 0.0)
#define start_timer(tmr) (tmr -= gk_clockseconds())
#define stop_timer(tmr)  (tmr += gk_clockseconds())

