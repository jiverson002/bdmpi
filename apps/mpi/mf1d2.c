/*!
\file
\brief A parallel SGD matrix -factorization program using 2D distribution of 
       the matrix
\date Started 6/1/2013
\author George
*/

#define _LARGEFILE64_SOURCE
#include <GKlib.h>
#include <mpi.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>


/**************************************************************************/
/* data structures */
/**************************************************************************/
typedef struct {
  int npes, mype;
  MPI_Comm comm;
  char *filename;

  /* mf parameters */
  int nr, nc;
  int niters;
  int nfactors;
  float lrnrate;
  float alpha;
  float lambda;

  /* timers */
  double totalTmr;
  double loadTmr;
  double setupTmr;
  double compTmr;
  double commTmr;
} params_t;


/* distributed CSR */
typedef struct {
  /* the total number of rows/non-zeros and their overall distribution */
  int gnrows, gncols;
  size_t gnnz;

  gk_csr_t **mats;

  double rowsums; 
} dcsr_t;

/* mf model */
typedef struct {
  int nfactors;
  double mu;
  double *rb, *u;
  double **cbs, **vs;
} mf_t;



/**************************************************************************/
/* prototypes */
/**************************************************************************/
void SetupComms(params_t *params);
dcsr_t *LoadData(params_t *params);
gk_csr_t *CreateBlockMatrix(gk_csr_t *mat, int nr, int nc);
void CleanupData(params_t *params, dcsr_t *dmat);
mf_t *ComputeFactorsRR(params_t *params, dcsr_t *dmat);
void WriteFactors(params_t *params, dcsr_t *dmat, mf_t *model);


/**************************************************************************/
/**************************************************************************/
int main(int argc, char **argv)
{
  params_t *params;
  dcsr_t *dmat;
  mf_t *model;
  double *prvec;
  MPI_Status status;
  double max, current;

  setbuf(stdout, NULL);
  setbuf(stderr, NULL);

  MPI_Init(&argc, &argv);

  params = (params_t *)gk_malloc(sizeof(params_t), "params");
  memset(params, 0, sizeof(params_t));

  params->comm = MPI_COMM_WORLD;
  MPI_Comm_size(params->comm, &(params->npes));
  MPI_Comm_rank(params->comm, &(params->mype));


  if (argc != 6) {
    if (params->mype == 0)
      fprintf(stderr, "Usage: %s filename nr-copies nc-copies niters nfactors\n", argv[0]);
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  params->filename = gk_strdup(argv[1]);
  params->nr       = atoi(argv[2]);
  params->nc       = atoi(argv[3]);
  params->niters   = atoi(argv[4]);
  params->nfactors = atoi(argv[5]);
  params->lrnrate  = 0.001;
  params->alpha    = 0.01;
  params->lambda   = 0.01;

  gk_clearwctimer(params->totalTmr);
  gk_clearwctimer(params->loadTmr);
  gk_clearwctimer(params->setupTmr);
  gk_clearwctimer(params->compTmr);
  gk_clearwctimer(params->commTmr);

  MPI_Barrier(params->comm);
  MPI_Barrier(params->comm);
  gk_startwctimer(params->totalTmr);

  gk_startwctimer(params->loadTmr);
  dmat = LoadData(params);
  gk_stopwctimer(params->loadTmr);

  gk_startwctimer(params->compTmr);
  model = ComputeFactorsRR(params, dmat);
  gk_stopwctimer(params->compTmr);

  CleanupData(params, dmat);

  MPI_Barrier(params->comm);
  MPI_Barrier(params->comm);
  gk_stopwctimer(params->totalTmr);

  /* print timing stats */
  current = gk_getwctimer(params->loadTmr);
  MPI_Reduce(&current, &max, 1, MPI_DOUBLE, MPI_MAX, 0, params->comm);
  if (params->mype == 0)
    printf("  loadTmr:  %10.4lf\n", max);

  current = gk_getwctimer(params->setupTmr);
  MPI_Reduce(&current, &max, 1, MPI_DOUBLE, MPI_MAX, 0, params->comm);
  if (params->mype == 0)
    printf(" setupTmr:  %10.4lf\n", max);

  current = gk_getwctimer(params->compTmr);
  MPI_Reduce(&current, &max, 1, MPI_DOUBLE, MPI_MAX, 0, params->comm);
  if (params->mype == 0)
    printf("  compTmr:  %10.4lf\n", max);

  current = gk_getwctimer(params->commTmr);
  MPI_Reduce(&current, &max, 1, MPI_DOUBLE, MPI_MAX, 0, params->comm);
  if (params->mype == 0)
    printf("  commTmr:  %10.4lf\n", max);

  current = gk_getwctimer(params->totalTmr);
  MPI_Reduce(&current, &max, 1, MPI_DOUBLE, MPI_MAX, 0, params->comm);
  if (params->mype == 0)
    printf(" totalTmr:  %10.4lf\n", max);

  MPI_Finalize();

  return EXIT_SUCCESS;
}


/**************************************************************************/
/*! Reads a sparse matrix in binary CSR format, one process at a time. 
    \returns the local portion of the matrix.
*/
/**************************************************************************/
dcsr_t *LoadData(params_t *params)
{
  int mype=params->mype, npes=params->npes;
  size_t i, p, gnnz, lnnz;
  ssize_t rsize;
  int fd, gnrows, gncols, lnrows;
  dcsr_t *dmat=NULL;
  gk_csr_t *mat;

  if (mype == 0) {
    if (!gk_fexists(params->filename)) 
      errexit("File %s does not exist!\n", params->filename);
  }

  dmat = (dcsr_t *)gk_malloc(sizeof(dcsr_t), "dmat");
  memset(dmat, 0, sizeof(dcsr_t));

  mat = gk_csr_Read(params->filename, GK_CSR_FMT_BINROW, 1, 0);

  /* create the block matrices */
  dmat->mats = (gk_csr_t **)gk_malloc(npes*sizeof(gk_csr_t *), "mats");
  for (p=0; p<npes; p++)
    dmat->mats[p] = CreateBlockMatrix(mat, params->nr, params->nc);

  gk_csr_Free(&mat);

  dmat->gnrows = npes*dmat->mats[0]->nrows;
  dmat->gncols = npes*dmat->mats[0]->ncols;
  dmat->gnnz   = npes*npes*dmat->mats[0]->rowptr[dmat->mats[0]->nrows];

  dmat->rowsums = npes*gk_fsum(dmat->mats[0]->rowptr[dmat->mats[0]->nrows], dmat->mats[0]->rowval, 1);

  return dmat;
}


/**************************************************************************/
/*! This function creates a matrix by pasting together nr*nc copies of the
    input matrix. */
/**************************************************************************/
gk_csr_t *CreateBlockMatrix(gk_csr_t *mat, int nr, int nc)
{
  ssize_t ir, ic, i, j, k, nrows, ncols, bnrows, bncols, nnz, bnnz;
  ssize_t *rowptr, *browptr;
  int *rowind, *browind;
  float *rowval, *browval;
  gk_csr_t *bmat;

  nrows  = mat->nrows;
  ncols  = mat->ncols;
  rowptr = mat->rowptr;
  rowind = mat->rowind;
  rowval = mat->rowval;
  nnz    = rowptr[nrows];


  bmat    = gk_csr_Create();
  bnnz    = nr*nc*nnz;
  bnrows  = bmat->nrows = nr*nrows;
  bncols  = bmat->ncols = nr*ncols;
  browptr = bmat->rowptr = gk_zmalloc(bnrows+1, "browptr");
  browind = bmat->rowind = gk_imalloc(bnnz, "browind");
  browval = bmat->rowval = gk_fmalloc(bnnz, "browval");


  browptr[0] = bnnz= 0;
  for (ir=0; ir<nr; ir++) {
    for (i=0; i<nrows; i++) {
      for (j=rowptr[i]; j<rowptr[i+1]; j++) {
        for (ic=0; ic<nc; ic++, bnnz++) {
          browind[bnnz] = ic*ncols + rowind[j];
          browval[bnnz] = rowval[j];
        }
      }
      browptr[ir*nrows+i+1] = bnnz;
    }
  }

  return bmat;

}


/**************************************************************************/
/*! This function deallocates all the memory that was used */
/**************************************************************************/
void CleanupData(params_t *params, dcsr_t *dmat)
{
  int i;

  for (i=0; i<params->npes; i++)
    gk_csr_Free(&(dmat->mats[i]));

  gk_free((void **)&dmat->mats, &dmat, LTERM);

  return;
}


/**************************************************************************/
/*! This function performs SGD with row-wise randomized traversal */
/**************************************************************************/
mf_t *ComputeFactorsRR(params_t *params, dcsr_t *dmat)
{
  int npes=params->npes, mype=params->mype;
  size_t i, j, iter, block;
  int ir, ic, k, nrows, ncols, nfactors, ldu, ldv, cblock;
  int *rowind;
  ssize_t *rowptr;
  float *rowval;
  mf_t *model;
  double *u, *v, *rb, *cb, uir[params->nfactors];
  double lrnrate, lambda, alpha, mu, r, mae, rmse, gmae, grmse;
  MPI_Status status;
  MPI_Request requests[2];

  model = (mf_t *)gk_malloc(sizeof(mf_t), "model");
  memset(model, 0, sizeof(mf_t));

  nfactors = params->nfactors;
  lrnrate  = params->lrnrate;
  alpha    = params->alpha;
  lambda   = params->lambda;

  nrows  = dmat->mats[0]->nrows;
  ncols  = dmat->mats[0]->ncols;

  /* compute model->mu */
  MPI_Allreduce(&dmat->rowsums, &model->mu, 1, MPI_DOUBLE, MPI_SUM, params->comm);
  mu = model->mu = model->mu/dmat->gnnz;

  /* setup ldu/ldv */
  ldu = ldv = nfactors;

  /* setup cbs/vs */
  model->cbs = (double **)gk_malloc(sizeof(double **)*npes, "cbs");
  model->vs  = (double **)gk_malloc(sizeof(double **)*npes, "vs");
  memset(model->cbs, 0, sizeof(double **)*npes);
  memset(model->vs, 0, sizeof(double **)*npes);
  for (i=0; i<npes; i++)
    model->cbs[i] = model->vs[i] = NULL;

  gk_isrand(101+params->mype);

    
  /* get into the sgd iterations */
  rb = cb = u = v = NULL; 
  for (iter=0; iter<params->niters; iter++) {
    if (mype == 0)
      printf("Working on iteration: %zu\n", iter);

    mae = rmse = 0.0;
    for (block=0; block<npes; block++) {
      cblock = (mype+block)%npes;

      printf("[%3d] Computing %3d during %zu.%zu [ts: %d] SP.\n", 
          mype, cblock, iter, block, (int)time(NULL));

      rowptr = dmat->mats[cblock]->rowptr;
      rowind = dmat->mats[cblock]->rowind;
      rowval = dmat->mats[cblock]->rowval;

      /* if first time, allocate and initialize */
      if (iter == 0 && model->rb == NULL) { 
        rb = model->rb = gk_dmalloc(nrows, "rb");
        u  = model->u  = gk_dmalloc(nrows*nfactors, "u");
        for (ir=0; ir<nrows; ir++) {
          rb[ir] = (RandomInRange(20000)-10000)*0.000001;
          for (k=0; k<nfactors; k++)
            u[ir*ldu+k] = (RandomInRange(20000)-10000)*0.000001;
        }
      }

      /* if first time, allocate and initialize */
      if (iter == 0 && model->cbs[cblock] == NULL) { 
        cb = model->cbs[cblock] = gk_dmalloc(ncols, "cb");
        v  = model->vs[cblock]  = gk_dmalloc(ncols*nfactors, "v");
        for (ic=0; ic<ncols; ic++) {
          cb[ic] = (RandomInRange(20000)-10000)*0.000001;
          for (k=0; k<nfactors; k++)
            v[ic*ldv+k] = (RandomInRange(20000)-10000)*0.000001;
        }
      }
      cb = model->cbs[cblock];
      v  = model->vs[cblock];

      printf("[%3d] Computing %3d during %zu.%zu [ts: %d] P1. [nnz: %zd]\n", 
          mype, cblock, iter, block, (int)time(NULL), rowptr[nrows]);

      /* perform SGD with random row sampling with replacement */
      for (i=0; i<nrows; i++) {
        ir = RandomInRange(nrows);

        for (j=rowptr[ir]; j<rowptr[ir+1]; j++) {
          ic = rowind[j];

          /* save U[ir] since we are updating it */
          for (k=0; k<nfactors; k++)
            uir[k] = u[ir*ldu+k];

          /* compute the error of the prediction */
          r = rowval[j] - mu - rb[ir] - cb[ic];
          for (k=0; k<nfactors; k++)
            r -= uir[k]*v[ic*ldv+k];

          /* aggregate the local mae/rmse */
          mae += fabs(r);
          rmse += r*r;

          /* update the factors */
          rb[ir] += lrnrate*(r - alpha*rb[ir]);
          cb[ic] += lrnrate*(r - alpha*cb[ic]);

          for (k=0; k<nfactors; k++)
            u[ir*ldu+k] += lrnrate*(r*v[ic*ldv+k] - lambda*uir[k]);

          for (k=0; k<nfactors; k++)
            v[ic*ldv+k] += lrnrate*(r*uir[k] - lambda*v[ic*ldv+k]);
        }
      }

      printf("[%3d] Computing %3d during %zu.%zu [ts: %d] P2.\n", 
          mype, cblock, iter, block, (int)time(NULL));


      /* send col model data up using non-blocking send */
      MPI_Isend(cb, ncols, MPI_DOUBLE, (mype-1+npes)%npes, 1, params->comm, &requests[0]);
      MPI_Isend(v, ncols*nfactors, MPI_DOUBLE, (mype-1+npes)%npes, 1, params->comm, &requests[1]);

      printf("[%3d] Computing %3d during %zu.%zu [ts: %d] EP.\n", 
          mype, cblock, iter, block, (int)time(NULL));

      /* receive the column block from down */
      cblock = (cblock+1)%npes;
      model->cbs[cblock] = gk_dmalloc(ncols, "cb");
      model->vs[cblock]  = gk_dmalloc(ncols*nfactors, "v");
      MPI_Recv(model->cbs[cblock], ncols, MPI_DOUBLE, (mype+1)%npes, 1, params->comm, &status);
      MPI_Recv(model->vs[cblock], ncols*nfactors, MPI_DOUBLE, (mype+1)%npes, 1, params->comm, &status);

      MPI_Wait(&requests[0], MPI_STATUS_IGNORE);
      MPI_Wait(&requests[1], MPI_STATUS_IGNORE);
      cblock = (mype+block)%npes;
      gk_free((void **)&model->cbs[cblock], &model->vs[cblock], LTERM);
      cb = v = NULL;
    }

    MPI_Reduce(&mae, &gmae, 1, MPI_DOUBLE, MPI_SUM, 0, params->comm);
    MPI_Reduce(&rmse, &grmse, 1, MPI_DOUBLE, MPI_SUM, 0, params->comm);
    if (mype == 0) {
      gmae = gmae/dmat->gnnz;
      grmse = sqrt(grmse/dmat->gnnz);
      printf("Iter %4zu: MAE: %3.5lf RMSE: %3.5lf\n", iter, gmae, grmse);
    }
  }

  return model;
}
