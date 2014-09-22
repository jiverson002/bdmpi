/*!
\file
\brief A parallel SGD matrix -factorization program using 2D distribution of 
       the matrix
\date Started 6/1/2013
\author George
*/


#define _LARGEFILE64_SOURCE
#include <GKlib.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>


/**************************************************************************/
/* data structures */
/**************************************************************************/
typedef struct {
  int npes, nrowpes, ncolpes;
  char *filename;

  /* out-of-core files */
  char **mfiles;
  char **ufiles;
  char **vfiles;
  char **rbfiles;
  char **cbfiles;

  /* mf parameters */
  int type;
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

  /* some additional info that gets computed prior to storing the data to disk */
  int lnrows, lncols;
  size_t lnnz;
  double rvalsum;

  gk_csr_t *mat;
} dcsr_t;

/* mf model */
typedef struct {
  int nfactors;
  double mu;
  double *u, *v, *rb, *cb;
} mf_t;


/**************************************************************************/
/* prototypes */
/**************************************************************************/
dcsr_t *LoadData2D(params_t *params);
void CleanupData(params_t *params, dcsr_t *dmat);
mf_t *ComputeFactorsFR(params_t *params, dcsr_t *dmat);
mf_t *ComputeFactorsRR(params_t *params, dcsr_t *dmat);
void WriteFactors(params_t *params, dcsr_t *dmat, mf_t *model);


/**************************************************************************/
/**************************************************************************/
int main(int argc, char **argv)
{
  int pe;
  params_t *params;
  dcsr_t *dmat;
  mf_t *model;
  double *prvec;
  double max, current;
  char string[8192];

  setbuf(stdout, NULL);
  setbuf(stderr, NULL);

  params = (params_t *)gk_malloc(sizeof(params_t), "params");
  memset(params, 0, sizeof(params_t));


  if (argc != 7) {
    fprintf(stderr, "Usage: %s filename nrowpes ncolpes niters nfactors type(1:FR;2:RR)\n", argv[0]);
    return EXIT_FAILURE;
  }

  params->filename = gk_strdup(argv[1]);
  params->nrowpes  = atoi(argv[2]);
  params->ncolpes  = atoi(argv[3]);
  params->niters   = atoi(argv[4]);
  params->nfactors = atoi(argv[5]);
  params->type     = atoi(argv[6]);
  params->npes     = params->nrowpes*params->ncolpes;
  params->lrnrate  = 0.001;
  params->alpha    = 0.01;
  params->lambda   = 0.01;


  params->mfiles  = (char **)gk_malloc(params->npes*sizeof(char *), "mfiles");
  params->ufiles  = (char **)gk_malloc(params->npes*sizeof(char *), "ufiles");
  params->vfiles  = (char **)gk_malloc(params->npes*sizeof(char *), "vfiles");
  params->rbfiles = (char **)gk_malloc(params->npes*sizeof(char *), "rbfiles");
  params->cbfiles = (char **)gk_malloc(params->npes*sizeof(char *), "cbfiles");
  for (pe=0; pe<params->npes; pe++) {
    sprintf(string,  "m%d.%d", (int)getpid(), pe);
    params->mfiles[pe] = gk_strdup(string);
    sprintf(string,  "u%d.%d", (int)getpid(), pe);
    params->ufiles[pe] = gk_strdup(string);
    sprintf(string,  "v%d.%d", (int)getpid(), pe);
    params->vfiles[pe] = gk_strdup(string);
    sprintf(string,  "cb%d.%d", (int)getpid(), pe);
    params->cbfiles[pe] = gk_strdup(string);
    sprintf(string,  "rb%d.%d", (int)getpid(), pe);
    params->rbfiles[pe] = gk_strdup(string);
  }

  gk_clearwctimer(params->totalTmr);
  gk_clearwctimer(params->loadTmr);
  gk_clearwctimer(params->setupTmr);
  gk_clearwctimer(params->compTmr);
  gk_clearwctimer(params->commTmr);

  gk_startwctimer(params->totalTmr);

  gk_startwctimer(params->loadTmr);
  dmat = LoadData2D(params);
  gk_stopwctimer(params->loadTmr);

  gk_startwctimer(params->compTmr);
  if (params->type == 1)
    model = ComputeFactorsFR(params, dmat);
  else
    model = ComputeFactorsRR(params, dmat);
  gk_stopwctimer(params->compTmr);

  CleanupData(params, dmat);

  gk_stopwctimer(params->totalTmr);

  /* print timing stats */
  max = gk_getwctimer(params->loadTmr);
  printf("  loadTmr:  %10.4lf\n", max);

  max = gk_getwctimer(params->setupTmr);
  printf(" setupTmr:  %10.4lf\n", max);

  max = gk_getwctimer(params->compTmr);
  printf("  compTmr:  %10.4lf\n", max);

  max = gk_getwctimer(params->commTmr);
  printf("  commTmr:  %10.4lf\n", max);

  max = gk_getwctimer(params->totalTmr);
  printf(" totalTmr:  %10.4lf\n", max);


  return EXIT_SUCCESS;
}


/**************************************************************************/
/*! Reads a sparse matrix in binary CSR format, one process at a time. 
    This version just reads the entire input file and treats it as a local
    block with the rows/columns getting their "natural" positions in the
    overall matrix.
    \returns the local portion of the matrix.
*/
/**************************************************************************/
dcsr_t *LoadData2D(params_t *params)
{
  int pe, myrow, mycol, npes=params->npes, nrowpes=params->nrowpes, ncolpes=params->ncolpes;
  int lrank, lsize, token=1, lnrows, lncols;
  size_t i, j, iR, jR, lnnz;
  dcsr_t *dmat=NULL;


  if (!gk_fexists(params->filename)) 
    errexit("File %s does not exist!\n", params->filename);

  dmat = (dcsr_t *)gk_malloc(sizeof(dcsr_t), "dmat");
  memset(dmat, 0, sizeof(dcsr_t));

  dmat->rvalsum = 0.0;

  /* read the matrix */
  for (pe=0; pe<npes; pe++) {
    myrow = npes/ncolpes;
    mycol = npes%ncolpes;

    dmat->mat = gk_csr_Read(params->filename, GK_CSR_FMT_BINROW, 1, 0);
    lnrows = dmat->mat->nrows;
    lncols = dmat->mat->ncols;
    lnnz   = dmat->mat->rowptr[lnrows];

    dmat->gnrows = lnrows*nrowpes;
    dmat->gncols = lncols*ncolpes;
    dmat->gnnz   = lnnz*npes;

    dmat->rvalsum += gk_fsum(dmat->lnnz, dmat->mat->rowval, 1);

    printf("[%3d, %3d] dmat->gnrows/lnrows: %d/%d, dmat->gncols/lncols: %d/%d, dmat->gnnz/lnnz: %zu/%zu\n", 
        myrow, mycol, 
        dmat->gnrows, dmat->mat->nrows, 
        dmat->gncols, dmat->mat->ncols, 
        dmat->gnnz, dmat->mat->rowptr[dmat->mat->nrows]);

    gk_csr_Write(dmat->mat, params->mfiles[pe], GK_CSR_FMT_BINROW, 1, 0);
    gk_csr_Free(&(dmat->mat));
  }

  return dmat;
}


/**************************************************************************/
/*! This function deallocates all the memory that was used */
/**************************************************************************/
void CleanupData(params_t *params, dcsr_t *dmat)
{
  int pe;

  gk_csr_Free(&(dmat->mat));
  gk_free((void **)&dmat, LTERM);

  for (pe=0; pe<params->npes; pe++) {
    unlink(params->mfiles[pe]);
    unlink(params->ufiles[pe]);
    unlink(params->vfiles[pe]);
    unlink(params->rbfiles[pe]);
    unlink(params->cbfiles[pe]);
  }

  return;
}


/**************************************************************************/
/*! Performs SGD with full randomized nnz traversal */
/**************************************************************************/
mf_t *ComputeFactorsFR(params_t *params, dcsr_t *dmat)
{
  int pe, myrow, mycol, 
      npes=params->npes, nrowpes=params->nrowpes, ncolpes=params->ncolpes;
  size_t i, j, iter, block, nsamples, nelmnts;
  int ir, ic, k, nrows, ncols, nfactors, ldu, ldv;
  int *rowind;
  ssize_t *rowptr;
  float *rowval;
  mf_t *model;
  double *u=NULL, *v=NULL, *rb=NULL, *cb=NULL, uir[params->nfactors];
  double lrnrate, lambda, alpha, mu, r, mae, rmse;


  model = (mf_t *)gk_malloc(sizeof(mf_t), "model");
  memset(model, 0, sizeof(mf_t));

  nfactors = params->nfactors;
  lrnrate  = params->lrnrate;
  alpha    = params->alpha;
  lambda   = params->lambda;
  nsamples = dmat->gnnz/npes;


  /* compute model->mu */
  mu = model->mu = dmat->rvalsum/dmat->gnnz;

  /* setup ldu/ldv */
  ldu = ldv = nfactors;

  gk_isrand(101);

    
  /* get into the sgd iterations */
  for (iter=0; iter<params->niters; iter++) {
    printf("Working on iteration: %zu\n", iter);

    mae = rmse = 0.0;
    for (pe=0; pe<npes; pe++) {
      myrow = pe/ncolpes;
      mycol = pe%ncolpes;

      printf("[%3d %3d] Computing during %zu [ts: %d] SP.\n",
          myrow, mycol, iter, (int)time(NULL));

      dmat->mat = gk_csr_Read(params->mfiles[pe], GK_CSR_FMT_BINROW, 1, 0);
      nrows  = dmat->mat->nrows;
      ncols  = dmat->mat->ncols;
      rowptr = dmat->mat->rowptr;
      rowind = dmat->mat->rowind;
      rowval = dmat->mat->rowval;

      /* if first time, allocate and initialize */
      if (mycol == 0) { 
        if (iter == 0) {
          rb = model->rb = gk_dmalloc(nrows, "rb");
          u  = model->u  = gk_dmalloc(nrows*nfactors, "u");
          for (ir=0; ir<nrows; ir++) {
            rb[ir] = (RandomInRange(20000)-10000)*0.000001;
            for (k=0; k<nfactors; k++)
              u[ir*ldu+k] = (RandomInRange(20000)-10000)*0.000001;
          }
        }
        else { /* read it from the previous iteration */
          rb = model->rb = gk_dreadfilebin(params->rbfiles[myrow], &nelmnts);
          u  = model->u  = gk_dreadfilebin(params->ufiles[myrow], &nelmnts);
          unlink(params->rbfiles[myrow]);
          unlink(params->ufiles[myrow]);
        }
      }

      /* if first time, allocate and initialize */
      if (iter == 0 && myrow == 0 ) { 
        cb = model->cb = gk_dmalloc(ncols, "cb");
        v  = model->v  = gk_dmalloc(ncols*nfactors, "v");
        for (ic=0; ic<ncols; ic++) {
          cb[ic] = (RandomInRange(20000)-10000)*0.000001;
          for (k=0; k<nfactors; k++)
            v[ic*ldv+k] = (RandomInRange(20000)-10000)*0.000001;
        }
      }
      else {
        cb = model->cb = gk_dreadfilebin(params->cbfiles[mycol], &nelmnts);
        v  = model->v  = gk_dreadfilebin(params->vfiles[mycol], &nelmnts);
        unlink(params->cbfiles[mycol]);
        unlink(params->vfiles[mycol]);
      }

      printf("[%3d %3d] Computing during %zu [ts: %d] P1.\n",
          myrow, mycol, iter, (int)time(NULL));

      GKASSERT(rb != NULL);
      GKASSERT(u != NULL);
      GKASSERT(cb != NULL);
      GKASSERT(v != NULL);

      /* perform SGD with random sampling with replacement */
      for (i=0; i<nsamples; i++) {
        ir = RandomInRange(nrows);
        if (rowptr[ir+1]-rowptr[ir] == 0)
          continue;

        j  = rowptr[ir] + RandomInRange(rowptr[ir+1]-rowptr[ir]);
        ic = rowind[j];

        gk_dcopy(nfactors, u+ir*ldu, uir); /* save U[ir] since we are updating it */

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

      printf("[%3d %3d] Computing during %zu [ts: %d] P2.\n",
          myrow, mycol, iter, (int)time(NULL));

      /* free the read-only data */
      gk_csr_Free(&(dmat->mat));

      if (mycol == ncolpes-1) { /* write the row models */
        gk_dwritefilebin(params->rbfiles[myrow], nrows, model->rb);
        gk_dwritefilebin(params->ufiles[myrow], nrows*nfactors, model->u);
        gk_free((void **)&model->rb, &model->u, LTERM);
        rb = u = NULL;
      }

      /* always write the column models as you move along the columns */
      gk_dwritefilebin(params->cbfiles[mycol], ncols, model->cb);
      gk_dwritefilebin(params->vfiles[mycol], ncols*nfactors, model->v);
      gk_free((void **)&model->cb, &model->v, LTERM);
      cb = v = NULL;
    }

    mae = mae/dmat->gnnz;
    rmse = sqrt(rmse/dmat->gnnz);
    printf("Iter %4zu: MAE: %3.5lf RMSE: %3.5lf\n", iter, mae, rmse);
  }

  return model;
}


/**************************************************************************/
/*! Performs SGD with row-wise randomized traversal */
/**************************************************************************/
mf_t *ComputeFactorsRR(params_t *params, dcsr_t *dmat)
{
  int pe, myrow, mycol, 
      npes=params->npes, nrowpes=params->nrowpes, ncolpes=params->ncolpes;
  size_t i, j, iter, block, nelmnts;
  int ir, ic, k, nrows, ncols, nfactors, ldu, ldv;
  int *rowind;
  ssize_t *rowptr;
  float *rowval;
  mf_t *model;
  double *u=NULL, *v=NULL, *rb=NULL, *cb=NULL, uir[params->nfactors];
  double lrnrate, lambda, alpha, mu, r, mae, rmse;


  model = (mf_t *)gk_malloc(sizeof(mf_t), "model");
  memset(model, 0, sizeof(mf_t));

  nfactors = params->nfactors;
  lrnrate  = params->lrnrate;
  alpha    = params->alpha;
  lambda   = params->lambda;


  /* compute model->mu */
  mu = model->mu = dmat->rvalsum/dmat->gnnz;

  /* setup ldu/ldv */
  ldu = ldv = nfactors;

  gk_isrand(101);

    
  /* get into the sgd iterations */
  for (iter=0; iter<params->niters; iter++) {
    printf("Working on iteration: %zu\n", iter);

    mae = rmse = 0.0;
    for (pe=0; pe<npes; pe++) {
      myrow = pe/ncolpes;
      mycol = pe%ncolpes;

      printf("[%3d %3d] Computing during %zu [ts: %d] SP.\n",
          myrow, mycol, iter, (int)time(NULL));

      dmat->mat = gk_csr_Read(params->mfiles[pe], GK_CSR_FMT_BINROW, 1, 0);
      nrows  = dmat->mat->nrows;
      ncols  = dmat->mat->ncols;
      rowptr = dmat->mat->rowptr;
      rowind = dmat->mat->rowind;
      rowval = dmat->mat->rowval;

      /* if first time, allocate and initialize */
      if (mycol == 0) { 
        if (iter == 0) {
          rb = model->rb = gk_dmalloc(nrows, "rb");
          u  = model->u  = gk_dmalloc(nrows*nfactors, "u");
          for (ir=0; ir<nrows; ir++) {
            rb[ir] = (RandomInRange(20000)-10000)*0.000001;
            for (k=0; k<nfactors; k++)
              u[ir*ldu+k] = (RandomInRange(20000)-10000)*0.000001;
          }
        }
        else { /* read it from the previous iteration */
          rb = model->rb = gk_dreadfilebin(params->rbfiles[myrow], &nelmnts);
          u  = model->u  = gk_dreadfilebin(params->ufiles[myrow], &nelmnts);
          unlink(params->rbfiles[myrow]);
          unlink(params->ufiles[myrow]);
        }
      }

      /* if first time, allocate and initialize */
      if (iter == 0 && myrow == 0 ) { 
        cb = model->cb = gk_dmalloc(ncols, "cb");
        v  = model->v  = gk_dmalloc(ncols*nfactors, "v");
        for (ic=0; ic<ncols; ic++) {
          cb[ic] = (RandomInRange(20000)-10000)*0.000001;
          for (k=0; k<nfactors; k++)
            v[ic*ldv+k] = (RandomInRange(20000)-10000)*0.000001;
        }
      }
      else {
        cb = model->cb = gk_dreadfilebin(params->cbfiles[mycol], &nelmnts);
        v  = model->v  = gk_dreadfilebin(params->vfiles[mycol], &nelmnts);
        unlink(params->cbfiles[mycol]);
        unlink(params->vfiles[mycol]);
      }

      printf("[%3d %3d] Computing during %zu [ts: %d] P1.\n",
          myrow, mycol, iter, (int)time(NULL));

      GKASSERT(rb != NULL);
      GKASSERT(u != NULL);
      GKASSERT(cb != NULL);
      GKASSERT(v != NULL);

      /* perform SGD with random row sampling with replacement */
      for (i=0; i<nrows; i++) {
        ir = RandomInRange(nrows);
        if (rowptr[ir+1]-rowptr[ir] == 0)
          continue;

        for (j=rowptr[ir]; j<rowptr[ir+1]; j++) {
          ic = rowind[j];

          gk_dcopy(nfactors, u+ir*ldu, uir); /* save U[ir] since we are updating it */

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

      printf("[%3d %3d] Computing during %zu [ts: %d] P2.\n",
          myrow, mycol, iter, (int)time(NULL));

      /* free the read-only data */
      gk_csr_Free(&(dmat->mat));

      if (mycol == ncolpes-1) { /* write the row models */
        gk_dwritefilebin(params->rbfiles[myrow], nrows, model->rb);
        gk_dwritefilebin(params->ufiles[myrow], nrows*nfactors, model->u);
        gk_free((void **)&model->rb, &model->u, LTERM);
        rb = u = NULL;
      }

      /* always write the column models as you move along the columns */
      gk_dwritefilebin(params->cbfiles[mycol], ncols, model->cb);
      gk_dwritefilebin(params->vfiles[mycol], ncols*nfactors, model->v);
      gk_free((void **)&model->cb, &model->v, LTERM);
      cb = v = NULL;
    }

    mae = mae/dmat->gnnz;
    rmse = sqrt(rmse/dmat->gnnz);
    printf("Iter %4zu: MAE: %3.5lf RMSE: %3.5lf\n", iter, mae, rmse);
  }

  return model;
}

