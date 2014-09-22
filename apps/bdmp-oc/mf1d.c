/*!
\file
\brief A parallel SGD matrix -factorization program using 2D distribution of 
       the matrix
\date Started 6/1/2013
\author George
*/

#define _LARGEFILE64_SOURCE
#include <GKlib.h>
#include <bdmpi.h>
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
  BDMPI_Comm comm;
  char *filename;

  /* out-of-core files */
  char mfile[8192];
  char ufile[8192];
  char vfile[8192];
  char rbfile[8192];
  char cbfile[8192];

  /* mf parameters */
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
  int *rowdist;

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
void SetupComms(params_t *params);
dcsr_t *LoadData1D(params_t *params);
void CleanupData(params_t *params, dcsr_t *dmat);
mf_t *ComputeFactors(params_t *params, dcsr_t *dmat);
void WriteFactors(params_t *params, dcsr_t *dmat, mf_t *model);


/**************************************************************************/
/**************************************************************************/
int main(int argc, char **argv)
{
  params_t *params;
  dcsr_t *dmat;
  mf_t *model;
  double *prvec;
  BDMPI_Status status;
  double max, current;

  setbuf(stdout, NULL);
  setbuf(stderr, NULL);

  BDMPI_Init(&argc, &argv);

  params = (params_t *)gk_malloc(sizeof(params_t), "params");
  memset(params, 0, sizeof(params_t));

  params->comm = BDMPI_COMM_WORLD;
  BDMPI_Comm_size(params->comm, &(params->npes));
  BDMPI_Comm_rank(params->comm, &(params->mype));


  if (argc != 4) {
    if (params->mype == 0)
      fprintf(stderr, "Usage: %s filename niters nfactors\n", argv[0]);
    BDMPI_Finalize();
    return EXIT_FAILURE;
  }

  params->filename = gk_strdup(argv[1]);
  params->niters   = atoi(argv[2]);
  params->nfactors = atoi(argv[3]);
  params->lrnrate  = 0.001;
  params->alpha    = 0.01;
  params->lambda   = 0.01;

  sprintf(params->mfile,  "m%d.%d", (int)getpid(), params->mype);
  sprintf(params->ufile,  "u%d.%d", (int)getpid(), params->mype);
  sprintf(params->vfile,  "v%d.%d", (int)getpid(), params->mype);
  sprintf(params->rbfile, "rb%d.%d", (int)getpid(), params->mype);
  sprintf(params->cbfile, "cb%d.%d", (int)getpid(), params->mype);

  gk_clearwctimer(params->totalTmr);
  gk_clearwctimer(params->loadTmr);
  gk_clearwctimer(params->setupTmr);
  gk_clearwctimer(params->compTmr);
  gk_clearwctimer(params->commTmr);

  BDMPI_Barrier(params->comm);
  BDMPI_Barrier(params->comm);
  gk_startwctimer(params->totalTmr);

  gk_startwctimer(params->loadTmr);
  dmat = LoadData1D(params);
  gk_stopwctimer(params->loadTmr);

  gk_startwctimer(params->compTmr);
  model = ComputeFactors(params, dmat);
  gk_stopwctimer(params->compTmr);

  //WriteFactors(params, dmat, model);

  CleanupData(params, dmat);

  BDMPI_Barrier(params->comm);
  BDMPI_Barrier(params->comm);
  gk_stopwctimer(params->totalTmr);

  /* print timing stats */
  current = gk_getwctimer(params->loadTmr);
  BDMPI_Reduce(&current, &max, 1, BDMPI_DOUBLE, BDMPI_MAX, 0, params->comm);
  if (params->mype == 0)
    printf("  loadTmr:  %10.4lf\n", max);

  current = gk_getwctimer(params->setupTmr);
  BDMPI_Reduce(&current, &max, 1, BDMPI_DOUBLE, BDMPI_MAX, 0, params->comm);
  if (params->mype == 0)
    printf(" setupTmr:  %10.4lf\n", max);

  current = gk_getwctimer(params->compTmr);
  BDMPI_Reduce(&current, &max, 1, BDMPI_DOUBLE, BDMPI_MAX, 0, params->comm);
  if (params->mype == 0)
    printf("  compTmr:  %10.4lf\n", max);

  current = gk_getwctimer(params->commTmr);
  BDMPI_Reduce(&current, &max, 1, BDMPI_DOUBLE, BDMPI_MAX, 0, params->comm);
  if (params->mype == 0)
    printf("  commTmr:  %10.4lf\n", max);

  current = gk_getwctimer(params->totalTmr);
  BDMPI_Reduce(&current, &max, 1, BDMPI_DOUBLE, BDMPI_MAX, 0, params->comm);
  if (params->mype == 0)
    printf(" totalTmr:  %10.4lf\n", max);

//DONE:
  BDMPI_Finalize();

  return EXIT_SUCCESS;
}


/**************************************************************************/
/*! Reads a sparse matrix in binary CSR format, one process at a time. 
    \returns the local portion of the matrix.
*/
/**************************************************************************/
dcsr_t *LoadData1D(params_t *params)
{
  int mype=params->mype, npes=params->npes, token=1;
  int lrank, lsize;
  size_t i, p, gnnz, lnnz;
  ssize_t rsize;
  int fd, gnrows, gncols, lnrows;
  dcsr_t *dmat=NULL;
  BDMPI_Status status;
  off64_t fpos;
  ssize_t *rowptr;

  BDMPI_Comm_lrank(params->comm, &lrank);
  BDMPI_Comm_lsize(params->comm, &lsize);

  if (mype == 0) {
    if (!gk_fexists(params->filename)) 
      errexit("File %s does not exist!\n", params->filename);
  }

  dmat = (dcsr_t *)gk_malloc(sizeof(dcsr_t), "dmat");
  memset(dmat, 0, sizeof(dcsr_t));

  dmat->rowdist = gk_imalloc(npes+1, "rowdist");

  /* root determines the rowdist array so that it balances the lnnz's */
  if (mype == 0) {
    if ((fd = open(params->filename, O_RDONLY)) == -1)
      errexit("Failed opeing the file %s. [%s]\n", params->filename, strerror(errno));
    if (gk_read(fd, &gnrows, sizeof(int)) != sizeof(int))
      errexit("Failed to read the nrows from file %s!\n", params->filename);
    if (gk_read(fd, &gncols, sizeof(int)) != sizeof(int))
      errexit("Failed to read the ncols from file %s!\n", params->filename);

    rowptr = gk_zmalloc(gnrows+1, "rowptr");
    if (gk_read(fd, rowptr, sizeof(ssize_t)*(gnrows+1)) != sizeof(ssize_t)*(gnrows+1))
      errexit("Failed to read the rowptr from file %s!\n", params->filename);
    close(fd);

    /* populate the rowdist */
    dmat->rowdist[0] = 0;
    for (i=0, p=0; p<npes; p++) {
      lnnz = rowptr[i] + (rowptr[gnrows] - rowptr[i] + npes - p - 1)/(npes-p);
      for (; i<gnrows; i++) {
        if (rowptr[i] >= lnnz)
          break;
      }
      dmat->rowdist[p+1] = i;
      //printf("%5zu %10zu %10d %10zu %10zu\n", p, i, gnrows, lnnz, rowptr[gnrows]);
    }

    gk_free((void **)&rowptr, LTERM);
  }

  /* broadcast rowdist */
  BDMPI_Bcast(dmat->rowdist, npes+1, BDMPI_INT, 0, params->comm);


  /* wait your turn */
  if (lrank != 0) 
    BDMPI_Recv(&token, 1, BDMPI_INT, mype-1, 1, params->comm, &status);

  if ((fd = open(params->filename, O_RDONLY)) == -1)
    errexit("Failed opeing the file %s. [%s]\n", params->filename, strerror(errno));
  if (gk_read(fd, &gnrows, sizeof(int)) != sizeof(int))
    errexit("Failed to read the nrows from file %s!\n", params->filename);
  if (gk_read(fd, &gncols, sizeof(int)) != sizeof(int))
    errexit("Failed to read the ncols from file %s!\n", params->filename);

  fpos = 2*sizeof(int) + gnrows*sizeof(ssize_t);
  if (lseek64(fd, fpos, SEEK_SET) == -1)
    errexit("Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
  if (gk_read(fd, &gnnz, sizeof(size_t)) != sizeof(size_t))
    errexit("Failed to read the gnnz from file %s!\n", params->filename);

  dmat->gnrows = gnrows;
  dmat->gncols = gncols;
  dmat->gnnz   = gnnz;

  lnrows = dmat->rowdist[mype+1]-dmat->rowdist[mype];

  dmat->mat = gk_csr_Create();
  dmat->mat->nrows = lnrows;
  dmat->mat->ncols = gncols;

  /* read the rowptr */
  dmat->mat->rowptr = gk_zmalloc(lnrows+1, "dmat->mat->rowptr");
  fpos = 2*sizeof(int) + dmat->rowdist[mype]*sizeof(ssize_t);
  if (lseek64(fd, fpos, SEEK_SET) == -1)
    gk_errexit(SIGERR, "Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
  if (gk_read(fd, dmat->mat->rowptr, sizeof(ssize_t)*(lnrows+1)) != sizeof(ssize_t)*(lnrows+1))
    gk_errexit(SIGERR, "Failed to read the rowptr from file %s!\n", params->filename);

  /* read the rowind */
  lnnz = dmat->mat->rowptr[lnrows]-dmat->mat->rowptr[0];
  dmat->mat->rowind = gk_imalloc(lnnz, "dmat->mat->rowind");
  fpos = 2*sizeof(int) + sizeof(ssize_t)*(gnrows+1) + sizeof(int)*dmat->mat->rowptr[0];
  if (lseek64(fd, fpos, SEEK_SET) == -1)
    gk_errexit(SIGERR, "Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
  if ((rsize = gk_read(fd, dmat->mat->rowind, sizeof(int)*lnnz)) != (ssize_t)(sizeof(int)*lnnz))
    gk_errexit(SIGERR, "Failed to read the rowind from file %s [%zd %zd]!\n", 
        params->filename, rsize, sizeof(float)*lnnz);

  /* read the rowval */
  lnnz = dmat->mat->rowptr[lnrows]-dmat->mat->rowptr[0];
  dmat->mat->rowval = gk_fmalloc(lnnz, "dmat->mat->rowval");
  fpos = 2*sizeof(int) + sizeof(ssize_t)*(gnrows+1) + sizeof(int)*gnnz + sizeof(float)*dmat->mat->rowptr[0];
  if (lseek64(fd, fpos, SEEK_SET) == -1)
    gk_errexit(SIGERR, "Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
  if ((rsize = gk_read(fd, dmat->mat->rowval, sizeof(float)*lnnz)) != (ssize_t)(sizeof(float)*lnnz))
    gk_errexit(SIGERR, "Failed to read the rowval from file %s [%zd %zd]!\n", 
        params->filename, rsize, sizeof(float)*lnnz);

  /* localize adjust rowptr */
  for (i=lnrows; i>0; i--)
    dmat->mat->rowptr[i] -= dmat->mat->rowptr[0];
  dmat->mat->rowptr[0] = 0;

  close(fd);

  if (lrank != lsize-1)
    BDMPI_Send(&token, 1, BDMPI_INT, mype+1, 1, params->comm);

  printf("[%3d] dmat->gnrows/lnrows: %d/%d, dmat->gncols/lncols: %d/%d, "
      "dmat->gnnz/lnnz: %zu/%zu [ts: %d]\n", 
      mype, 
      dmat->gnrows, dmat->mat->nrows, 
      dmat->gncols, dmat->mat->ncols, 
      dmat->gnnz, dmat->mat->rowptr[dmat->mat->nrows],(int)time(NULL));

  /* dmat->mat will not be needed anymore, it can be saved to disk, after
     extracting some key info that is needed about the data */
  dmat->lnrows  = dmat->mat->nrows;
  dmat->lncols  = dmat->mat->ncols;
  dmat->lnnz    = dmat->mat->rowptr[dmat->lnrows];
  dmat->rvalsum = gk_fsum(dmat->lnnz, dmat->mat->rowval, 1);

  gk_csr_Write(dmat->mat, params->mfile, GK_CSR_FMT_BINROW, 1, 0);
  gk_csr_Free(&(dmat->mat));

  return dmat;
}


/**************************************************************************/
/*! This function deallocates all the memory that was used */
/**************************************************************************/
void CleanupData(params_t *params, dcsr_t *dmat)
{

  gk_csr_Free(&(dmat->mat));
  gk_free((void **)&dmat->rowdist, &dmat, LTERM);

  unlink(params->mfile);
  unlink(params->ufile);
  unlink(params->vfile);
  unlink(params->rbfile);
  unlink(params->cbfile);

  return;
}


/**************************************************************************/
/*! This function computes the page-rank scores using a push approach */
/**************************************************************************/
mf_t *ComputeFactors(params_t *params, dcsr_t *dmat)
{
  int npes=params->npes, mype=params->mype;
  size_t i, j, iter, block, nelmnts;
  int ir, ic, k, nrows, ncols, nfactors, ldu, ldv;
  int *rowind;
  ssize_t *rowptr;
  float *rowval;
  mf_t *model;
  double *u, *v, *rb, *cb, uir[params->nfactors];
  double lrnrate, lambda, alpha, mu, r, mae, rmse, gmae, grmse;
  BDMPI_Status status;

  model = (mf_t *)gk_malloc(sizeof(mf_t), "model");
  memset(model, 0, sizeof(mf_t));

  nfactors = params->nfactors;
  lrnrate  = params->lrnrate;
  alpha    = params->alpha;
  lambda   = params->lambda;

  nrows  = dmat->lnrows;
  ncols  = dmat->lncols;

  /* compute model->mu */
  BDMPI_Allreduce(&dmat->rvalsum, &model->mu, 1, BDMPI_DOUBLE, BDMPI_SUM, params->comm);
  mu = model->mu = model->mu/dmat->gnnz;

  /* setup ldu/ldv */
  ldu = ldv = nfactors;

  gk_isrand(101+params->mype);

    
  /* get into the sgd iterations */
  rb = cb = u = v = NULL; 
  for (iter=0; iter<params->niters; iter++) {
    if (mype == 0)
      printf("Working on iteration: %zu\n", iter);

    for (block=0; block<npes; block++) {
      /* if active at this iteration */
      if (mype == block) {
        printf("[%3d] Computing during %zu.%zu [ts: %d] SP.\n", 
            mype, iter, block, (int)time(NULL));

        /* load the data */
        dmat->mat = gk_csr_Read(params->mfile, GK_CSR_FMT_BINROW, 1, 0);
        rowptr = dmat->mat->rowptr;
        rowind = dmat->mat->rowind;
        rowval = dmat->mat->rowval;

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
        if (model->rb == NULL) {
          rb = model->rb = gk_dreadfilebin(params->rbfile, &nelmnts);
          u  = model->u  = gk_dreadfilebin(params->ufile, &nelmnts);
          unlink(params->rbfile);
          unlink(params->ufile);
        }

        /* if first time, allocate and initialize */
        if (iter == 0 && model->cb == NULL) { 
          cb = model->cb = gk_dmalloc(ncols, "cb");
          v  = model->v  = gk_dmalloc(ncols*nfactors, "v");
          for (ic=0; ic<ncols; ic++) {
            cb[ic] = (RandomInRange(20000)-10000)*0.000001;
            for (k=0; k<nfactors; k++)
              v[ic*ldv+k] = (RandomInRange(20000)-10000)*0.000001;
          }
        }
        /*
        if (model->cb == NULL) {
          cb = model->cb = gk_dreadfilebin(params->cbfile, &nelmnts);
          v  = model->v  = gk_dreadfilebin(params->vfile, &nelmnts);
          unlink(params->cbfile);
          unlink(params->vfile);
        }
        */

        printf("[%3d] Computing during %zu.%zu [ts: %d] P1.\n", 
            mype, iter, block, (int)time(NULL));

        /* perform SGD with random sampling with replacement */
        mae = rmse = 0.0;
        for (i=0; i<nrows; i++) {
          ir = RandomInRange(nrows);

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

        printf("[%3d] Computing during %zu.%zu [ts: %d] P2.\n", 
            mype, iter, block, (int)time(NULL));

        /* free the read-only data */
        gk_csr_Free(&(dmat->mat));

        /* store row model to disk */
        gk_dwritefilebin(params->rbfile, nrows, rb);
        gk_dwritefilebin(params->ufile, nrows*nfactors, u);
        gk_free((void **)&model->rb, &model->u, LTERM);
        rb = u = NULL;

        /* send col model data to down */
        //printf("[%3d] Sending col to [%3d]\n", mype, (mype+1)%npes);
        BDMPI_Send(cb, ncols, BDMPI_DOUBLE, (block+1)%npes, 1, params->comm);
        BDMPI_Send(v, ncols*nfactors, BDMPI_DOUBLE, (block+1)%npes, 1, params->comm);

        gk_free((void **)&model->cb, &model->v, LTERM);
        cb = v = NULL;

        printf("[%3d] Computing during %zu.%zu [ts: %d] EP.\n", 
            mype, iter, block, (int)time(NULL));
      }

      /* see if you are receiving the column model data */
      if (mype == (block+1)%npes) {
        //printf("[%3d] Receiving col from [%3d]\n", mype, block);
        cb = model->cb = gk_dmalloc(ncols, "cb");
        BDMPI_Recv(cb, ncols, BDMPI_DOUBLE, block, 1, params->comm, &status);
        //gk_dwritefilebin(params->cbfile, ncols, cb);
        //gk_free((void **)&cb, LTERM);

        v = model->v = gk_dmalloc(ncols*nfactors, "v");
        BDMPI_Recv(v, ncols*nfactors, BDMPI_DOUBLE, block, 1, params->comm, &status);
        //gk_dwritefilebin(params->vfile, ncols*nfactors, v);
        //gk_free((void **)&v, LTERM);
      }
    }

    //BDMPI_Barrier(params->comm);
    BDMPI_Reduce(&mae, &gmae, 1, BDMPI_DOUBLE, BDMPI_SUM, 0, params->comm);
    BDMPI_Reduce(&rmse, &grmse, 1, BDMPI_DOUBLE, BDMPI_SUM, 0, params->comm);
    if (mype == 0) {
      gmae = gmae/dmat->gnnz;
      grmse = sqrt(grmse/dmat->gnnz);
      printf("Iter %4zu: MAE: %3.5lf RMSE: %3.5lf\n", iter, gmae, grmse);
    }
  }

  return model;
}


#ifdef XXX
/**************************************************************************/
/*! Writes the page-rank vector. It just let each process write its portion
    to the file in a round-robin fashion.
*/
/**************************************************************************/
void WriteFactors(params_t *params, dcsr_t *dmat, mf_t *model)
{
  int npes=params->npes, mype=params->mype, dummy=0;
  size_t i;
  BDMPI_Status status;
  FILE *fpout;
  char outfile[1024];

  sprintf(outfile, "%s.prvec", params->filename);

  if (mype == 0) {
    fpout = gk_fopen(outfile, "w", "outfile");
    for (i=0; i<dmat->mat->nrows; i++)
      fprintf(fpout, "%.8le\n", prvec[i]);
    gk_fclose(fpout);
    if (mype+1 < npes)
      BDMPI_Send(&dummy, 1, BDMPI_INT, mype+1, 1, params->comm);
  }
  else {
    BDMPI_Recv(&dummy, 1, BDMPI_INT, mype-1, 1, params->comm, &status);
    fpout = gk_fopen(outfile, "a", "outfile");
    for (i=0; i<dmat->mat->nrows; i++)
      fprintf(fpout, "%lf\n", prvec[i]);
    gk_fclose(fpout);
    if (mype+1 < npes)
      BDMPI_Send(&mype, 1, BDMPI_INT, mype+1, 1, params->comm);
  }
}
#endif


