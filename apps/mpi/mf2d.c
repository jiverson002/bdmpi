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
#include <fcntl.h>
#include <unistd.h>


/**************************************************************************/
/* data structures */
/**************************************************************************/
typedef struct {
  int npes, mype, nrowpes, ncolpes, myrow, mycol;
  MPI_Comm comm, rowcomm, colcomm;
  char *filename;

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
  int *rowdist, *coldist;

  /* sendinfo */
  int nsend;
  int *scounts;
  int *sdispls;
  int *sinds;

  /* recvinfo */
  int nrecv;
  int *rcounts;
  int *rdispls;
  int *rinds;

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
dcsr_t *LoadData2D(params_t *params);
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
      fprintf(stderr, "Usage: %s filename nrowpes ncolpes niters nfactors\n", argv[0]);
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  params->filename = gk_strdup(argv[1]);
  params->nrowpes  = atoi(argv[2]);
  params->ncolpes  = atoi(argv[3]);
  params->niters   = atoi(argv[4]);
  params->nfactors = atoi(argv[5]);
  params->lrnrate  = 0.001;
  params->alpha    = 0.01;
  params->lambda   = 0.01;

  if (params->npes != params->nrowpes*params->ncolpes) {
    fprintf(stderr, "The number of processors specified does not match nrowpes*ncolpes.\n");
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  gk_clearwctimer(params->totalTmr);
  gk_clearwctimer(params->loadTmr);
  gk_clearwctimer(params->setupTmr);
  gk_clearwctimer(params->compTmr);
  gk_clearwctimer(params->commTmr);

  gk_startwctimer(params->totalTmr);

  gk_startwctimer(params->setupTmr);
  SetupComms(params);
  gk_stopwctimer(params->setupTmr);

  gk_startwctimer(params->loadTmr);
  dmat = LoadData2D(params);
  gk_stopwctimer(params->loadTmr);

  gk_startwctimer(params->compTmr);
  model = ComputeFactors(params, dmat);
  gk_stopwctimer(params->compTmr);

  //WriteFactors(params, dmat, model);

  CleanupData(params, dmat);

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

//DONE:
  MPI_Finalize();

  return EXIT_SUCCESS;
}


/**************************************************************************/
/*! Creates the row/column communicators and PE ids within. */
/**************************************************************************/
void SetupComms(params_t *params)
{
  int nrowpes, ncolpes, myrow, mycol;

  myrow = params->mype / params->ncolpes;
  mycol = params->mype % params->ncolpes;

  MPI_Comm_split(params->comm, myrow, 1, &params->rowcomm);
  MPI_Comm_split(params->comm, mycol, 1, &params->colcomm);

  MPI_Comm_rank(params->rowcomm, &params->mycol);
  MPI_Comm_rank(params->colcomm, &params->myrow);
  MPI_Comm_size(params->rowcomm, &ncolpes);
  MPI_Comm_size(params->colcomm, &nrowpes);

  printf("[%3d] [nrowpes x ncolpes] = [%d %d]/[%d %d], [myrow, mycol] = [%d %d]/[%d %d]\n",
      params->mype, 
      params->nrowpes, params->ncolpes, nrowpes, ncolpes,
      myrow, mycol, params->myrow, params->mycol); 
}
  

/**************************************************************************/
/*! Reads a sparse matrix in binary CSR format, one process at a time. 
    \returns the local portion of the matrix.
*/
/**************************************************************************/
dcsr_t *LoadData2D(params_t *params)
{
  int mype=params->mype, npes=params->npes, 
      myrow=params->myrow, mycol=params->mycol,
      nrowpes=params->nrowpes, ncolpes=params->ncolpes,
      token=1;
  size_t i, j, gnnz, lnnz, nnz;
  int fd, gnrows, gncols, lnrows, lncols, firstcol, lastcol;
  dcsr_t *dmat=NULL;
  ssize_t *rowptr;
  int *rowind;
  float *rowval;
  MPI_Status status;
  off64_t fpos;

  if (mype == 0) {
    if (!gk_fexists(params->filename)) 
      errexit("File %s does not exist!\n", params->filename);
  }

  dmat = (dcsr_t *)gk_malloc(sizeof(dcsr_t), "dmat");
  memset(dmat, 0, sizeof(dcsr_t));

  /* wait your turn */
  if (mype != 0) 
    MPI_Recv(&token, 1, MPI_INT, mype-1, 1, params->comm, &status);

  if ((fd = open(params->filename, O_RDONLY)) == -1)
    errexit("Failed opeing the file %s. [%s]\n", params->filename, strerror(errno));
  if (read(fd, &gnrows, sizeof(int)) != sizeof(int))
    errexit("Failed to read the nrows from file %s!\n", params->filename);
  if (read(fd, &gncols, sizeof(int)) != sizeof(int))
    errexit("Failed to read the ncols from file %s!\n", params->filename);

  fpos = 2*sizeof(int) + gnrows*sizeof(ssize_t);
  if (lseek64(fd, fpos, SEEK_SET) == -1)
    errexit("Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
  if (read(fd, &gnnz, sizeof(size_t)) != sizeof(size_t))
    errexit("Failed to read the gnnz from file %s!\n", params->filename);

  dmat->gnrows = gnrows;
  dmat->gncols = gncols;
  dmat->gnnz   = gnnz;

  /* create and populate the rowdist/coldist */
  lnrows = (gnrows+nrowpes-1)/nrowpes;
  dmat->rowdist = gk_imalloc(nrowpes+1, "rowdist");
  dmat->rowdist[0] = 0;
  for (i=1; i<nrowpes; i++) 
    dmat->rowdist[i] = dmat->rowdist[i-1] + lnrows;
  dmat->rowdist[nrowpes] = gnrows;
  lnrows = dmat->rowdist[myrow+1]-dmat->rowdist[myrow];

  lncols = (gncols+ncolpes-1)/ncolpes;
  dmat->coldist = gk_imalloc(ncolpes+1, "coldist");
  dmat->coldist[0] = 0;
  for (i=1; i<ncolpes; i++) 
    dmat->coldist[i] = dmat->coldist[i-1] + lncols;
  dmat->coldist[ncolpes] = gncols;
  lncols = dmat->coldist[mycol+1]-dmat->coldist[mycol];


  dmat->mat = gk_csr_Create();
  dmat->mat->nrows = lnrows;
  dmat->mat->ncols = lncols;

  /* read the rowptr */
  rowptr = gk_zmalloc(lnrows+1, "rowptr");
  fpos = 2*sizeof(int) + dmat->rowdist[myrow]*sizeof(ssize_t);
  if (lseek64(fd, fpos, SEEK_SET) == -1)
    gk_errexit(SIGERR, "Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
  if (read(fd, rowptr, sizeof(ssize_t)*(lnrows+1)) != sizeof(ssize_t)*(lnrows+1))
    gk_errexit(SIGERR, "Failed to read the rowptr from file %s!\n", params->filename);

  /* read the rowind */
  nnz = rowptr[lnrows]-rowptr[0];
  rowind = gk_imalloc(nnz, "rowind");
  fpos = 2*sizeof(int) + sizeof(ssize_t)*(gnrows+1) + sizeof(int)*rowptr[0];
  if (lseek64(fd, fpos, SEEK_SET) == -1)
    gk_errexit(SIGERR, "Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
  if (read(fd, rowind, sizeof(int)*nnz) != sizeof(int)*nnz)
    gk_errexit(SIGERR, "Failed to read the rowind from file %s!\n", params->filename);

  /* read the rowval */
  nnz = rowptr[lnrows]-rowptr[0];
  rowval = gk_fmalloc(nnz, "rowval");
  fpos = 2*sizeof(int) + sizeof(ssize_t)*(gnrows+1) + sizeof(int)*gnnz + sizeof(float)*rowptr[0];
  if (lseek64(fd, fpos, SEEK_SET) == -1)
    gk_errexit(SIGERR, "Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
  if (read(fd, rowval, sizeof(float)*nnz) != sizeof(float)*nnz)
    gk_errexit(SIGERR, "Failed to read the rowval from file %s!\n", params->filename);

  /* localize adjust rowptr */
  for (i=lnrows; i>0; i--)
    rowptr[i] -= rowptr[0];
  rowptr[0] = 0;

  close(fd);

  /* extract the submatrix that you need from rowptr/rowind/rowval */
  firstcol = dmat->coldist[mycol];
  lastcol  = dmat->coldist[mycol+1];
  for (lnnz=0, i=0; i<lnrows; i++) {
    for (j=rowptr[i]; j<rowptr[i+1]; j++) {
      if (rowind[j] >= firstcol && rowind[j] < lastcol)
        lnnz++;
    }
  }
  dmat->mat->rowptr = rowptr;
  dmat->mat->rowind = gk_imalloc(lnnz, "dmat->mat->rowind");
  dmat->mat->rowval = gk_fmalloc(lnnz, "dmat->mat->rowval");

  for (lnnz=0, i=0; i<lnrows; i++) {
    for (j=rowptr[i]; j<rowptr[i+1]; j++) {
      if (rowind[j] >= firstcol && rowind[j] < lastcol) {
        dmat->mat->rowind[lnnz] = rowind[j]-firstcol;
        dmat->mat->rowval[lnnz] = rowval[j];
        lnnz++;
      }
    }
    rowptr[i] = lnnz;
  }
  SHIFTCSR(i, lnrows, rowptr);

  gk_free((void **)&rowind, &rowval, LTERM);


  /* tell the next processor to have a go */
  if (mype != npes-1)
    MPI_Send(&token, 1, MPI_INT, mype+1, 1, params->comm);

  printf("[%3d, %3d] dmat->gnrows/lnrows: %d/%d, dmat->gncols/lncols: %d/%d, dmat->gnnz/lnnz: %zu/%zu\n", 
      myrow, mycol, 
      dmat->gnrows, dmat->mat->nrows, 
      dmat->gncols, dmat->mat->ncols, 
      dmat->gnnz, dmat->mat->rowptr[dmat->mat->nrows]);

  return dmat;
}


/**************************************************************************/
/*! This function deallocates all the memory that was used */
/**************************************************************************/
void CleanupData(params_t *params, dcsr_t *dmat)
{

  MPI_Comm_free(&params->rowcomm);
  MPI_Comm_free(&params->colcomm);

  gk_csr_Free(&(dmat->mat));
  gk_free((void **)&dmat->rowdist, &dmat->coldist,
                   &dmat->scounts, &dmat->sdispls, 
                   &dmat->rcounts, &dmat->rdispls, &dmat->rinds, 
                   &dmat, LTERM);

  return;
}


/**************************************************************************/
/*! This function computes the page-rank scores using a push approach */
/**************************************************************************/
mf_t *ComputeFactors(params_t *params, dcsr_t *dmat)
{
  int nrowpes=params->nrowpes, ncolpes=params->ncolpes,
      myrow=params->myrow, mycol=params->mycol;
  size_t i, j, iter, block, nsamples;
  int ir, ic, k, nrows, ncols, nfactors, ldu, ldv, nlarge;
  int *rowind;
  ssize_t *rowptr;
  float *rowval;
  mf_t *model;
  double *u, *v, *rb, *cb, uir[params->nfactors];
  double lrnrate, lambda, alpha, mu, r, mae, rmse, gmae, grmse;
  MPI_Status status;

  model = (mf_t *)gk_malloc(sizeof(mf_t), "model");
  memset(model, 0, sizeof(mf_t));

  nlarge = gk_max(nrowpes, ncolpes);

  nfactors = params->nfactors;
  lrnrate  = params->lrnrate;
  alpha    = params->alpha;
  lambda   = params->lambda;
  nsamples = dmat->gnnz/params->npes;

  nrows  = dmat->mat->nrows;
  ncols  = dmat->mat->ncols;
  rowptr = dmat->mat->rowptr;
  rowind = dmat->mat->rowind;
  rowval = dmat->mat->rowval;

  /* compute model->mu */
  mu = gk_fsum(rowptr[nrows], rowval, 1);
  MPI_Allreduce(&mu, &model->mu, 1, MPI_DOUBLE, MPI_SUM, params->comm);
  mu = model->mu = model->mu/dmat->gnnz;

  /* allocate memory for the model */
  rb = model->rb = gk_dmalloc(nrows, "rb");
  cb = model->cb = gk_dmalloc(ncols, "cb");
  u  = model->u  = gk_dmalloc(nrows*nfactors, "u");
  v  = model->v  = gk_dmalloc(ncols*nfactors, "v");

  /* setup ldu/ldv */
  ldu = ldv = nfactors;

  /* initialize the row/col models at the processors that are hit first */
  gk_isrand(101+params->mype);
  if (mycol == (myrow >= ncolpes ? 0 : myrow)) {
    //printf("[%3d %3d] Initializing row model\n", myrow, mycol);
    for (ir=0; ir<nrows; ir++) {
      rb[ir] = (RandomInRange(20000)-10000)*0.000001;
      for (k=0; k<nfactors; k++)
        u[ir*ldu+k] = (RandomInRange(20000)-10000)*0.000001;
    }
  }
  if (myrow == (mycol >= nrowpes ? nrowpes-1 : mycol)) {
    //printf("[%3d %3d] Initializing col model\n", myrow, mycol);
    for (ic=0; ic<ncols; ic++) {
      cb[ic] = (RandomInRange(20000)-10000)*0.000001;
      for (k=0; k<nfactors; k++)
        v[ic*ldv+k] = (RandomInRange(20000)-10000)*0.000001;
    }
  }

  //MPI_Barrier(params->comm);
    
  /* get into the sgd iterations */
  for (iter=0; iter<params->niters; iter++) {
    for (block=0; block<nlarge; block++) {
      if (mycol == (myrow+block)%nlarge) { /* active at this iteration */
        //printf("[%3d %3d] Computing during %zu.%zu\n", myrow, mycol, iter, block);

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

          /* update the factors */
          rb[ir] += lrnrate*(r - alpha*rb[ir]);
          cb[ic] += lrnrate*(r - alpha*cb[ic]);

          for (k=0; k<nfactors; k++)
            u[ir*ldu+k] += lrnrate*(r*v[ic*ldv+k] - lambda*uir[k]);

          for (k=0; k<nfactors; k++)
            v[ic*ldv+k] += lrnrate*(r*uir[k] - lambda*v[ic*ldv+k]);
        }


        /* send row model data to right */
        //printf("[%3d %3d] Sending row to [%3d %3d]\n", myrow, mycol, myrow, (mycol+1)%ncolpes);
        MPI_Send(rb, nrows, MPI_DOUBLE, (mycol+1)%ncolpes, 1, params->rowcomm);
        MPI_Send(u, nrows*nfactors, MPI_DOUBLE, (mycol+1)%ncolpes, 1, params->rowcomm);

        /* send col model data to up */
        //printf("[%3d %3d] Sending col to [%3d %3d]\n", myrow, mycol, (myrow+nrowpes-1)%nrowpes, mycol);
        MPI_Send(cb, ncols, MPI_DOUBLE, (myrow+nrowpes-1)%nrowpes, 1, params->colcomm);
        MPI_Send(v, ncols*nfactors, MPI_DOUBLE, (myrow+nrowpes-1)%nrowpes, 1, params->colcomm);
      }

      /* see if you are receiving the row model data */
      if ((mycol+ncolpes-1)%ncolpes == (myrow+block)%nlarge) { 
        //printf("[%3d %3d] Receiving row from [%3d %3d]\n", myrow, mycol, myrow, (mycol+ncolpes-1)%ncolpes);
        MPI_Recv(rb, nrows, MPI_DOUBLE, (mycol+ncolpes-1)%ncolpes, 1, params->rowcomm, &status);
        MPI_Recv(u, nrows*nfactors, MPI_DOUBLE, (mycol+ncolpes-1)%ncolpes, 1, params->rowcomm, &status);
      }

      /* see if you are receiving the column model data */
      if (mycol == ((myrow+1)%nrowpes+block)%nlarge) { 
        //printf("[%3d %3d] Receiving col from [%3d %3d]\n", myrow, mycol, (myrow+1)%nrowpes, mycol);
        MPI_Recv(cb, ncols, MPI_DOUBLE, (myrow+1)%nrowpes, 1, params->colcomm, &status);
        MPI_Recv(v, ncols*nfactors, MPI_DOUBLE, (myrow+1)%nrowpes, 1, params->colcomm, &status);
      }

      //MPI_Barrier(params->comm);
    }

    /* compute a training RMSE */
    MPI_Bcast(rb, nrows, MPI_DOUBLE, (myrow < ncolpes ? myrow : 0), params->rowcomm);
    MPI_Bcast(u, nrows*nfactors, MPI_DOUBLE, (myrow < ncolpes ? myrow : 0), params->rowcomm);

    MPI_Bcast(cb, ncols, MPI_DOUBLE, (mycol < nrowpes ? mycol : nrowpes-1), params->colcomm);
    MPI_Bcast(v, ncols*nfactors, MPI_DOUBLE, (mycol < nrowpes ? mycol : nrowpes-1), params->colcomm);
    
    mae = rmse = 0.0;
    for (ir=0; ir<nrows; ir++) {
      for (j=rowptr[ir]; j<rowptr[ir+1]; j++) {
        ic = rowind[j];
        r = rowval[j] - mu - rb[ir] - cb[ic];
        for (k=0; k<nfactors; k++)
          r -= u[ir*ldu+k]*v[ic*ldv+k];
        mae += fabs(r);
        rmse += r*r;
      }
    }
    MPI_Reduce(&mae, &gmae, 1, MPI_DOUBLE, MPI_SUM, 0, params->comm);
    MPI_Reduce(&rmse, &grmse, 1, MPI_DOUBLE, MPI_SUM, 0, params->comm);
    if (params->mype == 0) {
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
  MPI_Status status;
  FILE *fpout;
  char outfile[1024];

  sprintf(outfile, "%s.prvec", params->filename);

  if (mype == 0) {
    fpout = gk_fopen(outfile, "w", "outfile");
    for (i=0; i<dmat->mat->nrows; i++)
      fprintf(fpout, "%.8le\n", prvec[i]);
    gk_fclose(fpout);
    if (mype+1 < npes)
      MPI_Send(&dummy, 1, MPI_INT, mype+1, 1, params->comm);
  }
  else {
    MPI_Recv(&dummy, 1, MPI_INT, mype-1, 1, params->comm, &status);
    fpout = gk_fopen(outfile, "a", "outfile");
    for (i=0; i<dmat->mat->nrows; i++)
      fprintf(fpout, "%lf\n", prvec[i]);
    gk_fclose(fpout);
    if (mype+1 < npes)
      MPI_Send(&mype, 1, MPI_INT, mype+1, 1, params->comm);
  }
}
#endif


