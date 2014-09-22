/*!
\file
\brief A parallel spherical k-means program
\date Started 4/20/2013
\author George
*/


#include <GKlib.h>
#include <bdmpi.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>


/**************************************************************************/
/* data structures */
/**************************************************************************/
typedef struct {
  int npes, mype, nthreads;
  BDMPI_Comm comm;
  int nclusters;
  int ntrials, niters;
  char *filename;

  /* the total number of rows and their overall distribution */
  int nrows, ncols;

  /* temporary files */
  char mfile[8192];
  char cfile[8192];
  char cpfile[8192];
  char bpfile[8192];

  /* timers */
  double totalTmr;
  double compTmr;
  double commTmr;
} params_t;


/**************************************************************************/
/* prototypes */
/**************************************************************************/
int LoadData(params_t *params);
void WriteClustering(params_t *params);
void PreprocessData(params_t *params);
void ClusterData(params_t *params, int nrows);
void ComputeClusteringStatistics(params_t *params);


/**************************************************************************/
/**************************************************************************/
int main(int argc, char **argv)
{
  params_t *params;
  int nrows;
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

  if (argc != 6) {
    if (params->mype == 0)
      fprintf(stderr, "Usage: %s filename #clusters #trials #iters #threads\n", argv[0]);

    BDMPI_Finalize();
    return EXIT_FAILURE;
  }

  params->filename  = strdup(argv[1]);
  params->nclusters = atoi(argv[2]);
  params->ntrials   = atoi(argv[3]);
  params->niters    = atoi(argv[4]);
  params->nthreads  = atoi(argv[5]);
  sprintf(params->mfile, "m%d.%d", (int)getpid(), params->mype);
  sprintf(params->cfile, "c%d.%d", (int)getpid(), params->mype);
  sprintf(params->cpfile, "cp%d.%d", (int)getpid(), params->mype);
  sprintf(params->bpfile, "bp%d.%d", (int)getpid(), params->mype);

  omp_set_num_threads(params->nthreads);

  gk_clearwctimer(params->totalTmr);
  gk_clearwctimer(params->compTmr);
  gk_clearwctimer(params->commTmr);

  BDMPI_Barrier(params->comm);
  BDMPI_Barrier(params->comm);
  gk_startwctimer(params->totalTmr);

  nrows = LoadData(params);

  PreprocessData(params);

  srand(params->mype+101);

  printf("[%03d] timestamp01: %zu\n", params->mype, (size_t)time(NULL));
  ClusterData(params, nrows);
  printf("[%03d] timestamp02: %zu\n", params->mype, (size_t)time(NULL));

  WriteClustering(params);

  BDMPI_Barrier(params->comm);
  BDMPI_Barrier(params->comm);
  gk_stopwctimer(params->totalTmr);

  /* print timing stats */
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


  BDMPI_Finalize();

  return EXIT_SUCCESS;
}


/**************************************************************************/
/*! Reads a sparse matrix in headerless CSR format. The same matrix is read
    by everybody in order to simulate a large file clustering.
    \returns the local portion of the matrix.
*/
/**************************************************************************/
int LoadData0(params_t *params)
{
  int mype = params->mype, npes = params->npes;
  size_t i, k, l, p, lnlen;
  int ncols, lnrows;
  size_t nrows, nnz, lnnz, cnnz;
  ssize_t *rowptr;
  int *rowind, ival;
  float *rowval=NULL, fval;
  char *line=NULL, *head, *tail;
  FILE *fpin;
  gk_csr_t *mat=NULL;
  BDMPI_Status status;
  long int fpos;


  mat = gk_csr_Create();

  if (mype == 0) {
    if (!gk_fexists(params->filename)) 
      gk_errexit(SIGERR, "File %s does not exist!\n", params->filename);

    gk_getfilestats(params->filename, &nrows, &nnz, NULL, NULL);
    if (nnz%2 == 1)
      gk_errexit(SIGERR, "LoadData: The number of non-zeros [%zu] is not an even number.\n", nnz);

    params->nrows = nrows;
    nnz           = nnz/2;  /* account for the values */
  }

  /* send size info to everybody */
  BDMPI_Bcast(&params->nrows, 1, BDMPI_INT, 0, params->comm);
  BDMPI_Bcast(&nnz, sizeof(size_t), BDMPI_BYTE, 0, params->comm);

  /* wait your turn */
  if (mype != 0) 
    BDMPI_Recv(&fpos, sizeof(long int), BDMPI_BYTE, mype-1, 1, params->comm, &status);
  else
    fpos = 0;

  lnnz   = nnz;
  lnrows = params->nrows;
  if (mype == 0)
    printf("[%3d] lnrows: %d, lnnz: %zu\n", mype, lnrows, lnnz);

  rowptr = gk_zmalloc(lnrows+1, "LoadData: rowptr");
  rowind = gk_imalloc(lnnz, "LoadData: rowind");
  rowval = gk_fsmalloc(lnnz, 1.0, "LoadData: rowval");

  params->nrows *= npes;  /* duplicate the data to each pe */
  fpos = 0; /* in order to read everything */

  fpin = gk_fopen(params->filename, "r", "LoadData: fpin");
  fseek(fpin, fpos, SEEK_SET);

  /* read the file, one row at a time */
  cnnz = rowptr[0] = ncols = 0;
  for (i=0; i<lnrows; i++) {
    do {
      if (gk_getline(&line, &lnlen, fpin) == -1)
        gk_errexit(SIGERR, "Premature end of input file: file while reading row %zd\n", i);
    } while (line[0] == '%');
  
    head = line;
    tail = NULL;
  
    /* parse the row */
    while (1) {
      ival = (int)strtol(head, &tail, 0);
      if (tail == head) 
        break;
      head = tail;
      
      fval = strtof(head, &tail);
      if (tail == head)
        gk_errexit(SIGERR, "Value could not be found for column! Row:%zd, NNZ:%zd\n", i, cnnz);
      head = tail;

      if (cnnz == lnnz) { /* adjust the memory */
        lnnz = 1.2*lnnz;
        rowind = gk_irealloc(rowind, lnnz, "LoadData: rowind");
        rowval = gk_frealloc(rowval, lnnz, "LoadData: rowval");
      }

      rowind[cnnz] = ival;
      rowval[cnnz] = fval;
      cnnz++;

      ncols = (ncols < ival ? ival : ncols);
    }
    rowptr[i+1] = cnnz;
  }

  mat->nrows  = lnrows;
  mat->ncols  = ncols+1;
  mat->rowptr = rowptr;
  mat->rowind = rowind;
  mat->rowval = rowval;

  fpos = ftell(fpin);
  gk_fclose(fpin);
  gk_free((void **)&line, LTERM);

  gk_csr_Write(mat, params->mfile, GK_CSR_FMT_BINROW, 1, 0);
  gk_csr_Free(&mat);

  if (mype != npes-1)
    BDMPI_Send(&fpos, sizeof(long int), BDMPI_BYTE, mype+1, 1, params->comm);

  BDMPI_Allreduce(&ncols, &ncols, 1, BDMPI_INT, BDMPI_MAX, params->comm);
  params->ncols = ncols+1;
  
  if (mype == 0)
    printf("[%3d] params->nrows: %d, params->ncols: %d, cnnz: %zu\n", 
        mype, params->nrows, params->ncols, cnnz);

  return lnrows;
}


/**************************************************************************/
/*! Reads a sparse matrix in binary CSR format. The same matrix is read
    by everybody in order to simulate a large file clustering.
    \returns the local portion of the matrix.
*/
/**************************************************************************/
int LoadData(params_t *params)
{
  int mype=params->mype, npes=params->npes;
  int lrank, lsize;
  int lnrows, flag;
  size_t lnnz;
  gk_csr_t *mat=NULL;
  BDMPI_Status status;

  BDMPI_Comm_lrank(params->comm, &lrank);
  BDMPI_Comm_lsize(params->comm, &lsize);

  if (mype == 0) {
    if (!gk_fexists(params->filename)) 
      gk_errexit(SIGERR, "File %s does not exist!\n", params->filename);
  }

  /* wait your turn */
  if (lrank != 0) 
    BDMPI_Recv(&flag, 1, BDMPI_INT, mype-1, 1, params->comm, &status);

  mat = gk_csr_Read(params->filename, GK_CSR_FMT_BINROW, 1, 0);
  lnrows = mat->nrows;
  lnnz = mat->rowptr[mat->nrows];

  params->nrows = lnrows*npes;
  params->ncols = mat->ncols;

  if (mype == 0)
    printf("[%3d] lnrows: %d, lnnz: %zu\n", mype, mat->nrows, lnnz);

  BDMPI_Entercritical();
  gk_csr_Write(mat, params->mfile, GK_CSR_FMT_BINROW, 1, 0);
  BDMPI_Exitcritical();
  gk_csr_Free(&mat);

  if (lrank != lsize-1)
    BDMPI_Send(&flag, 1, BDMPI_INT, mype+1, 1, params->comm);

  if (mype == 0)
    printf("[%3d] params->nrows: %d, params->ncols: %d, tnnz: %zu\n", 
        mype, params->nrows, params->ncols, npes*lnnz);

  return lnrows;
}


/**************************************************************************/
/*! Writes a clustering vector. It just let each process write its portion
    to the file in a round-robin fashion.
*/
/**************************************************************************/
void WriteClustering(params_t *params)
{
  int npes=params->npes, mype=params->mype, dummy=0, *cvec;
  size_t i, nrows;
  BDMPI_Status status;
  FILE *fpout;
  char outfile[1024];

  sprintf(outfile, "%s.part.%d", params->filename, params->nclusters);

  if (mype == 0) {
    BDMPI_Entercritical();
    cvec = gk_i32readfilebin(params->bpfile, &nrows);
    BDMPI_Exitcritical();
    unlink(params->bpfile);

    fpout = gk_fopen(outfile, "w", "outfile");
    for (i=0; i<nrows; i++)
      fprintf(fpout, "%d\n", cvec[i]);
    gk_fclose(fpout);

    gk_free((void **)&cvec, LTERM);

    if (mype+1 < npes)
      BDMPI_Send(&dummy, 1, BDMPI_INT, mype+1, 1, params->comm);
  }
  else {
    BDMPI_Recv(&dummy, 1, BDMPI_INT, mype-1, 1, params->comm, &status);

    BDMPI_Entercritical();
    cvec = gk_i32readfilebin(params->bpfile, &nrows);
    BDMPI_Exitcritical();
    unlink(params->bpfile);

    fpout = gk_fopen(outfile, "a", "outfile");
    for (i=0; i<nrows; i++)
      fprintf(fpout, "%d\n", cvec[i]);
    gk_fclose(fpout);

    gk_free((void **)&cvec, LTERM);

    if (mype+1 < npes)
      BDMPI_Send(&mype, 1, BDMPI_INT, mype+1, 1, params->comm);
  }
}


/**************************************************************************/
/*! This function performs various pre-processing steps on the matrix. 
*/
/**************************************************************************/
void PreprocessData(params_t *params)
{
  int mype=params->mype;
  gk_csr_t *mat;

  BDMPI_Entercritical();
  mat = gk_csr_Read(params->mfile, GK_CSR_FMT_BINROW, 1, 0);
  BDMPI_Exitcritical();

  gk_csr_Normalize(mat, GK_CSR_ROW, 2);

  BDMPI_Entercritical();
  gk_csr_Write(mat, params->mfile, GK_CSR_FMT_BINROW, 1, 0);
  BDMPI_Exitcritical();
  gk_csr_Free(&mat);
}


/**************************************************************************/
/*! This function computes the k-way clustering solution.
*/
/**************************************************************************/
void ClusterData(params_t *params, int nrows)
{
  int npes=params->npes, mype=params->mype;
  size_t trial, iter, i, j, k, offset, nelmnts;
  int ncols, nclusters, *rowind, *cpart, *bcpart;
  int myrnum, grnum, p;
  int lnmoves, gnmoves;
  ssize_t *rowptr;
  float *rowval, *centers;
  float dnorms[params->nclusters], crval, bcrval;
  gk_csr_t *mat;
  BDMPI_Status status;

  nclusters = params->nclusters;
  ncols     = params->ncols;

  /* perform a number of random trials */
  for (bcrval=0.0, trial=0; trial<params->ntrials; trial++) {
    /* select the initial cluster seeds */
    myrnum = rand();

    BDMPI_Allreduce(&myrnum, &grnum, 1, BDMPI_INT, BDMPI_MAX, params->comm);
    if (myrnum == grnum) 
      myrnum = mype;
    else
      myrnum = npes;

    BDMPI_Allreduce(&myrnum, &grnum, 1, BDMPI_INT, BDMPI_MIN, params->comm);

    if (grnum == mype) { /* this pe will be selecting the initial centers */
      /* load the data */
      BDMPI_Entercritical();
      mat = gk_csr_Read(params->mfile, GK_CSR_FMT_BINROW, 1, 0);
      BDMPI_Exitcritical();
      rowptr = mat->rowptr;
      rowind = mat->rowind;
      rowval = mat->rowval;

      centers = gk_fsmalloc(nclusters*ncols, 0.0, "centers");

      /* pick the centers */
      myrnum = RandomInRange(nrows/nclusters);
      for (k=0; k<nclusters; k++) {
        i = ((k+1)*myrnum)%nrows;
        for (j=rowptr[i]; j<rowptr[i+1]; j++)
          centers[rowind[j]*nclusters+k] = rowval[j];
      }
      gk_csr_Free(&mat);

      BDMPI_Entercritical();
      gk_fwritefilebin(params->cfile, nclusters*ncols, centers);
      BDMPI_Exitcritical();
      gk_free((void **)&centers, LTERM);
    }

    /* get into the iterative refinement */
    cpart = gk_ismalloc(nrows, -1, "cpart");
    BDMPI_Entercritical();
    gk_i32writefilebin(params->cpfile, nrows, cpart); 
    BDMPI_Exitcritical();
    gk_free((void **)&cpart, LTERM);

    for (iter=0; iter<params->niters; iter++) {
      if (grnum == mype) { /* root sends the centers to everybody else */
        BDMPI_Entercritical();
        centers = gk_freadfilebin(params->cfile, &nelmnts);
        BDMPI_Exitcritical();

        for (p=0; p<npes; p++) {
          if (p != grnum) 
            BDMPI_Send(centers, ncols*nclusters, BDMPI_FLOAT, p, 1, params->comm);
        }
      }
      else { /* everybody else receives the data */
        centers = gk_fmalloc(nclusters*ncols, "centers");
        BDMPI_Recv(centers, ncols*nclusters, BDMPI_FLOAT, grnum, 1, params->comm, &status);
      }

      printf("[%03d]%04zu.%04zu.0 ts: %d\n", mype, trial, iter, (int)time(NULL));

      /* assign each local row to the closest cluster */
      gk_startwctimer(params->compTmr);

      BDMPI_Entercritical();
      mat = gk_csr_Read(params->mfile, GK_CSR_FMT_BINROW, 1, 0);
      rowptr = mat->rowptr;
      rowind = mat->rowind;
      rowval = mat->rowval;

      cpart = gk_i32readfilebin(params->cpfile, &nelmnts);
      BDMPI_Exitcritical();

      lnmoves = 0;
#pragma omp parallel default(none),\
                     shared(i, nrows, nclusters, rowptr, rowind, rowval, centers, cpart),\
                     private(j, k, offset),\
                     reduction(+:lnmoves)
      {
        float sims[nclusters];

#pragma omp for schedule(dynamic,32)
        for (i=0; i<nrows; i++) {
          for (k=0; k<nclusters; k++) 
            sims[k] = 0.0;
          for (j=rowptr[i]; j<rowptr[i+1]; j++) {
            offset = rowind[j]*nclusters;
            for (k=0; k<nclusters; k++) 
              sims[k] += rowval[j]*centers[offset+k];
          }
          k = gk_fargmax(nclusters, sims, 1);

          if (k != cpart[i])
            lnmoves++;
          cpart[i] = k;
        }
      }

      /* compute the new local centers */
      gk_fset(nclusters*ncols, 0.0, centers);
      for (i=0; i<nrows; i++) {
        k = cpart[i];
        for (j=rowptr[i]; j<rowptr[i+1]; j++) 
          centers[rowind[j]*nclusters+k] += rowval[j];
      }

      BDMPI_Entercritical();
      gk_i32writefilebin(params->cpfile, nrows, cpart); 
      BDMPI_Exitcritical();
      gk_free((void **)&cpart, LTERM);

      gk_csr_Free(&mat); /* done with the matrix */

      printf("[%03d]%04zu.%04zu.1 ts: %d\n", mype, trial, iter, (int)time(NULL));

      if (mype == grnum) {
        float *rcenters;

        /* write the centers */
        BDMPI_Entercritical();
        gk_fwritefilebin(params->cfile, nclusters*ncols, centers);
        BDMPI_Exitcritical();
        gk_free((void **)&centers, LTERM);

        /* get the data from everybody else */
        for (p=0; p<npes-1; p++) {
          rcenters = gk_fmalloc(nclusters*ncols, "rcenters");
          BDMPI_Recv(rcenters, ncols*nclusters, BDMPI_FLOAT, BDMPI_ANY_SOURCE, 2, params->comm, 
               &status);

          BDMPI_Entercritical();
          centers = gk_freadfilebin(params->cfile, &nelmnts);
          BDMPI_Exitcritical();
          for (i=0; i<nclusters*ncols; i++)
            centers[i] += rcenters[i];
          BDMPI_Entercritical();
          gk_fwritefilebin(params->cfile, nclusters*ncols, centers);
          BDMPI_Exitcritical();

          gk_free((void **)&rcenters, &centers, LTERM);
        }

        /* compute the new global centers */
        for (k=0; k<nclusters; k++)
          dnorms[k] = 0.0;

        BDMPI_Entercritical();
        centers = gk_freadfilebin(params->cfile, &nelmnts);
        BDMPI_Exitcritical();
        for (i=0; i<ncols; i++) {
          offset = i*nclusters;
          for (k=0; k<nclusters; k++) 
            dnorms[k] += centers[offset+k]*centers[offset+k];
        }

        for (crval=0.0, k=0; k<nclusters; k++) {
          if (dnorms[k] > 0) {
            crval += sqrt(dnorms[k]);
            dnorms[k] = 1.0/sqrt(dnorms[k]);
          }
        }

        for (i=0; i<ncols; i++) {
          offset = i*nclusters;
          for (k=0; k<nclusters; k++) 
            centers[offset+k] *= dnorms[k];
        }
            
        //printf("trial: %2zd; iter: %3zd; crval: %.8e; bcrval: %.8e\n", trial, iter, crval, bcrval);

        BDMPI_Entercritical();
        gk_fwritefilebin(params->cfile, nclusters*ncols, centers);
        BDMPI_Exitcritical();
        gk_free((void **)&centers, LTERM);
      }
      else {
        BDMPI_Send(centers, ncols*nclusters, BDMPI_FLOAT, grnum, 2, params->comm);
        gk_free((void **)&centers, LTERM);
      }

      gk_stopwctimer(params->compTmr);

      /* see if you are done refining */
      if (iter > 0) {
        BDMPI_Allreduce(&lnmoves, &gnmoves, 1, BDMPI_INT, BDMPI_SUM, params->comm);
        //printf("[%3d] lnmoves: %d, gnmoves: %d\n", mype, lnmoves, gnmoves);
        if (gnmoves == 0)
          break;
      }

      /* send the new crval to everybody */
      BDMPI_Bcast(&crval, 1, BDMPI_FLOAT, grnum, params->comm);

      if (crval > bcrval) {
        char cmd[8192];
        sprintf(cmd, "cp %s %s", params->cpfile, params->bpfile);
        GKASSERT(system(cmd) != -1);

        bcrval = crval;
      }
    }

    if (mype == 0)
      printf("[%3zu:%3zu] gnmoves: %8d; crval: %8.4e; bcrval: %8.4e\n",
             trial, iter, gnmoves, crval, bcrval);

    if (mype == grnum)
      unlink(params->cfile);
    unlink(params->cpfile);
  }

  /* compute the clustering statistics */
  //ComputeClusteringStatistics(params);

  unlink(params->mfile);
}

  
/**************************************************************************/
/*! This function prints final statistics for the clustering solution.
*/
/**************************************************************************/
void ComputeClusteringStatistics(params_t *params)
{
  int npes=params->npes, mype=params->mype, nclusters=params->nclusters, ncols=params->ncols;
  int nrows;
  size_t i, j, k, offset, nelmnts;
  int *cpart, *rowind, *pwgts;
  ssize_t *rowptr;
  float *rowval, *centers, *dnorms;
  float crval, tcrval;
  gk_csr_t *mat;

  centers = gk_fsmalloc(nclusters*ncols, 0.0, "centers");
  pwgts   = gk_ismalloc(nclusters, 0.0, "pwgts");

  /* load the data */
  mat = gk_csr_Read(params->mfile, GK_CSR_FMT_BINROW, 1, 0);
  nrows  = mat->nrows;
  rowptr = mat->rowptr;
  rowind = mat->rowind;
  rowval = mat->rowval;

  cpart = gk_i32readfilebin(params->bpfile, &nelmnts);

  /* compute the local centers and local partition weights */
  for (i=0; i<nrows; i++) {
    k = cpart[i];
    pwgts[k]++;
    for (j=rowptr[i]; j<rowptr[i+1]; j++) {
      centers[rowind[j]*nclusters+k] += rowval[j];
    }
  }

  gk_csr_Free(&mat);
  gk_free((void **)&cpart, LTERM);

  /* compute global centroids and partition weights */
  BDMPI_Reduce(pwgts, pwgts, nclusters, BDMPI_INT, BDMPI_SUM, 0, params->comm);
  BDMPI_Reduce(centers, centers, nclusters*ncols, BDMPI_FLOAT, BDMPI_SUM, 0, params->comm);

  if (mype == 0) {
    dnorms = gk_fsmalloc(nclusters, 0.0, "dnorms");
    for (i=0; i<ncols; i++) {
      offset = i*nclusters;
      for (k=0; k<nclusters; k++) 
        dnorms[k] += centers[offset+k]*centers[offset+k];
    }

    for (tcrval=0.0, k=0; k<nclusters; k++) {
      crval = (dnorms[k] > 0 ? sqrt(dnorms[k]) : 0.0);
      tcrval += crval;
      printf("Cluster: %4zu %6d %.4e\n", k, pwgts[k], crval);
    }
    printf("Overall: %.4e\n", tcrval);

    gk_free((void **)&dnorms, LTERM);
  }

  gk_free((void **)&centers, &pwgts, LTERM);

}

