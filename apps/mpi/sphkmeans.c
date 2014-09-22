/*!
\file
\brief A parallel spherical k-means program
\date Started 4/20/2013
\author George
*/


#include <GKlib.h>
#include <mpi.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>


/**************************************************************************/
/* data structures */
/**************************************************************************/
typedef struct {
  int npes, mype;
  MPI_Comm comm;
  int nclusters;
  char *filename;
  int ntrials;
  int niters;
  int nthreads;
  int ncopies;

  /* the total number of rows and their overall distribution */
  int nrows, ncols;
  int *rowdist;

  /* timers */
  double totalTmr;
  double compTmr;
  double commTmr;
} params_t;


typedef struct {
  int val;
  int loc;
} vlp_ii_t;


/**************************************************************************/
/* prototypes */
/**************************************************************************/
gk_csr_t *LoadData(params_t *params);
void WriteClustering(params_t *params, gk_csr_t *mat, int *cvec);
void PreprocessData(params_t *params, gk_csr_t *mat);
int *ClusterData(params_t *params, gk_csr_t *mat);
float ComputeClusteringStatistics(params_t *params, gk_csr_t *mat, int *cptr, 
         int *cind);


/**************************************************************************/
/**************************************************************************/
int main(int argc, char **argv)
{
  params_t *params;
  gk_csr_t *mat;
  int *cvec;
  MPI_Status status;
  double max, current;

  MPI_Init(&argc, &argv);

  params = (params_t *)gk_malloc(sizeof(params_t), "params");
  memset(params, 0, sizeof(params_t));

  params->comm = MPI_COMM_WORLD;
  MPI_Comm_size(params->comm, &(params->npes));
  MPI_Comm_rank(params->comm, &(params->mype));

  if (argc != 7) {
    if (params->mype == 0)
      fprintf(stderr, "Usage: %s filename #copies #clusters #trials #iters #threads\n", argv[0]);

    MPI_Finalize();
    return EXIT_FAILURE;
  }

  params->filename  = strdup(argv[1]);
  params->ncopies   = atoi(argv[2]);
  params->nclusters = atoi(argv[3]);
  params->ntrials   = atoi(argv[4]);
  params->niters    = atoi(argv[5]);
  params->nthreads  = atoi(argv[6]);

  printf("[%3d] nclusters: %d, ntrials: %d, niters: %d, nthreads: %d, ncopies: %d\n",
      params->mype, params->nclusters, params->ntrials, params->niters, params->nthreads,
      params->ncopies);

  omp_set_num_threads(params->nthreads);

  gk_clearwctimer(params->totalTmr);
  gk_clearwctimer(params->compTmr);
  gk_clearwctimer(params->commTmr);

  gk_startwctimer(params->totalTmr);
  mat = LoadData(params);

  PreprocessData(params, mat);

  srand(params->mype+101);

  cvec = ClusterData(params, mat);

  //WriteClustering(params, mat, cvec);

  gk_stopwctimer(params->totalTmr);

  /* print timing stats */
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
/*! Reads a sparse matrix in headerless CSR format, distributes it among 
    the PEs, and creates the params->rowdist array to reflect the distribution.
    The read is the last process to ensure that it gets the last chunk.
    \returns the local portion of the matrix.
*/
/**************************************************************************/
gk_csr_t *LoadData0(params_t *params)
{
  int mype = params->mype, npes = params->npes;
  size_t i, k, l, p, lnlen;
  int ncols, lnrows, cnrows;
  size_t nrows, nnz, lnnz, cnnz, rnnz;
  ssize_t *rowptr;
  int *rowind, ival;
  float *rowval=NULL, fval;
  char *line=NULL, *head, *tail;
  FILE *fpin;
  gk_csr_t *mat=NULL;
  MPI_Status status;


  mat = gk_csr_Create();

  if (mype == npes-1) {
    if (!gk_fexists(params->filename)) 
      gk_errexit(SIGERR, "File %s does not exist!\n", params->filename);

    gk_getfilestats(params->filename, &nrows, &nnz, NULL, NULL);
    if (nnz%2 == 1)
      gk_errexit(SIGERR, "LoadData: The number of non-zeros [%zu] is not an even number.\n", nnz);

    fpin = gk_fopen(params->filename, "r", "LoadData: fpin");

    params->nrows = nrows;
    params->rowdist = gk_imalloc(npes+1, "rowdist");

    nnz    = nnz/2;  /* account for the values */
    lnnz   = 1.1*nnz/npes;
    lnrows = (nrows+npes-1)/npes;
    rowptr = gk_zmalloc(lnrows+1, "LoadData: rowptr");
    rowind = gk_imalloc(lnnz, "LoadData: rowind");
    rowval = gk_fsmalloc(lnnz, 1.0, "LoadData: rowval");

    printf("nrows: %zu, nnz: %zu\n", nrows, nnz);

    /* read the file, one row at a time */
    p = ncols = params->rowdist[0] = 0;
    cnrows = cnnz = rnnz = rowptr[0] = 0;
    for (i=0; i<nrows; i++) {
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
      rowptr[++cnrows] = cnnz;

      //printf("cnrows: %d, lnrows: %d\n", cnrows, lnrows);

      if (cnrows == lnrows || i == nrows-1) { /* time to send the data away */
        rnnz += cnnz;

        printf("Sending to %zu [cnnz: %zu; rnnz: %zu]...\n", p, cnnz, rnnz);

        if (p != npes-1) {
          MPI_Send(&cnrows, 1, MPI_INT, p, 1, params->comm);
          MPI_Send(&cnnz, sizeof(size_t), MPI_BYTE, p, 2, params->comm);
          MPI_Send(rowptr, sizeof(ssize_t)*(cnrows+1), MPI_BYTE, p, 3, params->comm);
          MPI_Send(rowind, cnnz, MPI_INT, p, 4, params->comm);
          MPI_Send(rowval, cnnz, MPI_FLOAT, p, 5, params->comm);

          cnrows = cnnz = rowptr[0] = 0;
          params->rowdist[++p] = i+1;
        }
        else {
          mat->nrows  = cnrows;
          mat->ncols  = ncols+1;
          mat->rowptr = rowptr;
          mat->rowind = rowind;
          mat->rowval = rowval;
        }
      }
    }

    if (rnnz != nnz)
      gk_errexit(SIGERR, "LoadData: Something wrong with the number of nonzeros in "
                         "the input file. NNZ=%zd, ActualNNZ=%zd.\n", nnz, rnnz);

    gk_fclose(fpin);
    gk_free((void **)&line, LTERM);

    /*
    for (p=0; p<npes-1; p++) {
      MPI_Send(&mat->ncols, 1, MPI_INT, p, 6, params->comm);
      MPI_Send(&params->nrows, 1, MPI_INT, p, 7, params->comm);
    }
    */
  }
  else {
    MPI_Recv(&mat->nrows, 1, MPI_INT, npes-1, 1, params->comm, &status);
    MPI_Recv(&cnnz, sizeof(size_t), MPI_BYTE, npes-1, 2, params->comm, &status);

    printf("[%3d] mat->nrows: %d, cnnz: %zu\n", mype, mat->nrows, cnnz);

    mat->rowptr = gk_zmalloc(mat->nrows+1, "rowptr");
    mat->rowind = gk_imalloc(cnnz, "rowind");
    mat->rowval = gk_fmalloc(cnnz, "rowval");

    MPI_Recv(mat->rowptr, sizeof(ssize_t)*(mat->nrows+1), MPI_BYTE, npes-1, 3, params->comm, &status);
    MPI_Recv(mat->rowind, sizeof(int)*cnnz, MPI_BYTE, npes-1, 4, params->comm, &status);
    MPI_Recv(mat->rowval, sizeof(float)*cnnz, MPI_BYTE, npes-1, 5, params->comm, &status);

    /*
    MPI_Recv(&mat->ncols, 1, MPI_INT, npes-1, 6, params->comm, &status);
    MPI_Recv(&params->nrows, 1, MPI_INT, npes-1, 7, params->comm, &status);
    printf("[%3d] mat->ncols: %d, params->nrows: %d\n", mype, mat->ncols, params->nrows);
    */

  }

  MPI_Bcast(&mat->ncols, 1, MPI_INT, npes-1, params->comm);
  MPI_Bcast(&params->nrows, 1, MPI_INT, npes-1, params->comm);
  printf("[%3d] mat->ncols: %d, params->nrows: %d\n", mype, mat->ncols, params->nrows);

  params->ncols = mat->ncols;

  return mat;
}


/**************************************************************************/
/*! Reads a sparse matrix in headerless CSR format. The mype==0 starts the
    reading and then passes the current fpos to the next processor to continue
    the reading.
    \returns the local portion of the matrix.
*/
/**************************************************************************/
gk_csr_t *LoadData1(params_t *params)
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
  MPI_Status status;
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
  MPI_Bcast(&params->nrows, 1, MPI_INT, 0, params->comm);
  MPI_Bcast(&nnz, sizeof(size_t), MPI_BYTE, 0, params->comm);

  /* wait your turn */
  if (mype != 0) 
    MPI_Recv(&fpos, sizeof(long int), MPI_BYTE, mype-1, 1, params->comm, &status);
  else
    fpos = 0;

  lnnz   = 1.1*nnz/npes;
  lnrows = (params->nrows+npes-1)/npes;
  lnrows = gk_max(0, gk_min(lnrows, params->nrows - mype*lnrows));
  printf("[%3d] lnrows: %d, lnnz: %zu\n", mype, lnrows, lnnz);

  rowptr = gk_zmalloc(lnrows+1, "LoadData: rowptr");
  rowind = gk_imalloc(lnnz, "LoadData: rowind");
  rowval = gk_fsmalloc(lnnz, 1.0, "LoadData: rowval");

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
  mat->rowptr = rowptr;
  mat->rowind = rowind;
  mat->rowval = rowval;

  fpos = ftell(fpin);
  gk_fclose(fpin);
  gk_free((void **)&line, LTERM);

  if (mype != npes-1)
    MPI_Send(&fpos, sizeof(long int), MPI_BYTE, mype+1, 1, params->comm);

  MPI_Allreduce(&ncols, &mat->ncols, 1, MPI_INT, MPI_MAX, params->comm);
  mat->ncols++;
  params->ncols = mat->ncols;
  
  printf("[%3d] params->nrows: %d, params->ncols: %d, cnnz: %zu\n", 
      mype, params->nrows, params->ncols, cnnz);

  return mat;
}


/**************************************************************************/
/*! Reads a sparse matrix in headerless CSR format. The same matrix is read
    by everybody in order to simulate a large file clustering.
    \returns the local portion of the matrix.
*/
/**************************************************************************/
gk_csr_t *LoadData2(params_t *params)
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
  MPI_Status status;
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
  MPI_Bcast(&params->nrows, 1, MPI_INT, 0, params->comm);
  MPI_Bcast(&nnz, sizeof(size_t), MPI_BYTE, 0, params->comm);

  /* wait your turn */
  if (mype != 0) 
    MPI_Recv(&fpos, sizeof(long int), MPI_BYTE, mype-1, 1, params->comm, &status);
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
  mat->rowptr = rowptr;
  mat->rowind = rowind;
  mat->rowval = rowval;

  fpos = ftell(fpin);
  gk_fclose(fpin);
  gk_free((void **)&line, LTERM);

  if (mype != npes-1)
    MPI_Send(&fpos, sizeof(long int), MPI_BYTE, mype+1, 1, params->comm);

  MPI_Allreduce(&ncols, &mat->ncols, 1, MPI_INT, MPI_MAX, params->comm);
  mat->ncols++;
  params->ncols = mat->ncols;
  
  if (mype == 0)
    printf("[%3d] params->nrows: %d, params->ncols: %d, cnnz: %zu\n", 
        mype, params->nrows, params->ncols, cnnz);

  return mat;
}


/**************************************************************************/
/*! Reads a sparse matrix in binary CSR format. The same matrix is read
    by everybody in order to simulate a large file clustering.
    \returns the local portion of the matrix.
*/
/**************************************************************************/
gk_csr_t *LoadData3(params_t *params)
{
  int mype=params->mype, npes=params->npes;
  int lnrows, flag;
  size_t lnnz;
  gk_csr_t *mat=NULL;
  MPI_Status status;

  if (mype == 0) {
    if (!gk_fexists(params->filename)) 
      gk_errexit(SIGERR, "File %s does not exist!\n", params->filename);
  }

  /* wait your turn */
  if (mype != 0) 
    MPI_Recv(&flag, 1, MPI_INT, mype-1, 1, params->comm, &status);

  mat = gk_csr_Read(params->filename, GK_CSR_FMT_BINROW, 1, 0);
  lnrows = mat->nrows;
  lnnz = mat->rowptr[mat->nrows];

  params->nrows = lnrows*npes;
  params->ncols = mat->ncols;

  if (mype == 0)
    printf("[%3d] lnrows: %d, lnnz: %zu\n", mype, mat->nrows, lnnz);

  if (mype != npes-1)
    MPI_Send(&flag, 1, MPI_INT, mype+1, 1, params->comm);

  if (mype == 0)
    printf("[%3d] params->nrows: %d, params->ncols: %d, tnnz: %zu\n", 
        mype, params->nrows, params->ncols, npes*lnnz);

  return mat;
}


/**************************************************************************/
/*! Reads a sparse matrix in binary CSR format. The same matrix is read
    by everybody in order to simulate a large file clustering.
    \returns the local portion of the matrix.
*/
/**************************************************************************/
gk_csr_t *LoadData(params_t *params)
{
  int mype=params->mype, npes=params->npes;
  int ic, lnrows;
  size_t i, lnnz;
  gk_csr_t *mat=NULL, *tmat=NULL;
  MPI_Status status;

  if (mype == 0) {
    if (!gk_fexists(params->filename)) 
      gk_errexit(SIGERR, "File %s does not exist!\n", params->filename);
  }

  mat = gk_csr_Create();

  for (ic=0; ic<params->ncopies; ic++) {
    tmat = gk_csr_Read(params->filename, GK_CSR_FMT_BINROW, 1, 0);
    lnrows = tmat->nrows;
    lnnz   = tmat->rowptr[lnrows];

    if (ic == 0) {
      mat->nrows  = tmat->nrows*params->ncopies;
      mat->ncols  = tmat->ncols;
      mat->rowptr = gk_zmalloc(mat->nrows+1, "rowptr");
      mat->rowind = gk_imalloc(lnnz*params->ncopies, "rowind");
      mat->rowval = gk_fmalloc(lnnz*params->ncopies, "rowval");
    }

    /* append the new data */
    for (i=0; i<lnrows; i++)
      tmat->rowptr[i] = tmat->rowptr[i+1]-tmat->rowptr[i];

    gk_zcopy(lnrows, tmat->rowptr, mat->rowptr+ic*lnrows);
    gk_icopy(lnnz, tmat->rowind, mat->rowind+ic*lnnz);
    gk_fcopy(lnnz, tmat->rowval, mat->rowval+ic*lnnz);

    gk_csr_Free(&tmat);
  }
  MAKECSR(i, mat->nrows, mat->rowptr);
      
  params->nrows = mat->nrows*npes;
  params->ncols = mat->ncols;

  printf("[%3d] lnrows: %d, lnnz: %zu\n", mype, mat->nrows, mat->rowptr[mat->nrows]);

  return mat;
}


/**************************************************************************/
/*! Writes a clustering vector. It just let each process write its portion
    to the file in a round-robin fashion.
*/
/**************************************************************************/
void WriteClustering(params_t *params, gk_csr_t *mat, int *cvec)
{
  int npes=params->npes, mype=params->mype, dummy=0;
  size_t i;
  MPI_Status status;
  FILE *fpout;
  char outfile[1024];

  sprintf(outfile, "%s.part.%d", params->filename, params->nclusters);

  if (mype == 0) {
    fpout = gk_fopen(outfile, "w", "outfile");
    for (i=0; i<mat->nrows; i++)
      fprintf(fpout, "%d\n", cvec[i]);
    gk_fclose(fpout);
    if (mype+1 < npes)
      MPI_Send(&dummy, 1, MPI_INT, mype+1, 1, params->comm);
  }
  else {
    MPI_Recv(&dummy, 1, MPI_INT, mype-1, 1, params->comm, &status);
    fpout = gk_fopen(outfile, "a", "outfile");
    for (i=0; i<mat->nrows; i++)
      fprintf(fpout, "%d\n", cvec[i]);
    gk_fclose(fpout);
    if (mype+1 < npes)
      MPI_Send(&mype, 1, MPI_INT, mype+1, 1, params->comm);
  }
}


/**************************************************************************/
/*! This function performs various pre-processing steps on the matrix. 
*/
/**************************************************************************/
void PreprocessData(params_t *params, gk_csr_t *mat)
{
  gk_csr_Normalize(mat, GK_CSR_ROW, 2);
}


/**************************************************************************/
/*! This function computes the k-way clustering solution.
*/
/**************************************************************************/
int *ClusterData(params_t *params, gk_csr_t *mat)
{
  int npes=params->npes, mype=params->mype;
  size_t trial, iter, i, j, k, offset;
  int nrows, ncols, nclusters, rnum, *rowind, *cpart=NULL, *bcpart=NULL, *tptr;
  int lnmoves, gnmoves;
  ssize_t *rowptr;
  float *rowval, *centers=NULL, *ncenters=NULL;
  float dnorms[params->nclusters], crval, bcrval;
  vlp_ii_t lmaxloc, gmaxloc;

  nclusters = params->nclusters;

  nrows  = mat->nrows;
  ncols  = mat->ncols;
  rowptr = mat->rowptr;
  rowind = mat->rowind;
  rowval = mat->rowval;

  bcpart  = gk_imalloc(nrows, "bcpart");

  /* perform a number of random trials */
  for (bcrval=0.0, trial=0; trial<params->ntrials; trial++) {
    /* select the initial cluster seeds */
    lmaxloc.val = rand();
    lmaxloc.loc = mype;

    MPI_Allreduce(&lmaxloc, &gmaxloc, 1, MPI_2INT, MPI_MAXLOC, params->comm);

    if (mype == gmaxloc.loc) { /* this pe will be selecting the initial centers */
      centers = gk_fsmalloc(nclusters*ncols, 0.0, "centers");

      /* pick the centers */
      rnum = RandomInRange(nrows/nclusters);
      for (k=0; k<nclusters; k++) {
        i = ((k+1)*rnum)%nrows;
        for (j=rowptr[i]; j<rowptr[i+1]; j++)
          centers[rowind[j]*nclusters+k] = rowval[j];
      }
    }

    /* get into the iterative refinement */
    cpart = gk_ismalloc(nrows, -1, "cpart");
    for (iter=0; iter<params->niters; iter++) {
      if (mype == gmaxloc.loc)
        printf("Working on trial: %zu, iter: %zu\n", trial, iter);
      else
        centers = gk_fmalloc(nclusters*ncols, "centers");

      MPI_Bcast(centers, ncols*nclusters, MPI_FLOAT, gmaxloc.loc, params->comm);
      printf("[%03d]%04zu.%04zu.0 ts: %d\n", mype, trial, iter, (int)time(NULL));

      /* assign each local row to the closest cluster */
      gk_startwctimer(params->compTmr);

      lnmoves = 0;
#pragma omp parallel default(none),\
                     shared(nrows, nclusters, rowptr, rowind, rowval, centers, cpart),\
                     private(i, j, k, offset),\
                     reduction(+:lnmoves)
      {
        float sims[nclusters];

#pragma omp for schedule(static,1)
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
      gk_stopwctimer(params->compTmr);

      printf("[%03d]%04zu.%04zu.1 ts: %d\n", mype, trial, iter, (int)time(NULL));


      /* compute the new global centers */
      if (mype == gmaxloc.loc) 
        ncenters = gk_fmalloc(nclusters*ncols, "ncenters");

      MPI_Reduce(centers, ncenters, nclusters*ncols, MPI_FLOAT, MPI_SUM, 
          gmaxloc.loc, params->comm);

      if (mype == gmaxloc.loc) {
        gk_free((void **)&centers, LTERM);
        centers = ncenters;
        ncenters = NULL;

        for (k=0; k<nclusters; k++)
          dnorms[k] = 0.0;

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
      }
      else {
        gk_free((void **)&centers, LTERM);
      }

      /* see if you are done refining */
      if (iter > 0) {
        MPI_Allreduce(&lnmoves, &gnmoves, 1, MPI_INT, MPI_SUM, params->comm);
        if (gnmoves == 0)
          break;
      }
    }

    MPI_Bcast(&crval, 1, MPI_FLOAT, gmaxloc.loc, params->comm);
    if (crval > bcrval) {
      gk_SWAP(cpart, bcpart, tptr);
      bcrval = crval;
    }

    gk_free((void **)&cpart, LTERM);

    if (mype == gmaxloc.loc)
      printf("[%3zu:%3zu] gnmoves: %8d; crval: %8.4e; bcrval: %8.4e [ts: %d]\n",
             trial, iter, gnmoves, crval, bcrval, (int)time(NULL));
  }

  gk_free((void **)&centers, LTERM);


  return bcpart;
}

  
/**************************************************************************/
/*! This function prints final statistics for the clustering solution.
*/
/**************************************************************************/
float ComputeClusteringStatistics(params_t *params, gk_csr_t *mat, int *cptr, 
         int *cind)
{
  int npes=params->npes, mype=params->mype;
  size_t i, ii, j, k;
  int ncols, ncrows, *rowind, *cpart, *bcpart, *rowids;
  int lpwgt, gpwgt;
  ssize_t *rowptr;
  float *rowval, *lvec, *gvec;
  float dnorm, crval, tcrval;

  //printf("[%3d] Bisecting cluster %d\n", mype, cnum);

  ncols  = mat->ncols;
  rowptr = mat->rowptr;
  rowind = mat->rowind;
  rowval = mat->rowval;

  lvec = gk_fmalloc(ncols, "lvec");
  gvec = (mype == 0 ? gk_fmalloc(ncols, "gvec") : NULL);

  for (tcrval=0.0, k=0; k<params->nclusters; k++) {
    /* compute the new local composite vector */
    gk_fset(ncols, 0.0, lvec);
    for (ii=cptr[k]; ii<cptr[k+1]; ii++) {
      i = cind[ii];
      for (j=rowptr[i]; j<rowptr[i+1]; j++) 
        lvec[rowind[j]] += rowval[j];
    }

    /* compute the global composite vector and cluster size */
    lpwgt = cptr[k+1]-cptr[k];
    MPI_Reduce(lvec, gvec, ncols, MPI_FLOAT, MPI_SUM, 0, params->comm);
    MPI_Reduce(&lpwgt, &gpwgt, 1, MPI_INT, MPI_SUM, 0, params->comm);

    if (mype == 0) {
      dnorm = gk_fdot(ncols, gvec, 1, gvec, 1);
      crval = (dnorm > 0 ? sqrt(dnorm) : 0.0);
      tcrval += crval;
      printf("Cluster: %4zu %6d %.4e\n", k, gpwgt, crval);
    }
  }
  if (mype == 0)
    printf("Overall: %.4e\n", tcrval);

  gk_free((void **)&lvec, &gvec, LTERM);

  return tcrval;
}

