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
#include <fcntl.h>
#include <unistd.h>


/**************************************************************************/
/* data structures */
/**************************************************************************/
typedef struct {
  int npes, mype, nrowpes, ncolpes, myrow, mycol;
  BDMPI_Comm comm, rowcomm, colcomm, diagcomm;
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
  int type;


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
dcsr_t *LoadData2D(params_t *params);
void CleanupData(params_t *params, dcsr_t *dmat);
mf_t *ComputeFactors_FR(params_t *params, dcsr_t *dmat);
mf_t *ComputeFactors_RR(params_t *params, dcsr_t *dmat);
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


  if (argc != 7) {
    if (params->mype == 0)
      fprintf(stderr, "Usage: %s filename nrowpes ncolpes niters nfactors type[1:FR;2:RR]\n", argv[0]);
    BDMPI_Finalize();
    return EXIT_FAILURE;
  }

  params->filename = gk_strdup(argv[1]);
  params->nrowpes  = atoi(argv[2]);
  params->ncolpes  = atoi(argv[3]);
  params->niters   = atoi(argv[4]);
  params->nfactors = atoi(argv[5]);
  params->type     = atoi(argv[6]);
  params->lrnrate  = 0.001;
  params->alpha    = 0.01;
  params->lambda   = 0.01;

  printf("[%3d] grid: %dx%d, niters: %d, nfactors: %d, type: %d\n",
            params->mype, params->nrowpes, params->ncolpes, params->niters,
            params->nfactors, params->type);

  sprintf(params->mfile,  "m%d.%d", (int)getpid(), params->mype);
  sprintf(params->ufile,  "u%d.%d", (int)getpid(), params->mype);
  sprintf(params->vfile,  "v%d.%d", (int)getpid(), params->mype);
  sprintf(params->rbfile, "rb%d.%d", (int)getpid(), params->mype);
  sprintf(params->cbfile, "cb%d.%d", (int)getpid(), params->mype);

  if (params->npes != params->nrowpes*params->ncolpes) {
    fprintf(stderr, "The number of processors specified does not match nrowpes*ncolpes.\n");
    BDMPI_Finalize();
    return EXIT_FAILURE;
  }

  gk_clearwctimer(params->totalTmr);
  gk_clearwctimer(params->loadTmr);
  gk_clearwctimer(params->setupTmr);
  gk_clearwctimer(params->compTmr);
  gk_clearwctimer(params->commTmr);

  BDMPI_Barrier(params->comm);
  BDMPI_Barrier(params->comm);
  gk_startwctimer(params->totalTmr);

  gk_startwctimer(params->setupTmr);
  SetupComms(params);
  gk_stopwctimer(params->setupTmr);

  gk_startwctimer(params->loadTmr);
  dmat = LoadData2D(params);
  gk_stopwctimer(params->loadTmr);

  gk_startwctimer(params->compTmr);
  if (params->type == 1)
    model = ComputeFactors_FR(params, dmat);
  else
    model = ComputeFactors_RR(params, dmat);
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
/*! Creates the row/column communicators and PE ids within. */
/**************************************************************************/
void SetupComms(params_t *params)
{
  int nrowpes, ncolpes, myrow, mycol;

  myrow = params->mype / params->ncolpes;
  mycol = params->mype % params->ncolpes;

  BDMPI_Comm_split(params->comm, myrow, 1, &params->rowcomm);
  BDMPI_Comm_split(params->comm, mycol, 1, &params->colcomm);

  BDMPI_Comm_rank(params->rowcomm, &params->mycol);
  BDMPI_Comm_rank(params->colcomm, &params->myrow);
  BDMPI_Comm_size(params->rowcomm, &ncolpes);
  BDMPI_Comm_size(params->colcomm, &nrowpes);

  BDMPI_Comm_split(params->comm, (mycol-myrow+params->ncolpes)%params->ncolpes, 
      1, &params->diagcomm);

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
dcsr_t *LoadData2D0(params_t *params)
{
  int mype=params->mype, npes=params->npes, 
      myrow=params->myrow, mycol=params->mycol,
      nrowpes=params->nrowpes, ncolpes=params->ncolpes,
      lrank, lsize, token=1;
  size_t i, j, gnnz, lnnz, nnz;
  int fd, gnrows, gncols, lnrows, lncols, firstcol, lastcol;
  dcsr_t *dmat=NULL;
  ssize_t *rowptr;
  int *rowind;
  float *rowval;
  BDMPI_Status status;
  off64_t fpos;

  BDMPI_Comm_lrank(params->comm, &lrank);
  BDMPI_Comm_lsize(params->comm, &lsize);

  if (mype == 0) {
    if (!gk_fexists(params->filename)) 
      errexit("File %s does not exist!\n", params->filename);
  }

  dmat = (dcsr_t *)gk_malloc(sizeof(dcsr_t), "dmat");
  memset(dmat, 0, sizeof(dcsr_t));

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
  if (gk_read(fd, rowptr, sizeof(ssize_t)*(lnrows+1)) != sizeof(ssize_t)*(lnrows+1))
    gk_errexit(SIGERR, "Failed to read the rowptr from file %s!\n", params->filename);

  /* read the rowind */
  nnz = rowptr[lnrows]-rowptr[0];
  rowind = gk_imalloc(nnz, "rowind");
  fpos = 2*sizeof(int) + sizeof(ssize_t)*(gnrows+1) + sizeof(int)*rowptr[0];
  if (lseek64(fd, fpos, SEEK_SET) == -1)
    gk_errexit(SIGERR, "Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
  if (gk_read(fd, rowind, sizeof(int)*nnz) != sizeof(int)*nnz)
    gk_errexit(SIGERR, "Failed to read the rowind from file %s!\n", params->filename);

  /* read the rowval */
  nnz = rowptr[lnrows]-rowptr[0];
  rowval = gk_fmalloc(nnz, "rowval");
  fpos = 2*sizeof(int) + sizeof(ssize_t)*(gnrows+1) + sizeof(int)*gnnz + sizeof(float)*rowptr[0];
  if (lseek64(fd, fpos, SEEK_SET) == -1)
    gk_errexit(SIGERR, "Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
  if (gk_read(fd, rowval, sizeof(float)*nnz) != sizeof(float)*nnz)
    gk_errexit(SIGERR, "Failed to read the rowval from file %s!\n", params->filename);

  /* localize adjust rowptr */
  for (i=lnrows; i>0; i--)
    rowptr[i] -= rowptr[0];
  rowptr[0] = 0;

  close(fd);

  /* tell the next processor to have a go */
  if (lrank != lsize-1)
    BDMPI_Send(&token, 1, BDMPI_INT, mype+1, 1, params->comm);

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

  printf("[%3d, %3d] dmat->gnrows/lnrows: %d/%d, dmat->gncols/lncols: %d/%d, dmat->gnnz/lnnz: %zu/%zu\n", 
      myrow, mycol, 
      dmat->gnrows, dmat->mat->nrows, 
      dmat->gncols, dmat->mat->ncols, 
      dmat->gnnz, dmat->mat->rowptr[dmat->mat->nrows]);

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
/*! Reads a sparse matrix in binary CSR format, one process at a time. 
    \returns the local portion of the matrix.
*/
/**************************************************************************/
dcsr_t *LoadData2D1(params_t *params)
{
  int mype=params->mype, npes=params->npes, 
      myrow=params->myrow, mycol=params->mycol,
      nrowpes=params->nrowpes, ncolpes=params->ncolpes,
      lrank, lsize, token=1;
  size_t i, j, iR, jR, gnnz, lnnz, nnz, maxnnz, offset, N=10000000;
  int fd, gnrows, gncols, lnrows, lncols, firstcol, lastcol;
  dcsr_t *dmat=NULL;
  ssize_t *rowptr;
  int *rowind, *indbuf;
  float *rowval, *valbuf;
  BDMPI_Status status;
  off64_t fpos;

  BDMPI_Comm_lrank(params->comm, &lrank);
  BDMPI_Comm_lsize(params->comm, &lsize);

  if (mype == 0) {
    if (!gk_fexists(params->filename)) 
      errexit("File %s does not exist!\n", params->filename);
  }

  dmat = (dcsr_t *)gk_malloc(sizeof(dcsr_t), "dmat");
  memset(dmat, 0, sizeof(dcsr_t));

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


  /* read the rowptr */
  rowptr = gk_zmalloc(lnrows+1, "rowptr");
  fpos = 2*sizeof(int) + dmat->rowdist[myrow]*sizeof(ssize_t);
  if (lseek64(fd, fpos, SEEK_SET) == -1)
    gk_errexit(SIGERR, "Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
  if (gk_read(fd, rowptr, sizeof(ssize_t)*(lnrows+1)) != sizeof(ssize_t)*(lnrows+1))
    gk_errexit(SIGERR, "Failed to read the rowptr from file %s!\n", params->filename);


  /* get a reasonable upper bound on the number of local non-zeros */
  maxnnz = 1.1*(rowptr[lnrows]-rowptr[0])/ncolpes;

  /* create and allocate memory for the local matrix */
  dmat->mat         = gk_csr_Create();
  dmat->mat->nrows  = lnrows;
  dmat->mat->ncols  = lncols;
  dmat->mat->rowptr = rowptr;
  dmat->mat->rowind = gk_i32malloc(maxnnz, "dmat->mat->rowind");
  dmat->mat->rowval = gk_fmalloc(maxnnz, "dmat->mat->rowval");

  firstcol = dmat->coldist[mycol];
  lastcol  = dmat->coldist[mycol+1];

  /* read the rows, a few at a time and extract what you need */
  indbuf = gk_imalloc(N, "rowbuf");
  valbuf = gk_fmalloc(N, "valbuf");

  for (lnnz=0, iR=0; iR<lnrows;) {
    /* determine chunk size */
    for (nnz=0, jR=iR; jR<lnrows; jR++) {
      if (rowptr[jR+1]-rowptr[jR] + nnz >= N) 
        break;
      nnz += rowptr[jR+1]-rowptr[jR];
    }
    //printf("[%05d] Reading %zu out of %d\n", mype, jR-iR, lnrows);

    /* read the chunk */
    fpos = 2*sizeof(int) + sizeof(ssize_t)*(gnrows+1) + sizeof(int)*rowptr[0+iR];
    if (lseek64(fd, fpos, SEEK_SET) == -1)
      gk_errexit(SIGERR, "Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
    if (gk_read(fd, indbuf, sizeof(int)*nnz) != sizeof(int)*nnz)
      gk_errexit(SIGERR, "Failed to read the rowind from file %s!\n", params->filename);

    fpos = 2*sizeof(int) + sizeof(ssize_t)*(gnrows+1) + sizeof(int)*gnnz + sizeof(float)*rowptr[0+iR];
    if (lseek64(fd, fpos, SEEK_SET) == -1)
      gk_errexit(SIGERR, "Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
    if (gk_read(fd, valbuf, sizeof(float)*nnz) != sizeof(float)*nnz)
      gk_errexit(SIGERR, "Failed to read the rowval from file %s!\n", params->filename);


    offset = rowptr[iR];
    indbuf -= offset;
    valbuf -= offset;

    for (i=iR; i<jR; i++) {
      if (lnnz+(rowptr[i+1]-rowptr[i]) >= maxnnz) {
        maxnnz += .25*maxnnz;
        dmat->mat->rowind = gk_i32realloc(dmat->mat->rowind, maxnnz, "rowind");
        dmat->mat->rowval = gk_frealloc(dmat->mat->rowval, maxnnz, "rowval");
      }

      for (j=rowptr[i]; j<rowptr[i+1]; j++) {
        if (indbuf[j] >= firstcol && indbuf[j] < lastcol) {
          dmat->mat->rowind[lnnz] = indbuf[j] - firstcol;
          dmat->mat->rowval[lnnz] = valbuf[j];
          lnnz++;
        }
      }
      rowptr[i] = lnnz;
    }
    indbuf += offset;
    valbuf += offset;

    iR = jR;
  }
  SHIFTCSR(i, lnrows, rowptr);

  gk_free((void **)&indbuf, &valbuf, LTERM);

  close(fd);


  /* tell the next processor to have a go */
  if (lrank != lsize-1)
    BDMPI_Send(&token, 1, BDMPI_INT, mype+1, 1, params->comm);


  printf("[%3d, %3d] dmat->gnrows/lnrows: %d/%d, dmat->gncols/lncols: %d/%d, dmat->gnnz/lnnz: %zu/%zu\n", 
      myrow, mycol, 
      dmat->gnrows, dmat->mat->nrows, 
      dmat->gncols, dmat->mat->ncols, 
      dmat->gnnz, dmat->mat->rowptr[dmat->mat->nrows]);

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
/*! Reads a sparse matrix in binary CSR format, one process at a time. 
    This version just reads the entire input file and treats it as a local
    block with the rows/columns getting their "natural" positions in the
    overall matrix.
    \returns the local portion of the matrix.
*/
/**************************************************************************/
dcsr_t *LoadData2D(params_t *params)
{
  int mype=params->mype, npes=params->npes, 
      myrow=params->myrow, mycol=params->mycol,
      nrowpes=params->nrowpes, ncolpes=params->ncolpes,
      lrank, lsize, token=1,
      lnrows, lncols;
  size_t i, j, iR, jR, lnnz;
  dcsr_t *dmat=NULL;
  BDMPI_Status status;

  BDMPI_Comm_lrank(params->comm, &lrank);
  BDMPI_Comm_lsize(params->comm, &lsize);

  if (mype == 0) {
    if (!gk_fexists(params->filename)) 
      errexit("File %s does not exist!\n", params->filename);
  }

  dmat = (dcsr_t *)gk_malloc(sizeof(dcsr_t), "dmat");
  memset(dmat, 0, sizeof(dcsr_t));

  /* wait your turn */
  if (lrank != 0) 
    BDMPI_Recv(&token, 1, BDMPI_INT, mype-1, 1, params->comm, &status);

  /* read the matrix */
  dmat->mat = gk_csr_Read(params->filename, GK_CSR_FMT_BINROW, 1, 0);
  lnrows = dmat->mat->nrows;
  lncols = dmat->mat->ncols;
  lnnz   = dmat->mat->rowptr[lnrows];

  dmat->gnrows = lnrows*nrowpes;
  dmat->gncols = lncols*ncolpes;
  dmat->gnnz   = lnnz*npes;


  /* create and populate the rowdist/coldist */
  dmat->rowdist = gk_imalloc(nrowpes+1, "rowdist");
  for (i=0; i<=nrowpes; i++) 
    dmat->rowdist[i] = i*lnrows;
  
  dmat->coldist = gk_imalloc(ncolpes+1, "coldist");
  for (i=0; i<=ncolpes; i++) 
    dmat->coldist[i] = i*lncols;

  /* tell the next processor to have a go */
  if (lrank != lsize-1)
    BDMPI_Send(&token, 1, BDMPI_INT, mype+1, 1, params->comm);


  printf("[%3d, %3d] dmat->gnrows/lnrows: %d/%d, dmat->gncols/lncols: %d/%d, dmat->gnnz/lnnz: %zu/%zu\n", 
      myrow, mycol, 
      dmat->gnrows, dmat->mat->nrows, 
      dmat->gncols, dmat->mat->ncols, 
      dmat->gnnz, dmat->mat->rowptr[dmat->mat->nrows]);

  /* dmat->mat will not be needed anymore, it can be saved to disk, after
     extracting some key info that is needed about the data */
  dmat->lnrows  = dmat->mat->nrows;
  dmat->lncols  = dmat->mat->ncols;
  dmat->lnnz    = dmat->mat->rowptr[dmat->lnrows];
  dmat->rvalsum = gk_fsum(dmat->lnnz, dmat->mat->rowval, 1);

  BDMPI_Entercritical();
  gk_csr_Write(dmat->mat, params->mfile, GK_CSR_FMT_BINROW, 1, 0);
  BDMPI_Exitcritical();

  gk_csr_Free(&(dmat->mat));

  return dmat;
}


/**************************************************************************/
/*! This function deallocates all the memory that was used */
/**************************************************************************/
void CleanupData(params_t *params, dcsr_t *dmat)
{

  BDMPI_Comm_free(&params->rowcomm);
  BDMPI_Comm_free(&params->colcomm);
  BDMPI_Comm_free(&params->diagcomm);

  gk_csr_Free(&(dmat->mat));
  gk_free((void **)&dmat->rowdist, &dmat->coldist,
                   &dmat->scounts, &dmat->sdispls, 
                   &dmat->rcounts, &dmat->rdispls, &dmat->rinds, 
                   &dmat, LTERM);

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
mf_t *ComputeFactors_FR(params_t *params, dcsr_t *dmat)
{
  int nrowpes=params->nrowpes, ncolpes=params->ncolpes,
      myrow=params->myrow, mycol=params->mycol;
  size_t i, j, iter, block, nsamples, nelmnts;
  int ir, ic, k, nrows, ncols, nfactors, ldu, ldv, nlarge;
  int *rowind;
  ssize_t *rowptr;
  float *rowval;
  mf_t *model;
  double *u, *v, *rb, *cb, uir[params->nfactors];
  double lrnrate, lambda, alpha, mu, r, mae, rmse, gmae, grmse;
  BDMPI_Status status;


  model = (mf_t *)gk_malloc(sizeof(mf_t), "model");
  memset(model, 0, sizeof(mf_t));

  nlarge = gk_max(nrowpes, ncolpes);

  nfactors = params->nfactors;
  lrnrate  = params->lrnrate;
  alpha    = params->alpha;
  lambda   = params->lambda;
  nsamples = dmat->gnnz/params->npes;

  nrows  = dmat->lnrows;
  ncols  = dmat->lncols;

  /* compute model->mu */
  BDMPI_Allreduce(&dmat->rvalsum, &model->mu, 1, BDMPI_DOUBLE, BDMPI_SUM, params->comm);
  mu = model->mu = model->mu/dmat->gnnz;

  /* setup ldu/ldv */
  ldu = ldv = nfactors;

  gk_isrand(101+params->mype);

  /* get into the sgd iterations */
  for (iter=0; iter<params->niters; iter++) {
    if (params->mype == 0)
      printf("Working on iteration: %zu\n", iter);

    for (block=0; block<nlarge; block++) {
      /* if active at this iteration */
      if (mycol == (myrow+block)%nlarge) { /* active at this iteration */
        //BDMPI_Barrier(params->diagcomm);

        printf("[%3d %3d] Computing during %zu.%zu [ts: %d] SP.\n",
            myrow, mycol, iter, block, (int)time(NULL));

        /* load the data */
        BDMPI_Entercritical();
        dmat->mat = gk_csr_Read(params->mfile, GK_CSR_FMT_BINROW, 1, 0);
        rowptr = dmat->mat->rowptr;
        rowind = dmat->mat->rowind;
        rowval = dmat->mat->rowval;

        BDMPI_mlock(rowptr, sizeof(ssize_t)*(nrows+1));
        BDMPI_mlock(rowind, sizeof(int32_t)*rowptr[nrows]);
        BDMPI_mlock(rowval, sizeof(int32_t)*rowval[nrows]);
        BDMPI_Exitcritical();

        /* if first time, allocate and initialize */
        if (iter == 0 && model->rb == NULL) { 
          rb = model->rb = gk_dmalloc(nrows, "rb");
          u  = model->u  = gk_dmalloc(nrows*nfactors, "u");
          BDMPI_mlock(rb, nrows*sizeof(double));
          BDMPI_mlock(u, nrows*nfactors*sizeof(double));
          for (ir=0; ir<nrows; ir++) {
            rb[ir] = (RandomInRange(20000)-10000)*0.000001;
            for (k=0; k<nfactors; k++)
              u[ir*ldu+k] = (RandomInRange(20000)-10000)*0.000001;
          }
        }
        if (model->rb == NULL) {
          BDMPI_Entercritical();
          rb = model->rb = gk_dreadfilebin(params->rbfile, &nelmnts);
          u  = model->u  = gk_dreadfilebin(params->ufile, &nelmnts);
          BDMPI_mlock(rb, nrows*sizeof(double));
          BDMPI_mlock(u, nrows*nfactors*sizeof(double));
          BDMPI_Exitcritical();
          unlink(params->rbfile);
          unlink(params->ufile);
        }

        /* if first time, allocate and initialize */
        if (iter == 0 && model->cb == NULL) { 
          cb = model->cb = gk_dmalloc(ncols, "cb");
          v  = model->v  = gk_dmalloc(ncols*nfactors, "v");
          BDMPI_mlock(cb, ncols*sizeof(double));
          BDMPI_mlock(v, ncols*nfactors*sizeof(double));
          for (ic=0; ic<ncols; ic++) {
            cb[ic] = (RandomInRange(20000)-10000)*0.000001;
            for (k=0; k<nfactors; k++)
              v[ic*ldv+k] = (RandomInRange(20000)-10000)*0.000001;
          }
        }
        if (model->cb == NULL) {
          BDMPI_Entercritical();
          cb = model->cb = gk_dreadfilebin(params->cbfile, &nelmnts);
          v  = model->v  = gk_dreadfilebin(params->vfile, &nelmnts);
          BDMPI_mlock(cb, ncols*sizeof(double));
          BDMPI_mlock(v, ncols*nfactors*sizeof(double));
          BDMPI_Exitcritical();
          unlink(params->cbfile);
          unlink(params->vfile);
        }

        printf("[%3d %3d] Computing during %zu.%zu [ts: %d] P1.\n",
            myrow, mycol, iter, block, (int)time(NULL));

        /* perform SGD with random sampling with replacement */
        mae = rmse = 0.0;
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

        printf("[%3d %3d] Computing during %zu.%zu [ts: %d] P2.\n",
            myrow, mycol, iter, block, (int)time(NULL));

        /* free the read-only data */
        BDMPI_munlock(rowptr, sizeof(ssize_t)*(nrows+1));
        BDMPI_munlock(rowind, sizeof(int32_t)*rowptr[nrows]);
        BDMPI_munlock(rowval, sizeof(int32_t)*rowval[nrows]);
        gk_csr_Free(&(dmat->mat));


        /* send row model data to right */
        //printf("[%3d %3d] Sending row to [%3d %3d]\n", myrow, mycol, myrow, (mycol+1)%ncolpes);
        BDMPI_Send(rb, nrows, BDMPI_DOUBLE, (mycol+1)%ncolpes, 1, params->rowcomm);
        BDMPI_Send(u, nrows*nfactors, BDMPI_DOUBLE, (mycol+1)%ncolpes, 2, params->rowcomm);
        BDMPI_munlock(rb, nrows*sizeof(double));
        BDMPI_munlock(u, nrows*nfactors*sizeof(double));
        gk_free((void **)&model->rb, &model->u, LTERM);
        rb = u = NULL;

        /* send col model data to up */
        //printf("[%3d %3d] Sending col to [%3d %3d]\n", myrow, mycol, (myrow+nrowpes-1)%nrowpes, mycol);
        BDMPI_Send(cb, ncols, BDMPI_DOUBLE, (myrow+nrowpes-1)%nrowpes, 3, params->colcomm);
        BDMPI_Send(v, ncols*nfactors, BDMPI_DOUBLE, (myrow+nrowpes-1)%nrowpes, 4, params->colcomm);
        BDMPI_munlock(cb, ncols*sizeof(double));
        BDMPI_munlock(v, ncols*nfactors*sizeof(double));
        gk_free((void **)&model->cb, &model->v, LTERM);
        cb = v = NULL;

        printf("[%3d %3d] Computing during %zu.%zu [ts: %d] EP.\n",
            myrow, mycol, iter, block, (int)time(NULL));
      }

      /* see if you are receiving the row model data */
      if ((mycol+ncolpes-1)%ncolpes == (myrow+block)%nlarge) { 
        BDMPI_Probe((mycol+ncolpes-1)%ncolpes, 2, params->rowcomm, &status);
        //printf("[%3d %3d] Receiving row from [%3d %3d]\n", myrow, mycol, myrow, (mycol+ncolpes-1)%ncolpes);
        rb = model->rb = gk_dmalloc(nrows, "rb");
        u  = model->u  = gk_dmalloc(nrows*nfactors, "u");
        BDMPI_Recv(rb, nrows, BDMPI_DOUBLE, (mycol+ncolpes-1)%ncolpes, 1, params->rowcomm, &status);
        BDMPI_Recv(u, nrows*nfactors, BDMPI_DOUBLE, (mycol+ncolpes-1)%ncolpes, 2, params->rowcomm, &status);

        BDMPI_Entercritical();
        gk_dwritefilebin(params->rbfile, nrows, rb);
        gk_dwritefilebin(params->ufile, nrows*nfactors, u);
        BDMPI_Exitcritical();
        gk_free((void **)&model->rb, &model->u, LTERM);
        rb = u = NULL;
      }

      /* see if you are receiving the column model data */
      if (mycol == ((myrow+1)%nrowpes+block)%nlarge) { 
        BDMPI_Probe((myrow+1)%nrowpes, 4, params->colcomm, &status);
        //printf("[%3d %3d] Receiving col from [%3d %3d]\n", myrow, mycol, (myrow+1)%nrowpes, mycol);
        cb = model->cb = gk_dmalloc(ncols, "cb");
        v  = model->v  = gk_dmalloc(ncols*nfactors, "v");
        BDMPI_Recv(cb, ncols, BDMPI_DOUBLE, (myrow+1)%nrowpes, 3, params->colcomm, &status);
        BDMPI_Recv(v, ncols*nfactors, BDMPI_DOUBLE, (myrow+1)%nrowpes, 4, params->colcomm, &status);

        BDMPI_Entercritical();
        gk_dwritefilebin(params->cbfile, ncols, cb);
        gk_dwritefilebin(params->vfile, ncols*nfactors, v);
        BDMPI_Exitcritical();
        gk_free((void **)&model->cb, &model->v, LTERM);
        cb = v = NULL;
      }
    }

    BDMPI_Reduce(&mae, &gmae, 1, BDMPI_DOUBLE, BDMPI_SUM, 0, params->comm);
    BDMPI_Reduce(&rmse, &grmse, 1, BDMPI_DOUBLE, BDMPI_SUM, 0, params->comm);
    if (params->mype == 0) {
      gmae = gmae/dmat->gnnz;
      grmse = sqrt(grmse/dmat->gnnz);
      printf("Iter %4zu: MAE: %3.5lf RMSE: %3.5lf\n", iter, gmae, grmse);
    }
  }

  return model;
}


/**************************************************************************/
/*! This function computes the page-rank scores using a push approach */
/**************************************************************************/
mf_t *ComputeFactors_RR(params_t *params, dcsr_t *dmat)
{
  int nrowpes=params->nrowpes, ncolpes=params->ncolpes,
      myrow=params->myrow, mycol=params->mycol;
  size_t i, j, iter, block, nsamples, nelmnts;
  int ir, ic, k, nrows, ncols, nfactors, ldu, ldv, nlarge;
  int *rowind;
  ssize_t *rowptr;
  float *rowval;
  mf_t *model;
  double *u, *v, *rb, *cb, uir[params->nfactors];
  double lrnrate, lambda, alpha, mu, r, mae, rmse, gmae, grmse;
  BDMPI_Status status;


  model = (mf_t *)gk_malloc(sizeof(mf_t), "model");
  memset(model, 0, sizeof(mf_t));

  nlarge = gk_max(nrowpes, ncolpes);

  nfactors = params->nfactors;
  lrnrate  = params->lrnrate;
  alpha    = params->alpha;
  lambda   = params->lambda;
  nsamples = dmat->gnnz/params->npes;

  nrows  = dmat->lnrows;
  ncols  = dmat->lncols;

  /* compute model->mu */
  BDMPI_Allreduce(&dmat->rvalsum, &model->mu, 1, BDMPI_DOUBLE, BDMPI_SUM, params->comm);
  mu = model->mu = model->mu/dmat->gnnz;

  /* setup ldu/ldv */
  ldu = ldv = nfactors;

  gk_isrand(101+params->mype);

  //BDMPI_Barrier(params->comm);
    
  /* get into the sgd iterations */
  for (iter=0; iter<params->niters; iter++) {
    if (params->mype == 0)
      printf("Working on iteration: %zu\n", iter);

    for (block=0; block<nlarge; block++) {
      /* if active at this iteration */
      if (mycol == (myrow+block)%nlarge) { /* active at this iteration */
        //BDMPI_Barrier(params->diagcomm);

        printf("[%3d %3d] Computing during %zu.%zu [ts: %d] SP.\n",
            myrow, mycol, iter, block, (int)time(NULL));

        /* load the data */
        BDMPI_Entercritical();
        dmat->mat = gk_csr_Read(params->mfile, GK_CSR_FMT_BINROW, 1, 0);
        BDMPI_Exitcritical();
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
          BDMPI_Entercritical();
          rb = model->rb = gk_dreadfilebin(params->rbfile, &nelmnts);
          u  = model->u  = gk_dreadfilebin(params->ufile, &nelmnts);
          BDMPI_Exitcritical();
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
        if (model->cb == NULL) {
          BDMPI_Entercritical();
          cb = model->cb = gk_dreadfilebin(params->cbfile, &nelmnts);
          v  = model->v  = gk_dreadfilebin(params->vfile, &nelmnts);
          BDMPI_Exitcritical();
          unlink(params->cbfile);
          unlink(params->vfile);
        }

        printf("[%3d %3d] Computing during %zu.%zu [ts: %d] P1.\n",
            myrow, mycol, iter, block, (int)time(NULL));

        /* perform SGD with random sampling with replacement */
        mae = rmse = 0.0;
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

        printf("[%3d %3d] Computing during %zu.%zu [ts: %d] P2.\n",
            myrow, mycol, iter, block, (int)time(NULL));

        /* free the read-only data */
        gk_csr_Free(&(dmat->mat));


        /* send row model data to right */
        //printf("[%3d %3d] Sending row to [%3d %3d]\n", myrow, mycol, myrow, (mycol+1)%ncolpes);
        BDMPI_Send(rb, nrows, BDMPI_DOUBLE, (mycol+1)%ncolpes, 1, params->rowcomm);
        BDMPI_Send(u, nrows*nfactors, BDMPI_DOUBLE, (mycol+1)%ncolpes, 1, params->rowcomm);

        /* send col model data to up */
        //printf("[%3d %3d] Sending col to [%3d %3d]\n", myrow, mycol, (myrow+nrowpes-1)%nrowpes, mycol);
        BDMPI_Send(cb, ncols, BDMPI_DOUBLE, (myrow+nrowpes-1)%nrowpes, 1, params->colcomm);
        BDMPI_Send(v, ncols*nfactors, BDMPI_DOUBLE, (myrow+nrowpes-1)%nrowpes, 1, params->colcomm);

        gk_free((void **)&model->rb, &model->u, &model->cb, &model->v, LTERM);
        rb = cb = u = v = NULL;

        printf("[%3d %3d] Computing during %zu.%zu [ts: %d] EP.\n",
            myrow, mycol, iter, block, (int)time(NULL));
      }

      /* see if you are receiving the row model data */
      if ((mycol+ncolpes-1)%ncolpes == (myrow+block)%nlarge) { 
        //printf("[%3d %3d] Receiving row from [%3d %3d]\n", myrow, mycol, myrow, (mycol+ncolpes-1)%ncolpes);
        rb = model->rb = gk_dmalloc(nrows, "rb");
        u  = model->u  = gk_dmalloc(nrows*nfactors, "u");
        BDMPI_Recv(rb, nrows, BDMPI_DOUBLE, (mycol+ncolpes-1)%ncolpes, 1, params->rowcomm, &status);
        BDMPI_Recv(u, nrows*nfactors, BDMPI_DOUBLE, (mycol+ncolpes-1)%ncolpes, 1, params->rowcomm, &status);

        BDMPI_Entercritical();
        gk_dwritefilebin(params->rbfile, nrows, rb);
        gk_dwritefilebin(params->ufile, nrows*nfactors, u);
        BDMPI_Exitcritical();
        gk_free((void **)&model->rb, &model->u, LTERM);
        rb = u = NULL;
      }

      /* see if you are receiving the column model data */
      if (mycol == ((myrow+1)%nrowpes+block)%nlarge) { 
        //printf("[%3d %3d] Receiving col from [%3d %3d]\n", myrow, mycol, (myrow+1)%nrowpes, mycol);
        cb = model->cb = gk_dmalloc(ncols, "cb");
        v  = model->v  = gk_dmalloc(ncols*nfactors, "v");
        BDMPI_Recv(cb, ncols, BDMPI_DOUBLE, (myrow+1)%nrowpes, 1, params->colcomm, &status);
        BDMPI_Recv(v, ncols*nfactors, BDMPI_DOUBLE, (myrow+1)%nrowpes, 1, params->colcomm, &status);

        BDMPI_Entercritical();
        gk_dwritefilebin(params->cbfile, ncols, cb);
        gk_dwritefilebin(params->vfile, ncols*nfactors, v);
        BDMPI_Exitcritical();
        gk_free((void **)&model->cb, &model->v, LTERM);
        cb = v = NULL;
      }
    }

    //BDMPI_Barrier(params->comm);
    BDMPI_Reduce(&mae, &gmae, 1, BDMPI_DOUBLE, BDMPI_SUM, 0, params->comm);
    BDMPI_Reduce(&rmse, &grmse, 1, BDMPI_DOUBLE, BDMPI_SUM, 0, params->comm);
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
/*! TODO: Has not been written yet! 
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


