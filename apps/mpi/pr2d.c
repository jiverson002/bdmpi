/*!
\file
\brief A parallel page-rank program using 1D distribution of the graph
\date Started 5/29/2013
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
  int npes, nrowpes, ncolpes; 
  int mype, myrow, mycol;
  int lrank, lsize, rrank;
  MPI_Comm comm, rowcomm, colcomm;
  char *filename;
  int niters;

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


/**************************************************************************/
/* prototypes */
/**************************************************************************/
void SetupComms(params_t *params);
dcsr_t *LoadData(params_t *params);
void WritePR(params_t *params, dcsr_t *dmat, double *prvec);
void SetupData(params_t *params, dcsr_t *dmat);
void CleanupData(params_t *params, dcsr_t *dmat);
double *ComputePR(params_t *params, dcsr_t *dmat);


/**************************************************************************/
/**************************************************************************/
int main(int argc, char **argv)
{
  params_t *params;
  dcsr_t *dmat;
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

  if (argc != 5) {
    if (params->mype == 0)
      fprintf(stderr, "Usage: %s filename nrowpes ncolpes niters\n", argv[0]);

    MPI_Finalize();
    return EXIT_FAILURE;
  }

  params->filename = gk_strdup(argv[1]);
  params->nrowpes  = atoi(argv[2]);
  params->ncolpes  = atoi(argv[3]);
  params->niters   = atoi(argv[4]);

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
  dmat = LoadData(params);
  gk_stopwctimer(params->loadTmr);

  gk_startwctimer(params->setupTmr);
  SetupData(params, dmat);
  gk_stopwctimer(params->setupTmr);

  MPI_Barrier(MPI_COMM_WORLD);

  gk_startwctimer(params->compTmr);
  prvec = ComputePR(params, dmat);
  gk_stopwctimer(params->compTmr);

  //WritePR(params, dmat, prvec);

  CleanupData(params, dmat);

  gk_stopwctimer(params->totalTmr);

  /* print timing stats */
  current = gk_getwctimer(params->setupTmr);
  MPI_Reduce(&current, &max, 1, MPI_DOUBLE, MPI_MAX, 0, params->comm);
  if (params->mype == 0)
    printf(" setupTmr:  %10.4lf\n", max);

  current = gk_getwctimer(params->loadTmr);
  MPI_Reduce(&current, &max, 1, MPI_DOUBLE, MPI_MAX, 0, params->comm);
  if (params->mype == 0)
    printf("  loadTmr:  %10.4lf\n", max);

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

  /*
  MPI_Comm_lrank(params->comm, &params->lrank);
  MPI_Comm_lsize(params->comm, &params->lsize);
  MPI_Comm_rrank(params->comm, &params->rrank);
  */
  params->lrank = 0;
  params->lsize = 1;
  params->rrank = 0;

  printf("[%3d] [nr x nc] = [%d %d]/[%d %d], [myrow,mycol] = [%d %d]/[%d %d], local: [%d %d %d]\n",
      params->mype, 
      params->nrowpes, params->ncolpes, nrowpes, ncolpes,
      myrow, mycol, params->myrow, params->mycol,
      params->lrank, params->lsize, params->rrank); 
}
  

/**************************************************************************/
/*! Reads a sparse matrix in binary CSR format, one process at a time. 
    \returns the local portion of the matrix.
*/
/**************************************************************************/
dcsr_t *LoadData0(params_t *params)
{
  int mype=params->mype, npes=params->npes, token=1;
  size_t i, gnnz, lnnz;
  int fd, gnrows, gncols, lnrows;
  dcsr_t *dmat=NULL;
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

  /* create and populate the rowdist */
  lnrows = (gnrows+npes-1)/npes;
  dmat->rowdist = gk_imalloc(npes+1, "rowdist");
  dmat->rowdist[0] = 0;
  for (i=1; i<npes; i++) 
    dmat->rowdist[i] = dmat->rowdist[i-1] + lnrows;
  dmat->rowdist[npes] = gnrows;
  lnrows = dmat->rowdist[mype+1]-dmat->rowdist[mype];

  /*
  printf("[%3d] rowdist: ", mype);
  for (i=0; i<=npes; i++)
    printf("%d ", dmat->rowdist[i]);
  printf("\n");
  */

  dmat->mat = gk_csr_Create();
  dmat->mat->nrows = lnrows;
  dmat->mat->ncols = gncols;

  /* read the rowptr */
  dmat->mat->rowptr = gk_zmalloc(lnrows+1, "dmat->mat->rowptr");
  fpos = 2*sizeof(int) + dmat->rowdist[mype]*sizeof(ssize_t);
  if (lseek64(fd, fpos, SEEK_SET) == -1)
    gk_errexit(SIGERR, "Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
  if (read(fd, dmat->mat->rowptr, sizeof(ssize_t)*(lnrows+1)) != sizeof(ssize_t)*(lnrows+1))
    gk_errexit(SIGERR, "Failed to read the rowptr from file %s!\n", params->filename);

  /* read the rowind */
  lnnz = dmat->mat->rowptr[lnrows]-dmat->mat->rowptr[0];
  dmat->mat->rowind = gk_imalloc(lnnz, "dmat->mat->rowind");
  fpos = 2*sizeof(int) + sizeof(ssize_t)*(gnrows+1) + sizeof(int)*dmat->mat->rowptr[0];
  if (lseek64(fd, fpos, SEEK_SET) == -1)
    gk_errexit(SIGERR, "Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
  if (read(fd, dmat->mat->rowind, sizeof(int)*lnnz) != sizeof(int)*lnnz)
    gk_errexit(SIGERR, "Failed to read the rowind from file %s!\n", params->filename);

  /* read the rowval */
  lnnz = dmat->mat->rowptr[lnrows]-dmat->mat->rowptr[0];
  dmat->mat->rowval = gk_fmalloc(lnnz, "dmat->mat->rowval");
  fpos = 2*sizeof(int) + sizeof(ssize_t)*(gnrows+1) + sizeof(int)*gnnz + sizeof(float)*dmat->mat->rowptr[0];
  if (lseek64(fd, fpos, SEEK_SET) == -1)
    gk_errexit(SIGERR, "Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
  if (read(fd, dmat->mat->rowval, sizeof(float)*lnnz) != sizeof(float)*lnnz)
    gk_errexit(SIGERR, "Failed to read the rowval from file %s!\n", params->filename);

  /* localize adjust rowptr */
  for (i=lnrows; i>0; i--)
    dmat->mat->rowptr[i] -= dmat->mat->rowptr[0];
  dmat->mat->rowptr[0] = 0;

  close(fd);

  if (mype != npes-1)
    MPI_Send(&token, 1, MPI_INT, mype+1, 1, params->comm);

  printf("[%3d] dmat->gnrows/lnrows: %d/%d, dmat->gncols/lncols: %d/%d, dmat->gnnz/lnnz: %zu/%zu\n", 
      mype, 
      dmat->gnrows, dmat->mat->nrows, 
      dmat->gncols, dmat->mat->ncols, 
      dmat->gnnz, dmat->mat->rowptr[dmat->mat->nrows]);

  return dmat;
}


/**************************************************************************/
/*! Reads a sparse matrix in binary CSR format, one process at a time. 
    \returns the local portion of the matrix.
*/
/**************************************************************************/
dcsr_t *LoadData(params_t *params)
{
  int mype=params->mype, npes=params->npes, token=1;
  size_t i, p, gnnz, lnnz;
  int fd, gnrows, gncols, lnrows;
  dcsr_t *dmat=NULL;
  MPI_Status status;
  off64_t fpos;
  ssize_t *rowptr;

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
    if (read(fd, &gnrows, sizeof(int)) != sizeof(int))
      errexit("Failed to read the nrows from file %s!\n", params->filename);
    if (read(fd, &gncols, sizeof(int)) != sizeof(int))
      errexit("Failed to read the ncols from file %s!\n", params->filename);

    rowptr = gk_zmalloc(gnrows+1, "rowptr");
    if (read(fd, rowptr, sizeof(ssize_t)*(gnrows+1)) != sizeof(ssize_t)*(gnrows+1))
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
  MPI_Bcast(dmat->rowdist, npes+1, MPI_INT, 0, params->comm);


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

  lnrows = dmat->rowdist[mype+1]-dmat->rowdist[mype];

  dmat->mat = gk_csr_Create();
  dmat->mat->nrows = lnrows;
  dmat->mat->ncols = gncols;

  /* read the rowptr */
  dmat->mat->rowptr = gk_zmalloc(lnrows+1, "dmat->mat->rowptr");
  fpos = 2*sizeof(int) + dmat->rowdist[mype]*sizeof(ssize_t);
  if (lseek64(fd, fpos, SEEK_SET) == -1)
    gk_errexit(SIGERR, "Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
  if (read(fd, dmat->mat->rowptr, sizeof(ssize_t)*(lnrows+1)) != sizeof(ssize_t)*(lnrows+1))
    gk_errexit(SIGERR, "Failed to read the rowptr from file %s!\n", params->filename);

  /* read the rowind */
  lnnz = dmat->mat->rowptr[lnrows]-dmat->mat->rowptr[0];
  dmat->mat->rowind = gk_imalloc(lnnz, "dmat->mat->rowind");
  fpos = 2*sizeof(int) + sizeof(ssize_t)*(gnrows+1) + sizeof(int)*dmat->mat->rowptr[0];
  if (lseek64(fd, fpos, SEEK_SET) == -1)
    gk_errexit(SIGERR, "Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
  if (read(fd, dmat->mat->rowind, sizeof(int)*lnnz) != sizeof(int)*lnnz)
    gk_errexit(SIGERR, "Failed to read the rowind from file %s! [%s]\n", params->filename, strerror(errno));

  /* read the rowval */
  lnnz = dmat->mat->rowptr[lnrows]-dmat->mat->rowptr[0];
  dmat->mat->rowval = gk_fmalloc(lnnz, "dmat->mat->rowval");
  fpos = 2*sizeof(int) + sizeof(ssize_t)*(gnrows+1) + sizeof(int)*gnnz + sizeof(float)*dmat->mat->rowptr[0];
  if (lseek64(fd, fpos, SEEK_SET) == -1)
    gk_errexit(SIGERR, "Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
  if (read(fd, dmat->mat->rowval, sizeof(float)*lnnz) != sizeof(float)*lnnz)
    gk_errexit(SIGERR, "Failed to read the rowval from file %s!\n", params->filename);

  /* localize adjust rowptr */
  for (i=lnrows; i>0; i--)
    dmat->mat->rowptr[i] -= dmat->mat->rowptr[0];
  dmat->mat->rowptr[0] = 0;

  close(fd);

  if (mype != npes-1)
    MPI_Send(&token, 1, MPI_INT, mype+1, 1, params->comm);

  printf("[%3d] dmat->gnrows/lnrows: %d/%d, dmat->gncols/lncols: %d/%d, dmat->gnnz/lnnz: %zu/%zu\n", 
      mype, 
      dmat->gnrows, dmat->mat->nrows, 
      dmat->gncols, dmat->mat->ncols, 
      dmat->gnnz, dmat->mat->rowptr[dmat->mat->nrows]);

  return dmat;
}


/**************************************************************************/
/*! Writes the page-rank vector. It just let each process write its portion
    to the file in a round-robin fashion.
*/
/**************************************************************************/
void WritePR(params_t *params, dcsr_t *dmat, double *prvec)
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


/**************************************************************************/
/*! This function setups the various data-structures for communication and
    then proceeds to create a stochastic matrix for PR calculations.
*/
/**************************************************************************/
void SetupData(params_t *params, dcsr_t *dmat)
{
  int npes=params->npes, mype=params->mype;
  int nrows, firstrow, lastrow, nunique;
  size_t i, j, p, npairs;
  ssize_t *rowptr;
  int *rowind;
  float *rowval, sum;
  gk_ikv_t *pairs;


  nrows  = dmat->mat->nrows;
  rowptr = dmat->mat->rowptr;
  rowind = dmat->mat->rowind;
  rowval = dmat->mat->rowval;

  firstrow = dmat->rowdist[mype];
  lastrow  = dmat->rowdist[mype+1];


  /* normalize the weights of the matrix to make it stochastic */
  for (i=0; i<nrows; i++) {
    for (sum=0.0, j=rowptr[i]; j<rowptr[i+1]; j++)
      sum += rowval[j];
    if (sum > 0) {
      for (sum=1.0/sum, j=rowptr[i]; j<rowptr[i+1]; j++)
        rowval[j] *= sum;
    }
  }

  /* ============================================================ */
  /* SETUP SINFO */
  /* ============================================================ */

  /* determine the number of non-local non-zeros */
  for (npairs=0, i=0; i<rowptr[nrows]; i++) {
    if (rowind[i] < firstrow || rowind[i] >= lastrow)
      npairs++;
  }

  /* put them in a gk_ikv_t type so that you can renumber them afterwards */
  pairs = gk_ikvmalloc(npairs, "npairs");
  for (npairs=0, i=0; i<rowptr[nrows]; i++) {
    if (rowind[i] < firstrow || rowind[i] >= lastrow) {
      pairs[npairs].key = rowind[i];
      pairs[npairs].val = i;
      npairs++;
    }
    else { /* renumber the local index */
      rowind[i] -= firstrow;
    }
  }
  gk_ikvsorti(npairs, pairs);

  /* determine the unique remote indices and renumber them */
  rowind[pairs[0].val] = nrows;
  for (nunique=0, i=1; i<npairs; i++) {
    if (pairs[i-1].key != pairs[i].key)
      nunique++;
    rowind[pairs[i].val] = nrows + nunique;
  }
  nunique++;

  /* ============================================================= */
  /* dmat->mat will not be needed anymore, it can be saved to disk */
  /* ============================================================= */

  /* allocate memory for the sinfo (we implement a push-based algorithm) */
  dmat->nsend   = nunique;
  dmat->scounts = gk_imalloc(npes, "scounts");
  dmat->sdispls = gk_imalloc(npes+1, "sdispls");
  dmat->sinds   = gk_imalloc(nunique, "sinds");

  /* copy the unique indices into dmat->sinds */
  dmat->sinds[0] = pairs[0].key;
  for (nunique=1, i=1; i<npairs; i++) {
    if (pairs[i-1].key != pairs[i].key)
      dmat->sinds[nunique++] = pairs[i].key;
  }

  /* determine scounts/sdispls */
  dmat->sdispls[0] = 0;
  for (i=0, p=0; p<npes; p++) {
    for (; i<nunique; i++) {
      if (dmat->sinds[i] >= dmat->rowdist[p+1])
        break;
    }
    dmat->sdispls[p+1] = i;
    dmat->scounts[p] = dmat->sdispls[p+1]-dmat->sdispls[p];
  }

  gk_free((void **)&pairs, LTERM);


  /* ============================================================ */
  /* SETUP RINFO */
  /* ============================================================ */

  /* allocate memory for rcounts and perform an all-to-all to get the data */
  dmat->rcounts = gk_imalloc(npes, "rcounts");
  MPI_Alltoall(dmat->scounts, 1, MPI_INT, dmat->rcounts, 1, MPI_INT, params->comm);

  /* allocate memory for rdispls and fill it */
  dmat->rdispls = gk_imalloc(npes+1, "rdispls");
  dmat->rdispls[0] = 0;
  for (i=0; i<npes; i++)
    dmat->rdispls[i+1] = dmat->rdispls[i] + dmat->rcounts[i];
  dmat->nrecv = dmat->rdispls[npes];


  /* allocate memory for rinds and populate it via an all-to-all */
  dmat->rinds = gk_imalloc(dmat->nrecv, "rinds");
  MPI_Alltoallv(dmat->sinds, dmat->scounts, dmat->sdispls, MPI_INT, 
                dmat->rinds, dmat->rcounts, dmat->rdispls, MPI_INT, 
                params->comm);

  /* free sinds, as they will not be used again */
  gk_free((void **)&dmat->sinds, LTERM);

  /* localize the indices in dmat->rinds */
  for (i=0; i<dmat->nrecv; i++) {
    if (dmat->rinds[i] < firstrow || dmat->rinds[i] >= lastrow) {
      gk_errexit(SIGERR, "[%d] rinds[%zu]=%d is outside local range [%d %d]\n", 
          mype, i, dmat->rinds[i], firstrow, lastrow);
    }
    dmat->rinds[i] -= firstrow;
  }

  printf("[%3d] nsend: %d, nrecv: %d\n", mype, dmat->nsend, dmat->nrecv);

  return;
}


/**************************************************************************/
/*! This function deallocates all the memory that was used */
/**************************************************************************/
void CleanupData(params_t *params, dcsr_t *dmat)
{

  gk_csr_Free(&(dmat->mat));
  gk_free((void **)&dmat->rowdist,
                   &dmat->scounts, &dmat->sdispls, 
                   &dmat->rcounts, &dmat->rdispls, &dmat->rinds, 
                   &dmat, LTERM);

  return;
}


/**************************************************************************/
/*! This function computes the page-rank scores using a push approach */
/**************************************************************************/
double *ComputePR(params_t *params, dcsr_t *dmat)
{
  int npes=params->npes, mype=params->mype;
  size_t iter, i, j, nrows, nsend, nrecv;
  int *rowind, *rinds;
  ssize_t *rowptr;
  float *rowval;
  double *pr, *prnew, *prrecv;
  double lambda=.2, rprob, lrmsd, grmsd, lsinks, gsinks;

  nsend = dmat->nsend;
  nrecv = dmat->nrecv;
  rinds = dmat->rinds;

  nrows  = dmat->mat->nrows;
  rowptr = dmat->mat->rowptr;
  rowind = dmat->mat->rowind;
  rowval = dmat->mat->rowval;

  rprob = 1.0/dmat->gnrows;
  pr    = gk_dsmalloc(nrows, rprob, "pr");

  /* get into the PR iteration */
  for (iter=0; iter<=params->niters; iter++) {
    prnew = gk_dsmalloc(nrows+nsend, 0.0, "prnew");

    /* push random-walk scores to the outlinks */
    for (lsinks=0.0, i=0; i<nrows; i++) {
      if (rowptr[i+1]-rowptr[i] == 0) {
        lsinks += pr[i];
        continue;
      }

      for (j=rowptr[i]; j<rowptr[i+1]; j++) {
        if (rowind[j] < 0 || rowind[j] >= nrows+nsend)
          printf("[%3d] rowind[%zd]=%d is out of range %zd.\n", mype, j, rowind[j], nrows+nsend);
        prnew[rowind[j]] += pr[i]*rowval[j];
      }
    }

    gk_startwctimer(params->commTmr);
    /* get the overall gsinks across all processors */
    MPI_Allreduce(&lsinks, &gsinks, 1, MPI_DOUBLE, MPI_SUM, params->comm);

    /* get the partial PR scores computed by the other PEs for my rows and
       update local PR scores */
    prrecv = gk_dmalloc(nrecv, "prrecv");
    MPI_Alltoallv(prnew+nrows, dmat->scounts, dmat->sdispls, MPI_DOUBLE,
                  prrecv, dmat->rcounts, dmat->rdispls, MPI_DOUBLE,
                  params->comm);
    gk_stopwctimer(params->commTmr);

    /* shrink the size of prnew, as you do not need the last nsend entries */
    prnew = gk_drealloc(prnew, nrows, "drealloc: prnew");

    for (i=0; i<nrecv; i++)
      prnew[rinds[i]] += prrecv[i];

    gk_free((void **)&prrecv, LTERM);

    /* apply the restart condition */
    for (i=0; i<nrows; i++)
      prnew[i] = lambda*rprob + (1.0-lambda)*(gsinks*rprob+prnew[i]);

    if (iter%10 == 0) {
      /* compute the difference */
      for (lrmsd=0.0, i=0; i<nrows; i++) {
        lrmsd += (pr[i]-prnew[i])*(pr[i]-prnew[i]);
        pr[i] = prnew[i];
      }
      gk_free((void **)&prnew, LTERM);

      /* get the global rmsd across all processors */
      gk_startwctimer(params->commTmr);
      MPI_Allreduce(&lrmsd, &grmsd, 1, MPI_DOUBLE, MPI_SUM, params->comm);
      grmsd = sqrt(grmsd);
      gk_stopwctimer(params->commTmr);

      if (mype == 0)
        printf("Iter: %5zu, grmsd: %.6le\n", iter, grmsd);
    }
    else {
      for (i=0; i<nrows; i++) 
        pr[i] = prnew[i];
      gk_free((void **)&prnew, LTERM);
    }
  }

  return pr;
}

