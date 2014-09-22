/*!
\file
\brief A parallel page-rank program using 1D distribution of the graph
\date Started 5/29/2013
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
  int npes, mype;
  BDMPI_Comm comm;
  char *filename;
  int niters;
  int type;

  /* oc filenames */
  char mfile[1024];
  char pfile[1024];

  /* timers */
  double totalTmr;
  double setupTmr;
  double oocTmr;
  double compTmr;
  double commTmr;

} params_t;

/* distributed CSR */
typedef struct {
  /* the total number of rows/non-zeros and their overall distribution */
  int gnrows, gncols, lnrows;
  size_t gnnz;
  int *rowdist;

  /* sendinfo */
  int nsend;
  size_t *scounts;
  size_t *sdispls;
  int *sinds;

  /* recvinfo */
  int nrecv;
  size_t *rcounts;
  size_t *rdispls;
  int *rinds;

  gk_csr_t *mat;
} dcsr_t;


/**************************************************************************/
/* prototypes */
/**************************************************************************/
dcsr_t *LoadData(params_t *params);
void WritePR(params_t *params, dcsr_t *dmat, double *prvec);
void SetupData(params_t *params, dcsr_t *dmat);
void CleanupData(params_t *params, dcsr_t *dmat);
double *ComputePR_a2a(params_t *params, dcsr_t *dmat);
double *ComputePR_p2p(params_t *params, dcsr_t *dmat);


/**************************************************************************/
/**************************************************************************/
int main(int argc, char **argv)
{
  params_t *params;
  dcsr_t *dmat;
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
      fprintf(stderr, "Usage: %s filename niters type[1:a2a; 2:p2p]\n", argv[0]);

    BDMPI_Finalize();
    return EXIT_FAILURE;
  }

  params->filename = gk_strdup(argv[1]);
  params->niters   = atoi(argv[2]);
  params->type     = atoi(argv[3]);

  sprintf(params->mfile, "m%d.%d", (int)getpid(), params->mype);
  sprintf(params->pfile, "pr%d.%d", (int)getpid(), params->mype);

  gk_clearwctimer(params->totalTmr);
  gk_clearwctimer(params->setupTmr);
  gk_clearwctimer(params->oocTmr);
  gk_clearwctimer(params->compTmr);
  gk_clearwctimer(params->commTmr);

  BDMPI_Barrier(params->comm);
  BDMPI_Barrier(params->comm);
  gk_startwctimer(params->totalTmr);

  dmat = LoadData(params);

  gk_startwctimer(params->setupTmr);
  SetupData(params, dmat);
  gk_stopwctimer(params->setupTmr);

  gk_startwctimer(params->compTmr);
  if (params->type == 1)
    prvec = ComputePR_a2a(params, dmat);
  else
    prvec = ComputePR_p2p(params, dmat);
  gk_stopwctimer(params->compTmr);

  //WritePR(params, dmat, prvec);

  CleanupData(params, dmat);

  BDMPI_Barrier(params->comm);
  BDMPI_Barrier(params->comm);
  gk_stopwctimer(params->totalTmr);

  /* print timing stats */
  current = gk_getwctimer(params->setupTmr);
  BDMPI_Reduce(&current, &max, 1, BDMPI_DOUBLE, BDMPI_MAX, 0, params->comm);
  if (params->mype == 0)
    printf(" setupTmr:  %10.4lf\n", max);

  current = gk_getwctimer(params->compTmr);
  BDMPI_Reduce(&current, &max, 1, BDMPI_DOUBLE, BDMPI_MAX, 0, params->comm);
  if (params->mype == 0)
    printf("  compTmr:  %10.4lf\n", max);

  current = gk_getwctimer(params->oocTmr);
  BDMPI_Reduce(&current, &max, 1, BDMPI_DOUBLE, BDMPI_MAX, 0, params->comm);
  if (params->mype == 0)
    printf("   oocTmr:  %10.4lf\n", max);

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
/*! Reads a sparse matrix in binary CSR format, one process at a time. 
    \returns the local portion of the matrix.
*/
/**************************************************************************/
dcsr_t *LoadData(params_t *params)
{
  int mype=params->mype, npes=params->npes, token=1;
  int lrank, lsize;
  size_t i, p, gnnz, lnnz;
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
  dmat->lnrows = lnrows = dmat->rowdist[mype+1]-dmat->rowdist[mype];

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
  if (gk_read(fd, dmat->mat->rowind, sizeof(int)*lnnz) != sizeof(int)*lnnz)
    gk_errexit(SIGERR, "Failed to read the rowind from file %s!\n", params->filename);

#ifdef XXX
  /* read the rowval */
  lnnz = dmat->mat->rowptr[lnrows]-dmat->mat->rowptr[0];
  dmat->mat->rowval = gk_fmalloc(lnnz, "dmat->mat->rowval");
  fpos = 2*sizeof(int) + sizeof(ssize_t)*(gnrows+1) + sizeof(int)*gnnz + sizeof(float)*dmat->mat->rowptr[0];
  if (lseek64(fd, fpos, SEEK_SET) == -1)
    gk_errexit(SIGERR, "Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
  if (gk_read(fd, dmat->mat->rowval, sizeof(float)*lnnz) != sizeof(float)*lnnz)
    gk_errexit(SIGERR, "Failed to read the rowval from file %s!\n", params->filename);
#endif

  /* localize adjust rowptr */
  for (i=lnrows; i>0; i--)
    dmat->mat->rowptr[i] -= dmat->mat->rowptr[0];
  dmat->mat->rowptr[0] = 0;

  close(fd);

  if (lrank != lsize-1)
    BDMPI_Send(&token, 1, BDMPI_INT, mype+1, 1, params->comm);

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
  BDMPI_Status status;
  FILE *fpout;
  char outfile[1024];

  sprintf(outfile, "%s.prvec", params->filename);

  if (mype == 0) {
    fpout = gk_fopen(outfile, "w", "outfile");
    for (i=0; i<dmat->lnrows; i++)
      fprintf(fpout, "%.8le\n", prvec[i]);
    gk_fclose(fpout);
    if (mype+1 < npes)
      BDMPI_Send(&dummy, 1, BDMPI_INT, mype+1, 1, params->comm);
  }
  else {
    BDMPI_Recv(&dummy, 1, BDMPI_INT, mype-1, 1, params->comm, &status);
    fpout = gk_fopen(outfile, "a", "outfile");
    for (i=0; i<dmat->lnrows; i++)
      fprintf(fpout, "%lf\n", prvec[i]);
    gk_fclose(fpout);
    if (mype+1 < npes)
      BDMPI_Send(&mype, 1, BDMPI_INT, mype+1, 1, params->comm);
  }
}


/**************************************************************************/
/*! This function setups the various data-structures for communication and
    then proceeds to create a stochastic matrix for PR calculations.
*/
/**************************************************************************/
void SetupData0(params_t *params, dcsr_t *dmat)
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

#ifdef XXX
  /* normalize the weights of the matrix to make it stochastic */
  for (i=0; i<nrows; i++) {
    for (sum=0.0, j=rowptr[i]; j<rowptr[i+1]; j++)
      sum += rowval[j];
    if (sum > 0) {
      for (sum=1.0/sum, j=rowptr[i]; j<rowptr[i+1]; j++)
        rowval[j] *= sum;
    }
  }
#endif

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

  /* dmat->mat will not be needed anymore, it can be saved to disk */
  gk_csr_Write(dmat->mat, params->mfile, GK_CSR_FMT_BINROW, 0, 0);
  gk_csr_Free(&(dmat->mat));

  /* allocate memory for the sinfo (we implement a push-based algorithm) */
  dmat->nsend   = nunique;
  dmat->scounts = gk_zumalloc(npes, "scounts");
  dmat->sdispls = gk_zumalloc(npes+1, "sdispls");
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
  dmat->rcounts = gk_zumalloc(npes, "rcounts");
  BDMPI_Alltoall(dmat->scounts, 1, BDMPI_SIZE_T, dmat->rcounts, 1, BDMPI_SIZE_T, params->comm);

  /* allocate memory for rdispls and fill it */
  dmat->rdispls = gk_zumalloc(npes+1, "rdispls");
  dmat->rdispls[0] = 0;
  for (i=0; i<npes; i++)
    dmat->rdispls[i+1] = dmat->rdispls[i] + dmat->rcounts[i];
  dmat->nrecv = dmat->rdispls[npes];


  /* allocate memory for rinds and populate it via an all-to-all */
  dmat->rinds = gk_imalloc(dmat->nrecv, "rinds");
  BDMPI_Alltoallv(dmat->sinds, dmat->scounts, dmat->sdispls, BDMPI_INT, 
                dmat->rinds, dmat->rcounts, dmat->rdispls, BDMPI_INT, 
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
/*! This function setups the various data-structures for communication and
    then proceeds to create a stochastic matrix for PR calculations.
*/
/**************************************************************************/
void SetupData(params_t *params, dcsr_t *dmat)
{
  int npes=params->npes, mype=params->mype, p;
  int nrows, firstrow, lastrow, nunique, pnunique, pfirstrow, plastrow;
  size_t i, j, k;
  ssize_t *rowptr, *nrowptr;
  int *rowind, *rowdist, *marker;

  rowdist = dmat->rowdist;

  nrows  = dmat->mat->nrows;
  rowptr = dmat->mat->rowptr;
  rowind = dmat->mat->rowind;

  firstrow = rowdist[mype];
  lastrow  = rowdist[mype+1];

#ifdef XXX
  /* normalize the weights of the matrix to make it stochastic */
  for (i=0; i<nrows; i++) {
    for (sum=0.0, j=rowptr[i]; j<rowptr[i+1]; j++)
      sum += rowval[j];
    if (sum > 0) {
      for (sum=1.0/sum, j=rowptr[i]; j<rowptr[i+1]; j++)
        rowval[j] *= sum;
    }
  }
#endif


  /* ============================================================ */
  /* SETUP SINFO AND RENUMBER LOCAL INDICES */
  /* ============================================================ */

  /* sort the adjacency list of each vertex */
  for (i=0; i<nrows; i++) {
    if (rowptr[i+1]-rowptr[i] > 1)
      gk_isorti(rowptr[i+1]-rowptr[i], rowind+rowptr[i]);
  }


  /* allocate memory for the sendinfo */
  dmat->scounts = gk_zumalloc(npes, "scounts");
  dmat->sdispls = gk_zumalloc(npes+1, "sdispls");
  dmat->sinds   = gk_imalloc(rowdist[npes], "sinds");

  /* setup the nrowptr */
  nrowptr = gk_zmalloc(nrows+1, "nrowptr");
  gk_zcopy(nrows+1, rowptr, nrowptr);

  for (nunique=0, p=0; p<npes; p++) {
    pfirstrow = rowdist[p];
    plastrow  = rowdist[p+1];

    if (p == mype) {
      for (i=0; i<nrows; i++) {
        for (j=nrowptr[i]; j<rowptr[i+1]; j++) {
          if (rowind[j] < plastrow) 
            rowind[j] -= pfirstrow;
          else
            break;
        }
        nrowptr[i] = j;
      }
      dmat->scounts[p] = 0;
      dmat->sdispls[p] = nunique;
    }
    else {
      marker  = gk_ismalloc(plastrow-pfirstrow, -1, "marker");

      for (pnunique=0, i=0; i<nrows; i++) {
        for (j=nrowptr[i]; j<rowptr[i+1]; j++) {
          if (rowind[j] < plastrow) {
            k = rowind[j]-pfirstrow;
            if (marker[k] == -1) {
              dmat->sinds[nunique+pnunique] = rowind[j];
              marker[k] = pnunique++;
            }
            rowind[j] = nrows + nunique + marker[k];
          }
          else 
            break;
        }
        nrowptr[i] = j;
      }

      dmat->scounts[p] = pnunique;
      dmat->sdispls[p] = nunique;
      nunique += pnunique;

      gk_free((void **)&marker, LTERM);
    }
  }
  dmat->nsend = nunique;

  gk_free((void **)&nrowptr, LTERM);

  /* adjust the memory for sinds */
  dmat->sinds = gk_irealloc(dmat->sinds, nunique, "sinds");


  /* ============================================================= */
  /* dmat->mat will not be needed anymore, it can be saved to disk */
  /* ============================================================= */
  gk_csr_Write(dmat->mat, params->mfile, GK_CSR_FMT_BINROW, 0, 0);
  gk_csr_Free(&(dmat->mat));


  /* ============================================================ */
  /* SETUP RINFO */
  /* ============================================================ */

  /* allocate memory for rcounts and perform an all-to-all to get the data */
  dmat->rcounts = gk_zumalloc(npes, "rcounts");
  BDMPI_Alltoall(dmat->scounts, 1, BDMPI_SIZE_T, dmat->rcounts, 1, BDMPI_SIZE_T, params->comm);

  /* allocate memory for rdispls and fill it */
  dmat->rdispls = gk_zumalloc(npes+1, "rdispls");
  dmat->rdispls[0] = 0;
  for (i=0; i<npes; i++)
    dmat->rdispls[i+1] = dmat->rdispls[i] + dmat->rcounts[i];
  dmat->nrecv = dmat->rdispls[npes];


  /* allocate memory for rinds and populate it via an all-to-all */
  dmat->rinds = gk_imalloc(dmat->nrecv, "rinds");
  BDMPI_Alltoallv(dmat->sinds, dmat->scounts, dmat->sdispls, BDMPI_INT, 
                dmat->rinds, dmat->rcounts, dmat->rdispls, BDMPI_INT, 
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

  printf("[%3d] nsend: %d, nrecv: %d [ts: %d]\n", mype, dmat->nsend, 
      dmat->nrecv,(int)time(NULL));

  return;
}


/**************************************************************************/
/*! This function deallocates all the memory that was used */
/**************************************************************************/
void CleanupData(params_t *params, dcsr_t *dmat)
{

  gk_free((void **)&dmat->rowdist,
                   &dmat->scounts, &dmat->sdispls, 
                   &dmat->rcounts, &dmat->rdispls, &dmat->rinds, 
                   &dmat, LTERM);

  unlink(params->mfile);
  unlink(params->pfile);

  return;
}


/**************************************************************************/
/*! This function computes the page-rank scores using a push approach */
/**************************************************************************/
double *ComputePR_a2a(params_t *params, dcsr_t *dmat)
{
  int npes=params->npes, mype=params->mype;
  size_t iter, i, j, nrows, nsend, nrecv, nelmnts;
  int *rowind, *rinds;
  ssize_t *rowptr;
  float *rowval;
  double *pr, *prnew, *prrecv;
  double wgt, lambda=.2, rprob, lrmsd, grmsd, lsinks, gsinks=1.0;
  gk_csr_t *mat;

  nsend = dmat->nsend;
  nrecv = dmat->nrecv;
  rinds = dmat->rinds;

  nrows  = dmat->lnrows;


  rprob = 1.0/dmat->gnrows;
  pr    = gk_dsmalloc(nrows, rprob, "pr");

  /* get into the PR iteration */
  for (iter=0; iter<params->niters; iter++) {
    mat = gk_csr_Read(params->mfile, GK_CSR_FMT_BINROW, 0, 0); 
    rowptr = mat->rowptr;
    rowind = mat->rowind;
    rowval = mat->rowval;

    prnew = gk_dsmalloc(nrows+nsend, 0.0, "prnew");

    /* push random-walk scores to the outlinks */
    for (lsinks=0.0, i=0; i<nrows; i++) {
      if (rowptr[i+1]-rowptr[i] == 0) {
        lsinks += pr[i];
        continue;
      }
      else {
        wgt = pr[i]/(rowptr[i+1]-rowptr[i]);
        for (j=rowptr[i]; j<rowptr[i+1]; j++) {
          prnew[rowind[j]] += wgt;
        }
      }
    }
    gk_csr_Free(&mat);

    /* write and free the pr file */
    if (iter%10 == 0 || iter == params->niters-1) 
      gk_dwritefilebin(params->pfile, nrows, pr);
    gk_free((void **)&pr, LTERM);


    gk_startwctimer(params->commTmr);
    /* get the overall gsinks across all processors */
    if (gsinks > 0)
      BDMPI_Allreduce(&lsinks, &gsinks, 1, BDMPI_DOUBLE, BDMPI_SUM, params->comm);

    /* get the partial PR scores computed by the other PEs for my rows and
       update local PR scores */
    prrecv = gk_dmalloc(nrecv, "prrecv");
    BDMPI_Alltoallv(prnew+nrows, dmat->scounts, dmat->sdispls, BDMPI_DOUBLE,
                  prrecv, dmat->rcounts, dmat->rdispls, BDMPI_DOUBLE,
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

    if (iter%10 == 0 || iter == params->niters-1) {
      /* read back the old pr file */
      pr = gk_dreadfilebin(params->pfile, &nelmnts);
      if (nelmnts != nrows)
        errexit("The read pr file is not of proper length.\n");

      /* compute the difference */
      for (lrmsd=0.0, i=0; i<nrows; i++) 
        lrmsd += (pr[i]-prnew[i])*(pr[i]-prnew[i]);

      gk_free((void **)&pr, LTERM);

      /* get the global rmsd across all processors */
      gk_startwctimer(params->commTmr);
      BDMPI_Allreduce(&lrmsd, &grmsd, 1, BDMPI_DOUBLE, BDMPI_SUM, params->comm);
      grmsd = sqrt(grmsd);
      gk_stopwctimer(params->commTmr);

      if (mype == 0)
        printf("Iter: %5zu, grmsd: %.6le [gsinks: %.6le, rprob: %.6le]\n",
            iter, grmsd, gsinks, rprob);
    }
    pr = prnew;
  }

  return pr;
}


/**************************************************************************/
/*! This function computes the page-rank scores using a push approach */
/**************************************************************************/
double *ComputePR_p2p(params_t *params, dcsr_t *dmat)
{
  int npes=params->npes, mype=params->mype;
  size_t iter, i, j, p, nrows, nsend, nrecv, nelmnts;
  int *rowind, *rinds;
  ssize_t *rowptr;
  float *rowval;
  double *pr1, *pr2, *tpr;
  double wgt, lambda=.2, rprob, lrmsd, grmsd, lsinks, gsinks=1.0;
  gk_csr_t *mat;
  BDMPI_Status status;

  nsend = dmat->nsend;
  nrecv = dmat->nrecv;
  rinds = dmat->rinds;

  nrows  = dmat->lnrows;

  rprob = 1.0/dmat->gnrows;
  pr1   = gk_dsmalloc(nrows, rprob, "pr1");

  /* get into the PR iteration */
  for (iter=0; iter<params->niters; iter++) {
    gk_startwctimer(params->oocTmr);
    mat = gk_csr_Read(params->mfile, GK_CSR_FMT_BINROW, 0, 0); 
    gk_stopwctimer(params->oocTmr);

    rowptr = mat->rowptr;
    rowind = mat->rowind;
    rowval = mat->rowval;

    tpr = gk_dsmalloc(nrows+nsend, 0.0, "tpr");

    /* push random-walk scores to the outlinks */
    gk_startwctimer(params->compTmr);
    for (lsinks=0.0, i=0; i<nrows; i++) {
      if (rowptr[i+1]-rowptr[i] == 0) {
        lsinks += pr1[i];
        continue;
      }
      else {
        wgt = pr1[i]/(rowptr[i+1]-rowptr[i]);
        for (j=rowptr[i]; j<rowptr[i+1]; j++) {
          tpr[rowind[j]] += wgt;
        }
      }
    }
    gk_stopwctimer(params->compTmr);

    gk_csr_Free(&mat);

    /* write and free the pr1 file */
    if (iter%10 == 0 || iter == params->niters-1) {
      gk_startwctimer(params->oocTmr);
      gk_dwritefilebin(params->pfile, nrows, pr1);
      gk_stopwctimer(params->oocTmr);
    }
    gk_free((void **)&pr1, LTERM);


    /* get the overall gsinks across all processors */
    if (gsinks > 0)
      BDMPI_Allreduce(&lsinks, &gsinks, 1, BDMPI_DOUBLE, BDMPI_SUM, params->comm);

    /* send your partial results */
    tpr += nrows;
    GKASSERT(BDMPI_mlock(tpr, nsend*sizeof(double)) == 0);
    for (j=1; j<npes; j++) {
      p = (mype+j)%npes;
      BDMPI_Send(tpr+dmat->sdispls[p], dmat->scounts[p], BDMPI_DOUBLE, p, 2000, params->comm);
    }
    GKASSERT(BDMPI_munlock(tpr, nsend*sizeof(double)) == 0);
    tpr -= nrows;

    /* keep the partial results for your own vertices in tpr */
    pr2 = gk_dcopy(nrows, tpr, gk_dmalloc(nrows, "pr2"));
    gk_free((void **)&tpr, LTERM);

    /* receive and incorporate remote partial results */
    for (j=1; j<npes; j++) {
      BDMPI_Probe(BDMPI_ANY_SOURCE, 2000, params->comm, &status);
      p = status.BDMPI_SOURCE;

      tpr = gk_dmalloc(nrows, "tpr");
      BDMPI_Recv(tpr, nrows, BDMPI_DOUBLE, p, 2000, params->comm, &status);

      /* incorporate remote updates */
      gk_startwctimer(params->compTmr);
      rinds = dmat->rinds + dmat->rdispls[p];

      GKASSERT(BDMPI_mlock(rinds, dmat->rcounts[p]*sizeof(int)) == 0);
      GKASSERT(BDMPI_mlock(tpr, dmat->rcounts[p]*sizeof(double)) == 0);
      GKASSERT(BDMPI_mlock(pr2, nrows*sizeof(double)) == 0);

      for (i=0; i<dmat->rcounts[p]; i++)
        pr2[rinds[i]] += tpr[i];

      GKASSERT(BDMPI_munlock(rinds, dmat->rcounts[p]*sizeof(int)) == 0);
      GKASSERT(BDMPI_munlock(tpr, dmat->rcounts[p]*sizeof(double)) == 0);
      GKASSERT(BDMPI_munlock(pr2, nrows*sizeof(double)) == 0);

      gk_stopwctimer(params->compTmr);

      gk_free((void **)&tpr, LTERM);
    }

    /* apply the restart condition */
    gk_startwctimer(params->compTmr);
    GKASSERT(BDMPI_mlock(pr2, nrows*sizeof(double)) == 0);
    for (i=0; i<nrows; i++)
      pr2[i] = lambda*rprob + (1.0-lambda)*(gsinks*rprob+pr2[i]);
    GKASSERT(BDMPI_munlock(pr2, nrows*sizeof(double)) == 0);
    gk_stopwctimer(params->compTmr);

    if (iter%10 == 0 || iter == params->niters-1) {
      GKASSERT(BDMPI_mlock(pr2, nrows*sizeof(double)) == 0);

      pr1 = gk_dreadfilebin(params->pfile, &nelmnts);
      if (nelmnts != nrows)
        errexit("The read pr file is not of proper length.\n");

      /* compute the difference */
      for (lrmsd=0.0, i=0; i<nrows; i++) 
        lrmsd += (pr1[i]-pr2[i])*(pr1[i]-pr2[i]);

      GKASSERT(BDMPI_munlock(pr2, nrows*sizeof(double)) == 0);

      gk_free((void **)&pr1, LTERM);

      /* get the global rmsd across all processors */
      BDMPI_Allreduce(&lrmsd, &grmsd, 1, BDMPI_DOUBLE, BDMPI_SUM, params->comm);
      grmsd = sqrt(grmsd);

      if (mype == 0)
        printf("Iter: %5zu, grmsd: %.6le [ts: %d]\n", iter, grmsd, (int)time(NULL));
    }
    pr1 = pr2;
  }

  return pr1;
}

