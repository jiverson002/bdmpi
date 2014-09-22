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
  int niters;

  /* timers */
  double totalTmr;
  double setupTmr;
  double compTmr;
  double commTmr;
  double sort1Tmr;
  double sort2Tmr;
} params_t;

/* distributed CSR */
typedef struct {
  /* the total number of rows/non-zeros and their overall distribution */
  int gnrows, gncols;
  size_t gnnz;
  int *rowdist;

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

  if (argc != 3) {
    if (params->mype == 0)
      fprintf(stderr, "Usage: %s filename niters\n", argv[0]);

    MPI_Finalize();
    return EXIT_FAILURE;
  }

  params->filename = gk_strdup(argv[1]);
  params->niters   = atoi(argv[2]);

  gk_clearwctimer(params->totalTmr);
  gk_clearwctimer(params->setupTmr);
  gk_clearwctimer(params->compTmr);
  gk_clearwctimer(params->commTmr);
  gk_clearwctimer(params->sort1Tmr);
  gk_clearwctimer(params->sort2Tmr);

  MPI_Barrier(params->comm);
  gk_startwctimer(params->totalTmr);

  dmat = LoadData(params);

  gk_startwctimer(params->setupTmr);
  SetupData(params, dmat);
  gk_stopwctimer(params->setupTmr);

  prvec = ComputePR(params, dmat);

  //WritePR(params, dmat, prvec);

  CleanupData(params, dmat);

  MPI_Barrier(params->comm);
  gk_stopwctimer(params->totalTmr);

  /* print timing stats */
  current = gk_getwctimer(params->sort1Tmr);
  MPI_Reduce(&current, &max, 1, MPI_DOUBLE, MPI_MAX, 0, params->comm);
  if (params->mype == 0 && max>0)
    printf(" sort1Tmr:  %10.4lf\n", max);

  current = gk_getwctimer(params->sort2Tmr);
  MPI_Reduce(&current, &max, 1, MPI_DOUBLE, MPI_MAX, 0, params->comm);
  if (params->mype == 0 && max>0)
    printf(" sort2Tmr:  %10.4lf\n", max);

  current = gk_getwctimer(params->setupTmr);
  MPI_Reduce(&current, &max, 1, MPI_DOUBLE, MPI_MAX, 0, params->comm);
  if (params->mype == 0 && max>0)
    printf(" setupTmr:  %10.4lf\n", max);

  current = gk_getwctimer(params->compTmr);
  MPI_Reduce(&current, &max, 1, MPI_DOUBLE, MPI_MAX, 0, params->comm);
  if (params->mype == 0 && max>0)
    printf("  compTmr:  %10.4lf\n", max);

  current = gk_getwctimer(params->commTmr);
  MPI_Reduce(&current, &max, 1, MPI_DOUBLE, MPI_MAX, 0, params->comm);
  if (params->mype == 0 && max>0)
    printf("  commTmr:  %10.4lf\n", max);

  current = gk_getwctimer(params->totalTmr);
  MPI_Reduce(&current, &max, 1, MPI_DOUBLE, MPI_MAX, 0, params->comm);
  if (params->mype == 0 && max>0)
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
    }

    gk_free((void **)&rowptr, LTERM);
  }

  /* broadcast rowdist */
  MPI_Bcast(dmat->rowdist, npes+1, MPI_INT, 0, params->comm);


  /* read your portion of the data */
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
    errexit("Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
  if (gk_read(fd, dmat->mat->rowptr, sizeof(ssize_t)*(lnrows+1)) != sizeof(ssize_t)*(lnrows+1))
    errexit("Failed to read the rowptr from file %s!\n", params->filename);

  /* read the rowind */
  lnnz = dmat->mat->rowptr[lnrows]-dmat->mat->rowptr[0];
  dmat->mat->rowind = gk_imalloc(lnnz, "dmat->mat->rowind");
  fpos = sizeof(int)*2 + sizeof(ssize_t)*(gnrows+1) + sizeof(int)*dmat->mat->rowptr[0];
  if (lseek64(fd, fpos, SEEK_SET) == -1)
    errexit("Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
  if (gk_read(fd, dmat->mat->rowind, sizeof(int)*lnnz) != sizeof(int)*lnnz)
    errexit("Failed to read the rowind from file %s [%zd]!\n", params->filename, sizeof(int)*lnnz);


#ifdef XXX
  /* read the rowval */
  lnnz = dmat->mat->rowptr[lnrows]-dmat->mat->rowptr[0];
  dmat->mat->rowval = gk_fmalloc(lnnz, "dmat->mat->rowval");
  fpos = 2*sizeof(int) + sizeof(ssize_t)*(gnrows+1) + sizeof(int)*gnnz + sizeof(float)*dmat->mat->rowptr[0];
  if (lseek64(fd, fpos, SEEK_SET) == -1)
    gk_errexit(SIGERR, "Failed to lseek for %s. error: %s!\n", params->filename, strerror(errno));
  if (gk_read(fd, dmat->mat->rowval, sizeof(float)*lnnz) != (sizeof(float)*lnnz))
    gk_errexit(SIGERR, "Failed to read the rowval from file %s [%zd]!\n", 
        params->filename, sizeof(float)*lnnz);
#endif

  /* localize adjust rowptr */
  for (i=lnrows; i>0; i--)
    dmat->mat->rowptr[i] -= dmat->mat->rowptr[0];
  dmat->mat->rowptr[0] = 0;

  close(fd);

  printf("[%3d] dmat->gnrows/lnrows: %d/%d, dmat->gncols/lncols: %d/%d, "
      "dmat->gnnz/lnnz: %zu/%zu [ts: %d]\n", 
      mype, 
      dmat->gnrows, dmat->mat->nrows, 
      dmat->gncols, dmat->mat->ncols, 
      dmat->gnnz, dmat->mat->rowptr[dmat->mat->nrows],(int)time(NULL));

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
  gk_startwctimer(params->sort1Tmr);
  for (i=0; i<nrows; i++) {
    if (rowptr[i+1]-rowptr[i] > 1)
      gk_isorti(rowptr[i+1]-rowptr[i], rowind+rowptr[i]);
  }
  gk_stopwctimer(params->sort1Tmr);


  /* allocate memory for the sendinfo */
  dmat->scounts = gk_imalloc(npes, "scounts");
  dmat->sdispls = gk_imalloc(npes+1, "sdispls");
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

  printf("[%3d] nsend: %d, nrecv: %d [ts: %d]\n", mype, dmat->nsend, 
      dmat->nrecv,(int)time(NULL));

  return;
}


/**************************************************************************/
/*! This function deallocates all the memory that was used */
/**************************************************************************/
void CleanupData(params_t *params, dcsr_t *dmat)
{

  printf("[pe %2d] Clean.1 [ts: %d]\n",params->mype,(int)time(NULL));
  gk_csr_Free(&(dmat->mat));
  gk_free((void **)&dmat->rowdist,
                   &dmat->scounts, &dmat->sdispls, 
                   &dmat->rcounts, &dmat->rdispls, &dmat->rinds, 
                   &dmat, LTERM);

  printf("[pe %2d] Clean.2 [ts: %d]\n",params->mype,(int)time(NULL));
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
  double *pr, *prnew, *prrecv;
  double wgt, lambda=.2, rprob, lrmsd, grmsd, lsinks, gsinks=1.0;

  nsend = dmat->nsend;
  nrecv = dmat->nrecv;
  rinds = dmat->rinds;

  nrows  = dmat->mat->nrows;
  rowptr = dmat->mat->rowptr;
  rowind = dmat->mat->rowind;

  rprob = 1.0/dmat->gnrows;
  pr    = gk_dsmalloc(nrows, rprob, "pr");

  /* get into the PR iteration */
  for (iter=0; iter<params->niters; iter++) {
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

    /* free pr if you will not need it */
    if (!(iter%10 == 0 || iter == params->niters-1)) 
      gk_free((void **)&pr, LTERM);

    gk_startwctimer(params->commTmr);
    /* get the overall gsinks across all processors */
    if (gsinks > 0)
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

    if (iter%10 == 0 || iter == params->niters-1) {
      /* compute the difference */
      for (lrmsd=0.0, i=0; i<nrows; i++) 
        lrmsd += (pr[i]-prnew[i])*(pr[i]-prnew[i]);

      gk_free((void **)&pr, LTERM);

      /* get the global rmsd across all processors */
      gk_startwctimer(params->commTmr);
      MPI_Allreduce(&lrmsd, &grmsd, 1, MPI_DOUBLE, MPI_SUM, params->comm);
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
