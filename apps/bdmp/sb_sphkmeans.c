/*!
\file
\brief A parallel spherical k-means program
\date Started 4/20/2013
\author George
*/

#define _SVID_SOURCE

#if 0
int main() {}
#else
#include <GKlib.h>
#include <bdmpi.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/mman.h>
#include <unistd.h>
#include <malloc.h>
#include <pthread.h>



/**************************************************************************/
/* data structures */
/**************************************************************************/
typedef struct {
  int npes, mype, nthreads;
  BDMPI_Comm comm;
  int nclusters;
  int ntrials, niters;
  char *filename;
  int mlock;

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
void ComputeClusteringStatistics(params_t *params, gk_csr_t *mat, int *cpart);
void printInCoreInfo(char *msg, int mype, gk_csr_t *mat);


/**************************************************************************/
/**************************************************************************/
int main(int argc, char **argv)
{
  params_t *params;
  gk_csr_t *mat;
  int *cvec;
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
      fprintf(stderr, "Usage: %s filename #clusters #trials #iters #threads mlock[0;1]\n", argv[0]);

    goto DONE;
  }

  params->filename  = strdup(argv[1]);
  params->nclusters = atoi(argv[2]);
  params->ntrials   = atoi(argv[3]);
  params->niters    = atoi(argv[4]);
  params->nthreads  = atoi(argv[5]);
  params->mlock     = atoi(argv[6]);

  printf("[%3d] nclusters: %d, ntrials: %d, niters: %d, nthreads: %d, mlock: %d\n",
      params->mype, params->nclusters, params->ntrials, params->niters, params->nthreads,
      params->mlock);

  omp_set_num_threads(params->nthreads);

  gk_clearwctimer(params->totalTmr);
  gk_clearwctimer(params->compTmr);
  gk_clearwctimer(params->commTmr);

  BDMPI_Barrier(params->comm);
  BDMPI_Barrier(params->comm);
  gk_startwctimer(params->totalTmr);

  mat = LoadData(params);

  PreprocessData(params, mat);

  srand(params->mype+101);

  printf("[%03d] timestamp01: %zu\n", params->mype, (size_t)time(NULL));
  cvec = ClusterData(params, mat);
  printf("[%03d] timestamp02: %zu\n", params->mype, (size_t)time(NULL));

  //WriteClustering(params, mat, cvec);

  //sbsaveall();
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

DONE:
  BDMPI_Finalize();

  return EXIT_SUCCESS;
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
  int lnrows, flag;
  size_t lnnz;
  gk_csr_t *mat=NULL;
  BDMPI_Status status;

  if (!gk_fexists(params->filename)) 
    gk_errexit(SIGERR, "File %s does not exist!\n", params->filename);
  
  mat = gk_csr_Read(params->filename, GK_CSR_FMT_BINROW, 1, 0);
  lnrows = mat->nrows;
  lnnz = mat->rowptr[mat->nrows];

  params->nrows = lnrows*npes;
  params->ncols = mat->ncols;

  printf("[%3d] nrows: %d, ncols: %d, tnnz: %zu, lnrows: %d, lnnz: %zu [ts: %d]\n", 
      mype, params->nrows, params->ncols, npes*lnnz, lnrows, lnnz, (int)time(NULL));

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
  BDMPI_Status status;
  FILE *fpout;
  char outfile[1024];

  sprintf(outfile, "%s.part.%d", params->filename, params->nclusters);

  if (mype == 0) {
    fpout = gk_fopen(outfile, "w", "outfile");
    for (i=0; i<mat->nrows; i++)
      fprintf(fpout, "%d\n", cvec[i]);
    gk_fclose(fpout);
    if (mype+1 < npes)
      BDMPI_Send(&dummy, 1, BDMPI_INT, mype+1, 1, params->comm);
  }
  else {
    BDMPI_Recv(&dummy, 1, BDMPI_INT, mype-1, 1, params->comm, &status);
    fpout = gk_fopen(outfile, "a", "outfile");
    for (i=0; i<mat->nrows; i++)
      fprintf(fpout, "%d\n", cvec[i]);
    gk_fclose(fpout);
    if (mype+1 < npes)
      BDMPI_Send(&mype, 1, BDMPI_INT, mype+1, 1, params->comm);
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
  float *rowval, *centers=NULL;
  float dnorms[params->nclusters], crval, bcrval;
  vlp_ii_t lmaxloc, gmaxloc;

  nclusters = params->nclusters;

  nrows  = mat->nrows;
  ncols  = mat->ncols;
  rowptr = mat->rowptr;
  rowind = mat->rowind;
  rowval = mat->rowval;

  bcpart = gk_imalloc(nrows, "bcpart");

  /* perform a number of random trials */
  for (bcrval=0.0, trial=0; trial<params->ntrials; trial++) {
    /* select the initial cluster seeds */
    lmaxloc.val = rand();
    lmaxloc.loc = mype;

    BDMPI_Allreduce(&lmaxloc, &gmaxloc, 1, BDMPI_2INT, BDMPI_MAXLOC, params->comm);

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
      if (mype == gmaxloc.loc) {
        printf("Working on trial: %zu, iter: %zu\n", trial, iter);
      }
      else {
        centers = gk_fmalloc(nclusters*ncols, "centers");
      }

      BDMPI_Bcast(centers, ncols*nclusters, BDMPI_FLOAT, gmaxloc.loc, params->comm);
      printf("[%03d]%04zu.%04zu.0 ts: %d\n", mype, trial, iter, (int)time(NULL));

      if (params->mlock)
        GKWARN(BDMPI_mlockall(MCL_CURRENT) == 0);

      printf("[%03d]%04zu.%04zu.1 ts: %d\n", mype, trial, iter, (int)time(NULL));

      /* assign each local row to the closest cluster */
      gk_startwctimer(params->compTmr);

      //BDMPI_sbload(rowptr);
      //BDMPI_sbload(rowind);
      //BDMPI_sbload(rowval);
      //BDMPI_sbload(centers);
      //BDMPI_sbload(cpart);
      lnmoves = 0;
#pragma omp parallel default(none),\
                     shared(nrows, nclusters, rowptr, rowind, rowval, centers, cpart),\
                     private(i, j, k, offset),\
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
      gk_stopwctimer(params->compTmr);

      if (params->mlock)
        GKWARN(BDMPI_munlockall() == 0);

      printf("[%03d]%04zu.%04zu.2 ts: %d\n", mype, trial, iter, (int)time(NULL));

      /* compute the new global centers */
      BDMPI_Reduce(centers, centers, nclusters*ncols, BDMPI_FLOAT, BDMPI_SUM, 
          gmaxloc.loc, params->comm);

      if (mype == gmaxloc.loc) {
        if (params->mlock)
          GKWARN(BDMPI_mlock(centers, ncols*nclusters*sizeof(float)) == 0);

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
            
        if (params->mlock)
          GKWARN(BDMPI_munlock(centers, ncols*nclusters*sizeof(float)) == 0);

        //printf("trial: %2zd; iter: %3zd; crval: %.8e; bcrval: %.8e\n", trial, iter, crval, bcrval);
      }
      else {
        gk_free((void **)&centers, LTERM);
      }

      /* see if you are done refining */
      if (iter > 0) {
        BDMPI_Allreduce(&lnmoves, &gnmoves, 1, BDMPI_INT, BDMPI_SUM, params->comm);
        if (gnmoves == 0)
          break;
      }
    }

    BDMPI_Bcast(&crval, 1, BDMPI_FLOAT, gmaxloc.loc, params->comm);
    if (crval > bcrval) {
      gk_SWAP(cpart, bcpart, tptr);
      bcrval = crval;
    }

    gk_free((void **)&cpart, LTERM);

    if (mype == gmaxloc.loc)
      printf("[%3zu:%3zu] gnmoves: %8d; crval: %8.4e; bcrval: %8.4e [ts: %d]\n",
             trial, iter, gnmoves, crval, bcrval, (int)time(NULL));
  }

  if (centers != NULL) 
    gk_free((void **)&centers, LTERM);

  return bcpart;
}

  
/**************************************************************************/
/*! This function prints final statistics for the clustering solution. */
/**************************************************************************/
void ComputeClusteringStatistics(params_t *params, gk_csr_t *mat, int *cpart)
{
  int npes=params->npes, mype=params->mype, nclusters=params->nclusters;
  size_t i, j, k, offset;
  int nrows, ncols, *rowind, *pwgts;
  ssize_t *rowptr;
  float *rowval, *centers, *dnorms;
  float crval, tcrval;

  nrows  = mat->nrows;
  ncols  = mat->ncols;
  rowptr = mat->rowptr;
  rowind = mat->rowind;
  rowval = mat->rowval;

  centers = gk_fsmalloc(nclusters*ncols, 0.0, "centers");
  pwgts   = gk_ismalloc(nclusters, 0, "pwgts");

  /* compute the local centers and local partition weights */
  gk_fset(nclusters*ncols, 0.0, centers);
  for (i=0; i<nrows; i++) {
    k = cpart[i];
    pwgts[k]++;
    for (j=rowptr[i]; j<rowptr[i+1]; j++) 
      centers[rowind[j]*nclusters+k] += rowval[j];
  }

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


/**************************************************************************/
/*! Print residency information about the various elements of the matrix */
/**************************************************************************/
void printInCoreInfo(char *msg, int mype, gk_csr_t *mat)
{
  size_t i, vlen, pagesize, addr, size, offset, counts1[2], counts2[2], counts3[2];
  unsigned char *vec;
  struct rusage usage;
  char *ptr;

  pagesize = sysconf(_SC_PAGESIZE);
  vlen = mat->rowptr[mat->nrows]*sizeof(int)/pagesize + 10;
  vec = (unsigned char *)gk_malloc(sizeof(unsigned char *)*vlen, "vec");

  counts1[0] = counts1[1] = 0;
  counts2[0] = counts2[1] = 0;
  counts3[0] = counts3[1] = 0;

  ptr = (char *)mat->rowptr;
  size = sizeof(ssize_t)*(mat->nrows+1);
  addr = (size_t)ptr;
  offset = addr%pagesize;
  ptr -= offset;
  size += offset;
  if (mincore(ptr, size, vec) == -1)
    printf("mincore error for rowptr: %s [%zu %zu %zu]\n", strerror(errno), 
        pagesize, addr, addr%pagesize);
  else {
    vlen = (size+pagesize-1)/pagesize;
    for (i=0; i<vlen; i++)
      counts1[vec[i]&1]++;
  }

  ptr = (char *)mat->rowind;
  size = sizeof(int)*(mat->rowptr[mat->nrows]);
  addr = (size_t)ptr;
  offset = addr%pagesize;
  ptr -= offset;
  size += offset;
  if (mincore(ptr, size, vec) == -1)
    printf("mincore error for rowind: %s [%zu %zu %zu]\n", strerror(errno), 
        pagesize, addr, addr%pagesize);
  else {
    vlen = (size+pagesize-1)/pagesize;
    for (i=0; i<vlen; i++)
      counts2[vec[i]&1]++;
  }

  ptr = (char *)mat->rowval;
  size = sizeof(int)*(mat->rowptr[mat->nrows]);
  addr = (size_t)ptr;
  offset = addr%pagesize;
  ptr -= offset;
  size += offset;
  if (mincore(ptr, size, vec) == -1)
    printf("mincore error for rowval: %s [%zu %zu %zu]\n", strerror(errno), 
        pagesize, addr, addr%pagesize);
  else {
    vlen = (size+pagesize-1)/pagesize;
    for (i=0; i<vlen; i++)
      counts3[vec[i]&1]++;
  }

  gk_free((void **)&vec, LTERM);

  if (getrusage(RUSAGE_SELF, &usage) == -1) {
    printf("getrusage error: %s\n", strerror(errno));
  }
  else {
    printf("[%03d]%s [%5zu %5zu] [%5zu %5zu] [%5zu %5zu] mrss: %ld minf: %ld majf: %ld\n", 
        mype, msg, 
        counts1[0], counts1[1],
        counts2[0], counts2[1],
        counts3[0], counts3[1],
        usage.ru_maxrss, usage.ru_minflt, usage.ru_majflt);
  }
}
#endif
