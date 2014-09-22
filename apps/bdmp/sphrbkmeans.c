/*!
\file
\brief A parallel sphereical k-means program
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
  int npes, mype;
  BDMPI_Comm comm;
  int nclusters;
  char *filename;

  /* the total number of rows and their overall distribution */
  int nrows, ncols;
  int *rowdist;

  /* timers */
  double totalTmr;
  double compTmr;
  double commTmr;
} params_t;


/**************************************************************************/
/* prototypes */
/**************************************************************************/
gk_csr_t *LoadData(params_t *params);
void WriteClustering(params_t *params, gk_csr_t *mat, int *cvec);
void PreprocessData(params_t *params, gk_csr_t *mat);
int *ClusterData(params_t *params, gk_csr_t *mat);
float BisectCluster(params_t *params, gk_csr_t *mat, int nclusters, int *cptr, 
         int *cind, int cnum);
float ComputeClusteringStatistics(params_t *params, gk_csr_t *mat, int *cptr, 
         int *cind);


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

  if (argc != 3) {
    if (params->mype == 0)
      fprintf(stderr, "Usage: %s filename #clusters\n", argv[0]);

    BDMPI_Finalize();
    return EXIT_FAILURE;
  }

  params->filename  = strdup(argv[1]);
  params->nclusters = strtol(argv[2], NULL, 10);

  gk_clearwctimer(params->totalTmr);
  gk_clearwctimer(params->compTmr);
  gk_clearwctimer(params->commTmr);

  gk_startwctimer(params->totalTmr);
  mat = LoadData(params);

  PreprocessData(params, mat);

  srand(params->mype+101);

  cvec = ClusterData(params, mat);

  WriteClustering(params, mat, cvec);
  gk_free((void **)&cvec, LTERM);

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
gk_csr_t *LoadData(params_t *params)
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
  mat->rowptr = rowptr;
  mat->rowind = rowind;
  mat->rowval = rowval;

  fpos = ftell(fpin);
  gk_fclose(fpin);
  gk_free((void **)&line, LTERM);

  if (mype != npes-1)
    BDMPI_Send(&fpos, sizeof(long int), BDMPI_BYTE, mype+1, 1, params->comm);

  BDMPI_Allreduce(&ncols, &mat->ncols, 1, BDMPI_INT, BDMPI_MAX, params->comm);
  mat->ncols++;
  params->ncols = mat->ncols;
  
  if (mype == 0) {
    printf("[%3d] params->nrows: %d, params->ncols: %d, cnnz: %zu\n", 
        mype, params->nrows, params->ncols, cnnz);
  }
  //BDMPI_Barrier(params->comm);

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
  //printf("[%3d] [%.4lf] Start preprocess [%zd %zd]\n", params->mype, gk_WClockSeconds(), 
  //    mat->rowptr[1]-mat->rowptr[0], mat->rowptr[1000]-mat->rowptr[999]);
  gk_csr_Normalize(mat, GK_CSR_ROW, 2);
  //printf("[%3d] [%.4lf] End preprocess\n", params->mype, gk_WClockSeconds());
}


/**************************************************************************/
/*! This is the top-level routine for clustering.
*/
/**************************************************************************/
int *ClusterData(params_t *params, gk_csr_t *mat)
{
  size_t i, j, k, cnum;
  int *cvec, *cptr, *cind;
  int pwgts[params->nclusters];
  
  cptr = gk_imalloc(params->nclusters+1, "cptr");
  cind = gk_imalloc(mat->nrows, "cind");

  /* initialize the cptr/cind structure */
  cptr[0] = 0; 
  cptr[1] = mat->nrows;
  for (i=0; i<mat->nrows; i++)
    cind[i] = i;

  /* get into the repeated bisecting mode */
  for (k=1; k<params->nclusters; k++) {
    /* find the largest cluster */
    for (i=0; i<k; i++) 
      pwgts[i] = cptr[i+1]-cptr[i];

    BDMPI_Allreduce(pwgts, pwgts, k, BDMPI_INT, BDMPI_SUM, params->comm);

    cnum = gk_iargmax(k, pwgts, 1);

    BisectCluster(params, mat, k, cptr, cind, cnum);
  }

  ComputeClusteringStatistics(params, mat, cptr, cind);

  cvec = gk_imalloc(mat->nrows, "cvec");
  for (k=0; k<params->nclusters; k++) {
    for (j=cptr[k]; j<cptr[k+1]; j++)
      cvec[cind[j]] = k;
  }

  gk_free((void **)&cptr, &cind, LTERM);

  return cvec;
}


/**************************************************************************/
/*! This function bisects the cnum cluster.
*/
/**************************************************************************/
float BisectCluster(params_t *params, gk_csr_t *mat, int nclusters, int *cptr, 
         int *cind, int cnum)
{
  int npes=params->npes, mype=params->mype, dummy=0;
  size_t trial, iter, i, ii, j, k;
  int ncols, ncrows, *rowind, *cpart, *bcpart, *rowids;
  int myrnum, grnum;
  int npart, lnmoves, gnmoves;
  ssize_t *rowptr;
  float *rowval, *centers[2], *center, *dcenter;
  float dsim, dnorms[2], crval, bcrval;

  //printf("[%3d] Bisecting cluster %d\n", mype, cnum);

  ncols  = mat->ncols;
  rowptr = mat->rowptr;
  rowind = mat->rowind;
  rowval = mat->rowval;

  ncrows = cptr[cnum+1]-cptr[cnum];
  rowids = cind + cptr[cnum];

  dcenter    = gk_fmalloc(2*ncols, "dcenter");
  centers[0] = dcenter;
  centers[1] = dcenter+ncols;

  cpart  = gk_imalloc(ncrows, "cpart");
  bcpart = gk_imalloc(ncrows, "bcpart");


  /* perform a number of random trials */
  for (bcrval=0.0, trial=0; trial<10; trial++) {
    gk_iset(ncrows, -1, cpart);

    gk_startwctimer(params->commTmr);
    /* select the initial cluster seeds */
    if (ncrows > 10) 
      myrnum = rand();
    else
      myrnum = 0;

    BDMPI_Allreduce(&myrnum, &grnum, 1, BDMPI_INT, BDMPI_MAX, params->comm);
    if (myrnum == grnum) 
      myrnum = mype;
    else
      myrnum = npes;

    BDMPI_Allreduce(&myrnum, &grnum, 1, BDMPI_INT, BDMPI_MIN, params->comm);

    if (grnum == mype) { /* this pe will be selecting the initial centers */
      /* pick the first center */
      myrnum = RandomInRange(ncrows);
      gk_fset(ncols, 0.0, dcenter);
      i = rowids[myrnum];
      for (j=rowptr[i]; j<rowptr[i+1]; j++)
        dcenter[rowind[j]] = rowval[j];

      myrnum = (myrnum + RandomInRange(ncrows-3) + 1)%ncrows;
      i = rowids[myrnum];
      for (j=rowptr[i]; j<rowptr[i+1]; j++)
        dcenter[rowind[j]] -= rowval[j];
    }

    BDMPI_Bcast(dcenter, ncols, BDMPI_FLOAT, grnum, params->comm);
    gk_stopwctimer(params->commTmr);

    /* get into the iterative refinement */
    for (iter=0; iter<20; iter++) {
      /* assign each local row to the closest cluster */
      gk_startwctimer(params->compTmr);
      for (lnmoves=0, ii=0; ii<ncrows; ii++) {
        i = rowids[ii];
        for (dsim=0.0, j=rowptr[i]; j<rowptr[i+1]; j++) 
          dsim += rowval[j]*dcenter[rowind[j]];

        if ((npart = (dsim > 0 ? 0 : 1)) != cpart[ii])
          lnmoves++;
        cpart[ii] = npart;
      }

      /* compute the new local centers */
      gk_fset(ncols, 0.0, centers[0]);
      gk_fset(ncols, 0.0, centers[1]);
      for (ii=0; ii<ncrows; ii++) {
        center = centers[cpart[ii]];
        i = rowids[ii];
        for (j=rowptr[i]; j<rowptr[i+1]; j++) 
          center[rowind[j]] += rowval[j];
      }
      gk_stopwctimer(params->compTmr);

      /* see if you are done refining */
      if (iter > 0) {
        BDMPI_Allreduce(&lnmoves, &gnmoves, 1, BDMPI_INT, BDMPI_SUM, params->comm);
        //printf("[%3d] lnmoves: %d, gnmoves: %d\n", mype, lnmoves, gnmoves);
        if (gnmoves == 0)
          break;
      }


      /* compute the new global centers */
      gk_startwctimer(params->commTmr);
      BDMPI_Reduce(centers[0], centers[0], 2*ncols, BDMPI_FLOAT, BDMPI_SUM, 0, params->comm);

      if (mype == 0) {
        dnorms[0] = gk_fdot(ncols, centers[0], 1, centers[0], 1);
        dnorms[1] = gk_fdot(ncols, centers[1], 1, centers[1], 1);
        crval = (dnorms[0] > 0 ? sqrt(dnorms[0]) : 0.0) + 
                (dnorms[1] > 0 ? sqrt(dnorms[1]) : 0.0);

        if (dnorms[0] > 0)
          gk_fscale(ncols, 1.0/sqrt(dnorms[0]), centers[0], 1);
        if (dnorms[1] > 0)
          gk_fscale(ncols, 1.0/sqrt(dnorms[1]), centers[1], 1);
        gk_faxpy(ncols, -1.0, centers[1], 1, centers[0], 1);

        //printf("trial: %2zd; iter: %3zd; crval: %.8e; bcrval: %.8e\n", trial, iter, crval, bcrval);
      }

      BDMPI_Bcast(centers[0], ncols, BDMPI_FLOAT, 0, params->comm);
      BDMPI_Bcast(&crval, 1, BDMPI_FLOAT, 0, params->comm);
      gk_stopwctimer(params->commTmr);

      if (crval > bcrval) {
        gk_icopy(ncrows, cpart, bcpart);
        bcrval = crval;
      }
    }
  }

  /* bcpart has the clustering solution; update cptr/cind */
  for (i=nclusters; i>cnum; i--)
    cptr[i+1] = cptr[i];

  for (k=0, j=0, i=0; i<ncrows; i++) {
    if (bcpart[i] == 0)
      rowids[k++] = rowids[i];
    else
      cpart[j++] = rowids[i];
  }
  cptr[cnum+1] = cptr[cnum]+k;
  for (i=0; i<j; i++, k++)
    rowids[k] = cpart[i];
  
  gk_free((void **)&dcenter, &cpart, &bcpart, LTERM);

  if (mype == 0)
    printf("Bisecting %4d: crval: %.6e\n", cnum, bcrval);

  return bcrval;
}

  
/**************************************************************************/
/*! This function prints final statistics for the clustering solution.
*/
/**************************************************************************/
float ComputeClusteringStatistics(params_t *params, gk_csr_t *mat, int *cptr, 
         int *cind)
{
  int npes=params->npes, mype=params->mype, dummy=0;
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
    BDMPI_Reduce(lvec, gvec, ncols, BDMPI_FLOAT, BDMPI_SUM, 0, params->comm);
    BDMPI_Reduce(&lpwgt, &gpwgt, 1, BDMPI_INT, BDMPI_SUM, 0, params->comm);

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

