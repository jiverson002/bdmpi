/*!
\file
\brief An out-of-core multi-threades spherical k-means program
\date Started 5/17/2013
\author George
*/


#include <GKlib.h>
#include <omp.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>


/**************************************************************************/
/* data structures */
/**************************************************************************/
typedef struct {
  int npes, nthreads;
  int nclusters;
  char *filename;

  /* the total number of rows and their overall distribution */
  int nrows, ncols;
  size_t nnz;

  /* out-of-core files */
  char **mfiles;
  char **pfiles;
  char **bfiles;

  /* timers */
  double totalTmr;
  double loadTmr;
  double prepTmr;
  double compTmr;
  double commTmr;
} params_t;


/**************************************************************************/
/* prototypes */
/**************************************************************************/
void LoadData(params_t *params);
void CleanData(params_t *params);
void WriteClustering(params_t *params);
void PreprocessData(params_t *params);
float ClusterData(params_t *params);
float ComputeClusteringStatistics(params_t *params);


/**************************************************************************/
/**************************************************************************/
int main(int argc, char **argv)
{
  params_t *params;
  int i, *cvec;
  double max, current;
  char string[8192];

  params = (params_t *)gk_malloc(sizeof(params_t), "params");
  memset(params, 0, sizeof(params_t));

  if (argc != 5) {
    printf("Usage: %s npes nthreads filename #clusters\n", argv[0]);
    return EXIT_FAILURE;
  }

  params->npes      = atoi(argv[1]);
  params->nthreads  = atoi(argv[2]);
  params->filename  = strdup(argv[3]);
  params->nclusters = atoi(argv[4]);

  params->mfiles = (char **)gk_malloc(sizeof(char *)*params->npes, "mfiles");
  params->pfiles = (char **)gk_malloc(sizeof(char *)*params->npes, "pfiles");
  params->bfiles = (char **)gk_malloc(sizeof(char *)*params->npes, "bfiles");
  for (i=0; i<params->npes; i++) {
    sprintf(string, "m%d-%03d", (int)getpid(), i);
    params->mfiles[i] = gk_strdup(string);
    sprintf(string, "p%d-%03d", (int)getpid(), i);
    params->pfiles[i] = gk_strdup(string);
    sprintf(string, "b%d-%03d", (int)getpid(), i);
    params->bfiles[i] = gk_strdup(string);
  }

  omp_set_num_threads(params->nthreads);

  gk_clearwctimer(params->totalTmr);
  gk_clearwctimer(params->compTmr);
  gk_clearwctimer(params->commTmr);
  gk_clearwctimer(params->loadTmr);
  gk_clearwctimer(params->prepTmr);

  gk_startwctimer(params->totalTmr);

  gk_startwctimer(params->loadTmr);
  LoadData(params);
  gk_stopwctimer(params->loadTmr);

  gk_startwctimer(params->prepTmr);
  PreprocessData(params);
  gk_stopwctimer(params->prepTmr);

  ClusterData(params);

  //WriteClustering(params);

  CleanData(params);

  gk_stopwctimer(params->totalTmr);

  /* print timing stats */
  printf("  loadTmr:  %10.4lf\n", gk_getwctimer(params->loadTmr));
  printf("  prepTmr:  %10.4lf\n", gk_getwctimer(params->prepTmr));
  printf("  compTmr:  %10.4lf\n", gk_getwctimer(params->compTmr));
  printf("  commTmr:  %10.4lf\n", gk_getwctimer(params->commTmr));
  printf(" totalTmr:  %10.4lf\n", gk_getwctimer(params->totalTmr));

  return EXIT_SUCCESS;
}


/**************************************************************************/
/*! Reads a sparse matrix in headerless CSR format. The same matrix is read
    by everybody in order to simulate a large file clustering.
    \returns the local portion of the matrix. */
/**************************************************************************/
void LoadData(params_t *params)
{
  int pe, npes=params->npes;
  gk_csr_t *mat=NULL;

  params->nrows = 0;
  params->ncols = 0;
  params->nnz   = 0;
  for (pe=0; pe<npes; pe++) {
    mat = gk_csr_Read(params->filename, GK_CSR_FMT_BINROW, 1, 0);
    params->nrows += mat->nrows;
    params->nnz   += mat->rowptr[mat->nrows];
    params->ncols  = gk_max(params->ncols, mat->ncols);

    gk_csr_Write(mat, params->mfiles[pe], GK_CSR_FMT_BINROW, 1, 0);
    gk_csr_Free(&mat);
  }

  printf("params->nrows: %d, params->ncols: %d, params->nnz: %zu\n", 
        params->nrows, params->ncols, params->nnz);

}


/**************************************************************************/
/*! Writes a clustering vector. It just let each process write its portion
    to the file in a round-robin fashion. */
/**************************************************************************/
void WriteClustering(params_t *params)
{
  int pe, npes=params->npes;
  size_t i, nrows;
  int32_t *cvec;
  FILE *fpout;
  char outfile[1024];

  sprintf(outfile, "%s.part.%d", params->filename, params->nclusters);
  fpout = gk_fopen(outfile, "w", "outfile");

  for (pe=0; pe<npes; pe++) {
    cvec = gk_i32readfilebin(params->bfiles[pe], &nrows);

    for (i=0; i<nrows; i++)
      fprintf(fpout, "%d\n", cvec[i]);

    gk_free((void **)&cvec, LTERM);
  }

  gk_fclose(fpout);

}


/**************************************************************************/
/*! Removes any intermediate files. */
/**************************************************************************/
void CleanData(params_t *params)
{
  int pe, npes=params->npes;

  for (pe=0; pe<npes; pe++) {
    unlink(params->mfiles[pe]);
    unlink(params->pfiles[pe]);
    unlink(params->bfiles[pe]);
  }

}


/**************************************************************************/
/*! This function performs various pre-processing steps on the matrix. */
/**************************************************************************/
void PreprocessData(params_t *params)
{
  int pe, npes=params->npes;
  gk_csr_t *mat=NULL;

  for (pe=0; pe<npes; pe++) {
    mat = gk_csr_Read(params->mfiles[pe], GK_CSR_FMT_BINROW, 1, 0);

    gk_csr_Normalize(mat, GK_CSR_ROW, 2);

    gk_csr_Write(mat, params->mfiles[pe], GK_CSR_FMT_BINROW, 1, 0);
    gk_csr_Free(&mat);
  }
}


/**************************************************************************/
/*! This function computes the k-way clustering solution. */
/**************************************************************************/
float ClusterData(params_t *params)
{
  int pe, npes=params->npes;
  int nrows, ncols, nclusters, myrnum, lnmoves, gnmoves;
  size_t trial, iter, i, j, k, offset, nelmnts;
  int *rowind, *cpart;
  ssize_t *rowptr;
  float *rowval, *centers, *ncenters;
  float dnorms[params->nclusters], crval, bcrval;
  gk_csr_t *mat;
  char cmd[8192];

  nclusters = params->nclusters;
  ncols     = params->ncols;

  centers  = gk_fmalloc(nclusters*ncols, "centers");
  ncenters = gk_fmalloc(nclusters*ncols, "ncenters");


  /* perform a number of random trials */
  for (bcrval=0.0, trial=0; trial<1; trial++) {
    /* reset centroids */
    gk_fset(nclusters*ncols, 0.0, centers);

    /* select the pe that will select the cluster seeds */
    pe = RandomInRange(npes);

    gk_startwctimer(params->commTmr);
    mat = gk_csr_Read(params->mfiles[pe], GK_CSR_FMT_BINROW, 1, 0);
    nrows  = mat->nrows;
    rowptr = mat->rowptr;
    rowind = mat->rowind;
    rowval = mat->rowval;
    gk_stopwctimer(params->commTmr);

    /* pick the centers */
    myrnum = RandomInRange(nrows/nclusters);
    for (k=0; k<nclusters; k++) {
      i = ((k+1)*myrnum)%nrows;
      for (j=rowptr[i]; j<rowptr[i+1]; j++)
        centers[rowind[j]*nclusters+k] = rowval[j];
    }

    gk_csr_Free(&mat);


    /* get into the iterative refinement */
    for (crval=0.0, iter=0; iter<5; iter++) {
      /* clear the centers for the next iteration */
      gk_fset(nclusters*ncols, 0.0, ncenters);

      for (gnmoves=0, pe=0; pe<npes; pe++) {
        gk_startwctimer(params->commTmr);
        mat = gk_csr_Read(params->mfiles[pe], GK_CSR_FMT_BINROW, 1, 0);
        nrows  = mat->nrows;
        rowptr = mat->rowptr;
        rowind = mat->rowind;
        rowval = mat->rowval;

        if (iter == 0) {
          cpart = gk_ismalloc(nrows, -1, "cpart");
        }
        else {
          cpart = gk_i32readfilebin(params->pfiles[pe], &nelmnts);
          if (nelmnts != nrows) {
            printf("Read partition file did not match the #rows of the matrix.\n");
            exit(0);
          }
        }
        gk_stopwctimer(params->commTmr);


        /* assign each local row to the closest cluster */
        gk_startwctimer(params->compTmr);
        lnmoves = 0;
#pragma omp parallel default(none), \
                     shared(nrows, nclusters, rowptr, rowind, rowval, centers, cpart),\
                     private(i, j, k, offset),\
                     reduction(+:lnmoves)
        {
          float sims[nclusters];

#pragma omp for
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
        gnmoves += lnmoves;

        /* update the new centers */
        for (i=0; i<nrows; i++) {
          k = cpart[i];
          for (j=rowptr[i]; j<rowptr[i+1]; j++) 
            ncenters[rowind[j]*nclusters+k] += rowval[j];
        }
        gk_stopwctimer(params->compTmr);

        /* save the new solution */
        gk_startwctimer(params->commTmr);
        gk_i32writefilebin(params->pfiles[pe], nrows, cpart);

        gk_csr_Free(&mat);
        gk_free((void **)&cpart, LTERM);
        gk_stopwctimer(params->commTmr);
      }

      /* see if you are done refining */
      if (iter > 0 && gnmoves == 0)
        break;


      /* compute the quality of the solution and normalize the centroids */
      gk_startwctimer(params->compTmr);
      gk_fcopy(nclusters*ncols, ncenters, centers);
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
      gk_stopwctimer(params->compTmr);

      if (crval > bcrval) {
        bcrval = crval;
        for (pe=0; pe<npes; pe++) {
          sprintf(cmd, "cp %s %s", params->pfiles[pe], params->bfiles[pe]);
          if (system(cmd) == -1) 
            gk_errexit(SIGERR, "Failed on system(%s)\n", cmd);
        }
      }

      printf("[%3zu:%3zu] gnmoves: %6d; crval: %8.4e; bcrval: %8.4e\n", 
          trial, iter, gnmoves, crval, bcrval);
    }
  }

  gk_free((void **)&centers, &ncenters, LTERM);


  /* compute the clustering statistics */
  //ComputeClusteringStatistics(params);

  return bcrval;
}

  
/**************************************************************************/
/*! This function prints final statistics for the clustering solution. */
/**************************************************************************/
float ComputeClusteringStatistics(params_t *params)
{
  int pe, npes=params->npes, nclusters=params->nclusters;
  size_t i, ii, j, k, offset, nelmnts;
  int ncols, nrows, *rowind, *cpart;
  int pwgts[nclusters];
  ssize_t *rowptr;
  float *rowval, *centers;
  float dnorms[nclusters], crval, tcrval;
  gk_csr_t *mat;


  ncols   = params->ncols;
  centers = gk_fsmalloc(nclusters*ncols, 0.0, "centers");

  gk_iset(nclusters, 0, pwgts);

  for (pe=0; pe<npes; pe++) {
    gk_startwctimer(params->commTmr);
    mat = gk_csr_Read(params->mfiles[pe], GK_CSR_FMT_BINROW, 1, 0);
    nrows  = mat->nrows;
    rowptr = mat->rowptr;
    rowind = mat->rowind;
    rowval = mat->rowval;

    cpart = gk_i32readfilebin(params->bfiles[pe], &nelmnts);
    if (nelmnts != nrows) {
      printf("Read partition file did not match the #rows of the matrix.\n");
      exit(0);
    }
    gk_stopwctimer(params->commTmr);

    gk_startwctimer(params->compTmr);
    for (i=0; i<nrows; i++) {
      k = cpart[i];
      pwgts[k]++;
      for (j=rowptr[i]; j<rowptr[i+1]; j++) 
        centers[rowind[j]*nclusters+k] += rowval[j];
    }
    gk_stopwctimer(params->compTmr);

    gk_csr_Free(&mat);
    gk_free((void **)&cpart, LTERM);
  }

  /* compute the quality of the clusters */
  gk_startwctimer(params->compTmr);
  gk_fset(nclusters, 0.0, dnorms);
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
  gk_stopwctimer(params->compTmr);

  gk_free((void **)&centers, LTERM);

  return tcrval;
}

