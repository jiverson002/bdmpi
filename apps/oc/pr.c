/*!
\file
\brief A parallel page-rank program using 1D distribution of the graph
\date Started 5/29/2013
\author George
*/


#define _LARGEFILE64_SOURCE
#include <GKlib.h>
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
  int npes;
  char *filename;
  int niters;

  /* oc filenames */
  char **mfiles;

  /* timers */
  double totalTmr;
  double setupTmr;
  double compTmr;
  double commTmr;

} params_t;

/* distributed CSR */
typedef struct {
  /* the total number of rows/non-zeros and their overall distribution */
  int gnrows, gncols, lnrows;
  size_t gnnz;
  int *rowdist;

  gk_csr_t *mat;
} dcsr_t;


/**************************************************************************/
/* prototypes */
/**************************************************************************/
dcsr_t *LoadData(params_t *params);
void WritePR(params_t *params, dcsr_t *dmat, double *prvec);
void CleanupData(params_t *params, dcsr_t *dmat);
double *ComputePR(params_t *params, dcsr_t *dmat);


/**************************************************************************/
/**************************************************************************/
int main(int argc, char **argv)
{
  int pe;
  params_t *params;
  dcsr_t *dmat;
  double *prvec;
  double max, current;
  char string[8192];

  setbuf(stdout, NULL);
  setbuf(stderr, NULL);

  params = (params_t *)gk_malloc(sizeof(params_t), "params");
  memset(params, 0, sizeof(params_t));

  if (argc != 4) {
    fprintf(stderr, "Usage: %s npes filename niters\n", argv[0]);
    return EXIT_FAILURE;
  }

  params->npes     = atoi(argv[1]);
  params->filename = gk_strdup(argv[2]);
  params->niters   = atoi(argv[3]);

  params->mfiles = (char **)gk_malloc(params->npes*sizeof(char *), "mfiles");
  for (pe=0; pe<params->npes; pe++) {
    sprintf(string, "m%d.%d", (int)getpid(), pe);
    params->mfiles[pe] = gk_strdup(string);
  }

  gk_clearwctimer(params->totalTmr);
  gk_clearwctimer(params->setupTmr);
  gk_clearwctimer(params->compTmr);
  gk_clearwctimer(params->commTmr);

  gk_startwctimer(params->totalTmr);

  dmat = LoadData(params);

  gk_startwctimer(params->compTmr);
  prvec = ComputePR(params, dmat);
  gk_stopwctimer(params->compTmr);

  //WritePR(params, dmat, prvec);

  CleanupData(params, dmat);

  gk_stopwctimer(params->totalTmr);

  /* print timing stats */
  max = gk_getwctimer(params->setupTmr);
  printf(" setupTmr:  %10.4lf\n", max);

  max = gk_getwctimer(params->compTmr);
  printf("  compTmr:  %10.4lf\n", max);

  max = gk_getwctimer(params->commTmr);
  printf("  commTmr:  %10.4lf\n", max);

  max = gk_getwctimer(params->totalTmr);
  printf(" totalTmr:  %10.4lf\n", max);

  return EXIT_SUCCESS;
}


/**************************************************************************/
/*! Reads a sparse matrix in binary CSR format, one process at a time. 
    \returns the local portion of the matrix.
*/
/**************************************************************************/
dcsr_t *LoadData(params_t *params)
{
  int pe, npes=params->npes;
  size_t i, gnnz, lnnz;
  int fd, gnrows, gncols, lnrows;
  dcsr_t *dmat=NULL;
  off64_t fpos;
  ssize_t *rowptr;

  if (!gk_fexists(params->filename)) 
    errexit("File %s does not exist!\n", params->filename);

  dmat = (dcsr_t *)gk_malloc(sizeof(dcsr_t), "dmat");
  memset(dmat, 0, sizeof(dcsr_t));

  dmat->rowdist = gk_imalloc(npes+1, "rowdist");

  /* determine the rowdist array so that it balances the lnnz's */
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

  dmat->gnrows = gnrows;
  dmat->gncols = gncols;
  dmat->gnnz   = rowptr[gnrows];

  /* populate the rowdist */
  dmat->rowdist[0] = 0;
  for (i=0, pe=0; pe<npes; pe++) {
    lnnz = rowptr[i] + (rowptr[gnrows] - rowptr[i] + npes - pe - 1)/(npes-pe);
    for (; i<gnrows; i++) {
      if (rowptr[i] >= lnnz)
        break;
    }
    dmat->rowdist[pe+1] = i;
  }

  gk_free((void **)&rowptr, LTERM);


  for (pe=0; pe<npes; pe++) {
    if ((fd = open(params->filename, O_RDONLY)) == -1)
      errexit("Failed opeing the file %s. [%s]\n", params->filename, strerror(errno));

    lnrows = dmat->rowdist[pe+1]-dmat->rowdist[pe];

    dmat->mat = gk_csr_Create();
    dmat->mat->nrows = lnrows;
    dmat->mat->ncols = gncols;

    /* read the rowptr */
    dmat->mat->rowptr = gk_zmalloc(lnrows+1, "dmat->mat->rowptr");
    fpos = 2*sizeof(int) + dmat->rowdist[pe]*sizeof(ssize_t);
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
      errexit("Failed to read the rowind from file %s [%zd]!\n",
          params->filename, sizeof(int)*lnnz);

    /* localize adjust rowptr */
    for (i=lnrows; i>0; i--)
      dmat->mat->rowptr[i] -= dmat->mat->rowptr[0];
    dmat->mat->rowptr[0] = 0;

    close(fd);

    printf("[%3d] dmat->gnrows/lnrows: %d/%d, dmat->gncols/lncols: %d/%d, dmat->gnnz/lnnz: %zu/%zu\n", 
        pe, 
        dmat->gnrows, dmat->mat->nrows, 
        dmat->gncols, dmat->mat->ncols, 
        dmat->gnnz, dmat->mat->rowptr[dmat->mat->nrows]);

    gk_csr_Write(dmat->mat, params->mfiles[pe], GK_CSR_FMT_BINROW, 0, 0);
    gk_csr_Free(&(dmat->mat));
  }

  return dmat;
}


/**************************************************************************/
/*! Writes the page-rank vector. It just let each process write its portion
    to the file in a round-robin fashion.
*/
/**************************************************************************/
void WritePR(params_t *params, dcsr_t *dmat, double *prvec)
{
  int npes=params->npes;
  size_t i;
  FILE *fpout;
  char outfile[1024];

  sprintf(outfile, "%s.prvec", params->filename);

  fpout = gk_fopen(outfile, "w", "outfile");
  for (i=0; i<dmat->gnrows; i++)
    fprintf(fpout, "%.8le\n", prvec[i]);
  gk_fclose(fpout);
}


/**************************************************************************/
/*! This function deallocates all the memory that was used */
/**************************************************************************/
void CleanupData(params_t *params, dcsr_t *dmat)
{
  int pe;

  gk_free((void **)&dmat->rowdist, &dmat, LTERM);

  for (pe=0; pe<params->npes; pe++)
    unlink(params->mfiles[pe]);

  return;
}


/**************************************************************************/
/*! This function computes the page-rank scores using a push approach */
/**************************************************************************/
double *ComputePR(params_t *params, dcsr_t *dmat)
{
  int npes=params->npes, pe;
  size_t iter, i, j, firstrow, nrows, gnrows;
  int *rowind;
  ssize_t *rowptr;
  double *pr, *prnew;
  double wgt, lambda=.2, rprob, rmsd, gsinks;
  gk_csr_t *mat;

  gnrows  = dmat->gnrows;

  rprob = 1.0/gnrows;
  pr    = gk_dsmalloc(gnrows, rprob, "pr");

  /* get into the PR iteration */
  for (iter=0; iter<params->niters; iter++) {
    prnew = gk_dsmalloc(gnrows, 0.0, "prnew");

    for (gsinks=0.0, pe=0; pe<npes; pe++) {
      firstrow = dmat->rowdist[pe];

      mat = gk_csr_Read(params->mfiles[pe], GK_CSR_FMT_BINROW, 0, 0); 
      printf("[%2d] %d %zd \n", pe, mat->nrows, mat->rowptr[mat->nrows]);
      nrows  = mat->nrows;
      rowptr = mat->rowptr;
      rowind = mat->rowind;

      /* push random-walk scores to the outlinks */
      for (i=0; i<nrows; i++) {
        if (rowptr[i+1]-rowptr[i] == 0) {
          gsinks += pr[firstrow+i];
          continue;
        }
        else {
          wgt = pr[firstrow+i]/(rowptr[i+1]-rowptr[i]);
          for (j=rowptr[i]; j<rowptr[i+1]; j++) 
            prnew[rowind[j]] += wgt;
        }
      }
      gk_csr_Free(&mat);
    }

    /* apply the restart condition */
    for (i=0; i<gnrows; i++)
      prnew[i] = lambda*rprob + (1.0-lambda)*(gsinks*rprob+prnew[i]);


    if (1 || iter%10 == 0 || iter == params->niters) {
      /* compute the difference */
      for (rmsd=0.0, i=0; i<gnrows; i++) 
        rmsd += (pr[i]-prnew[i])*(pr[i]-prnew[i]);
      rmsd = sqrt(rmsd);

      printf("Iter: %5zu, rmsd: %.6le [gsinks: %.6le, rprob: %.6le]\n", 
          iter, rmsd, gsinks, rprob);
    }

    for (i=0; i<gnrows; i++) 
      pr[i] = prnew[i];

    gk_free((void **)&prnew, LTERM);
  }

  return pr;
}

