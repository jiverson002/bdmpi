/*!
\file
\brief A parallel samplesort program
\date Started 4/20/2013
\author George
*/


#include <GKlib.h>
#include <omp.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>



/**************************************************************************/
/*! Sorts the combined set of arrays of ints in increasing order. */
/**************************************************************************/
int samplesort(int npes, size_t nlocal)
{
  size_t p, i, j, k, nelmnts;
  int *elmnts, *mypicks, *allpicks;
  char filename[1024];
  omp_lock_t flocks[npes];

  printf("%zu Starting.\n", (size_t)time(NULL));

  /* allocate memory for allpicks to be shared by everybody */
  allpicks = gk_imalloc(npes*npes+2, "allppicks");

  for (p=0; p<npes; p++)
    omp_init_lock(&flocks[p]);


  /* sort the individual arrays and select splitters */
#pragma omp parallel default(none), \
                     shared(p, npes, allpicks), \
                     private(i, nelmnts, elmnts, mypicks, filename) 
  {
    FILE *fp;

#pragma omp for
    for (p=0; p<npes; p++) {
      /* read the data */
      sprintf(filename, "in-%zu", p);
      elmnts = gk_i32readfilebin(filename, &nelmnts);

      /* sort the local elements in increasing order */
      gk_isorti(nelmnts, elmnts);

      /* select the local npes-1 equally spaced elements and put them in allpicks */
      mypicks = allpicks + p*(npes-1);
      for (i=1; i<npes; i++) 
        mypicks[i-1] = elmnts[i*nelmnts/npes];

      /* write out the sorted data */
      fp = gk_fopen(filename, "wb", "fp-sort");
      fwrite(elmnts, sizeof(int), nelmnts, fp);
      gk_fclose(fp);

      gk_free((void **)&elmnts, LTERM);
    }
  }
  printf("%zu Done with initial sorting.\n", (size_t)time(NULL));


  /* sort all the picks */
  gk_isorti(npes*(npes-1), allpicks);

  /*
  for (i=0; i<npes*(npes-1); i++)
    printf("ap: %d\n", allpicks[i]);
  */

  /* Select the final splitters. Set the boundaries to simplify coding */
  mypicks = gk_imalloc(npes+1, "mypicks");
  for (i=1; i<npes; i++) 
    mypicks[i] = allpicks[i*(npes-1)];
  mypicks[0]    = INT_MIN;
  mypicks[npes] = INT_MAX;

  /*
  for (i=0; i<=npes; i++)
    printf("mp: %d\n", mypicks[i]);
  */


  /* go and split each local file based on the splitters */
#pragma omp parallel default(none), \
                     shared(p, npes, mypicks, flocks), \
                     private(elmnts, filename, i, j, nelmnts) 
  {
    size_t scount;
    FILE *fp;

#pragma omp for
    for (p=0; p<npes; p++) {
      /* read the data */
      sprintf(filename, "in-%zu", p);
      elmnts = gk_i32readfilebin(filename, &nelmnts);
      gk_rmpath(filename);

      /* compute the number of elements that belong to each bucket and append it. */
      //printf("For %zu\n", p);
      for (scount=0, j=0, i=0; i<nelmnts;) {
        if (elmnts[i] <= mypicks[j+1]) {
          //printf("  %d\n", elmnts[i]);
          scount++;
          i++;
        }
        else {
          /* write the data out */
          //printf("  --- %zu [scount: %zu] [%d]\n", j, scount, elmnts[i-scount]);
          if (scount > 0) {
            sprintf(filename, "out-%zu-%zu", p, j);
            omp_set_lock(&flocks[j]);
            fp = gk_fopen(filename, "ab", "fp-append");
            fwrite(elmnts+i-scount, sizeof(int), scount, fp);
            gk_fclose(fp);
            omp_unset_lock(&flocks[j]);
            scount = 0;
          }
          j++;
        }
      }
      /* write out the last block of data */
      //printf("  --- %zu [scount: %zu] [%d]\n", j, scount, elmnts[i-scount]);
      if (scount > 0) {
        sprintf(filename, "out-%zu-%zu", p, j);
        omp_set_lock(&flocks[j]);
        fp = gk_fopen(filename, "ab", "fp-append");
        fwrite(elmnts+i-scount, sizeof(int), scount, fp);
        gk_fclose(fp);
        omp_unset_lock(&flocks[j]);
      }

      gk_free((void **)&elmnts, LTERM);
    }
  }
  printf("%zu Done with splitting.\n", (size_t)time(NULL));



  /* go and merge the individual files */
#pragma omp parallel default(none), \
                     shared(p, npes, nlocal), \
                     private(i, j, nelmnts, elmnts, filename) 
  {
    size_t rcount, nmax;
    int *relmnts;
    FILE *fp;

#pragma omp for
    for (p=0; p<npes; p++) {
      rcount  = 0;
      nmax    = 5*nlocal/4;
      relmnts = gk_imalloc(nmax, "relmnts");

      for (j=0; j<npes; j++) {
        /* read the data */
        sprintf(filename, "out-%zu-%zu", j, p);
        if (gk_fexists(filename)) {
          elmnts = gk_i32readfilebin(filename, &nelmnts);
          gk_rmpath(filename);

          if (rcount+nelmnts < nmax) {
            nmax += 5*nelmnts/3;
            relmnts = gk_irealloc(relmnts, nmax, "relmnts");
          }
          memcpy(relmnts+rcount, elmnts, sizeof(int)*nelmnts);
          rcount += nelmnts;

          //printf("[%zu %zu] %zd %p %zu %zu\n", j, p, nelmnts, (void *)elmnts, nlocal, nmax);
          gk_free((void **)&elmnts, LTERM);
        }
      }

      gk_isorti(rcount, relmnts);
      sprintf(filename, "out-%zu", p);
      fp = gk_fopen(filename, "wb", "final out");
      fwrite(relmnts, sizeof(int), rcount, fp);
      gk_fclose(fp);
      gk_free((void **)&relmnts, LTERM);
    }
  }
  printf("%zu Done with merging.\n", (size_t)time(NULL));


  /* cleanup */
  for (p=0; p<npes; p++)
    omp_destroy_lock(&flocks[p]);

  gk_free((void **)&allpicks, &mypicks, LTERM);

  return 1;
}


/**************************************************************************/
/**************************************************************************/
int main(int argc, char **argv)
{
  int npes, nthreads, lastelmnt = 0;
  size_t i, p, nlocal, ntotal, nelmnts;
  int *elmnts;
  char filename[1024];
  FILE *fp;
  double tmr=0.0;

  if (argc < 4) {
    fprintf(stderr, "Usage: %s np nr nelems\n", argv[0]);
    return EXIT_FAILURE;
  }

  npes     = strtol(argv[1], NULL, 10);
  nthreads = strtol(argv[2], NULL, 10);
  nlocal   = strtol(argv[3], NULL, 10);

  printf("nlocal: %zd\n", nlocal);

  omp_set_num_threads(nthreads);

  gk_startwctimer(tmr);


#pragma omp parallel default(none), \
                     shared(p, npes, nlocal), \
                     private(i, elmnts, filename, fp)
  {
    unsigned int seed;
    int rnum;

#pragma omp for  
    for (p=0; p<npes; p++) {
      seed = p+7;

      elmnts = gk_imalloc(nlocal, "elmnts");
      for (i=0; i<nlocal; i++) 
        elmnts[i] = rand_r(&seed);

      sprintf(filename, "in-%zu", p);
      fp = gk_fopen(filename, "wb", "fpout");
      fwrite(elmnts, sizeof(int), nlocal, fp);
      gk_fclose(fp);

      gk_free((void **)&elmnts, LTERM);
    }
  }


  /* perform the sorting */
  if (samplesort(npes, nlocal) == 0) {
    printf("Samplesort returned with an error: npes: %d\n", npes);
    return EXIT_FAILURE;
  }


  /* the root will get all the data and write them to stdout */
  ntotal = 0;
  for (p=0; p<npes; p++) {
    sprintf(filename, "out-%zu", p);
    elmnts = gk_i32readfilebin(filename, &nelmnts);
    ntotal += nelmnts;

    if (p == 0) {
      lastelmnt = elmnts[nelmnts-1];
    }
    else {
      if (nelmnts > 0) {
        if (elmnts[0] < lastelmnt)
          printf("Inversion: %d %d [p: %zd]\n", lastelmnt, elmnts[0], p);
        lastelmnt = elmnts[nelmnts-1];
      }
    }

    gk_free((void **)&elmnts, LTERM);
    gk_rmpath(filename);
  }
  printf("%zu Done. [%zu]\n", (size_t)time(NULL), ntotal);

  gk_stopwctimer(tmr);

  printf("Total time: %8.3lf\n", gk_getwctimer(tmr));

  return EXIT_SUCCESS;
}



