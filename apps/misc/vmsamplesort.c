/*!
\file
\brief A virtual-memory based samplesort program for baseline comparison purposes.
\date Started 4/23/2013
\author George
*/


#include <GKlib.h>
#include <bdmpi.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>


/**************************************************************************/
/*! Sorts a distributed array of ints in increasing order. */
/**************************************************************************/
int samplesort(size_t *r_nlocal, int **r_elmnts, BDMPI_Comm icomm)
{
  int npes, mype, status;
  size_t i, j, k, nlocal, ntotal, nrecv, nrecv_min, nrecv_max;
  size_t *scounts, *rcounts, *sdispls, *rdispls;
  int *elmnts, *relmnts, *mypicks, *allpicks;
  BDMPI_Comm comm;

  nlocal = *r_nlocal;
  elmnts = *r_elmnts;

  /* duplicate the communicator and get my coordinates */
  BDMPI_Comm_dup(icomm, &comm);
  BDMPI_Comm_size(comm, &npes);
  BDMPI_Comm_rank(comm, &mype);

  status = (nlocal > npes ? 1 : 0);

  /* see if you can actually do the sort */
  BDMPI_Allreduce(&status, &status, 1, BDMPI_INT, BDMPI_PROD, comm);

  BDMPI_Barrier(comm);
  if (mype == 0) printf("%zu Starting.\n", (size_t)time(NULL));

  if (status == 0) {
    BDMPI_Comm_free(&comm);
    return 0;
  }
  

  /* get memory for the counts and displacements */
  scounts = (size_t *)gk_malloc(sizeof(size_t)*(npes+1), "scounts");
  rcounts = (size_t *)gk_malloc(sizeof(size_t)*(npes+1), "rcounts");
  sdispls = (size_t *)gk_malloc(sizeof(size_t)*(npes+1), "sdispls");
  rdispls = (size_t *)gk_malloc(sizeof(size_t)*(npes+1), "rdispls");

  /* get memory for the local splitters */
  mypicks  = gk_imalloc(npes+1, "mypicks");

  /* sort the local elements in increasing order */
  gk_isorti(nlocal, elmnts);

//  BDMPI_Barrier(comm);
//  if (mype == 0) printf("%zu Done with local sort.\n", (size_t)time(NULL));

  /* select the local npes-1 equally spaced elements */
  for (i=1; i<npes; i++) 
    mypicks[i-1] = elmnts[i*(nlocal/npes)];

  /* get memory for all splitters */
  if (mype == 0) 
    allpicks = gk_imalloc(npes*npes, "allppicks");

  /* gather the picks to the root processor */
  BDMPI_Gather(mypicks, npes-1, BDMPI_INT, allpicks, npes-1, BDMPI_INT, 0, comm);

  if (mype == 0) {
    /* sort all the picks */
    gk_isorti(npes*(npes-1), allpicks);

    /* Select the final splitters. Set the boundaries to simplify coding */
    for (i=1; i<npes; i++)
      mypicks[i] = allpicks[i*(npes-1)];
    mypicks[0]    = INT_MIN;
    mypicks[npes] = INT_MAX;

    gk_free((void **)&allpicks, LTERM);

    //for (i=0; i<=npes; i++)
    //  printf("s: %d\n", mypicks[i]);
  }

  /* broadcast the 2nd level picks to all the processors */
  BDMPI_Bcast(mypicks, npes+1, BDMPI_INT, 0, comm);
  if (mype == 0) printf("%zu Done with splitter bcast.\n", (size_t)time(NULL));

  /* compute the number of elements that belong to each bucket */
  for (i=0; i<npes; i++)
    scounts[i] = 0;
  for (j=0, i=0; i<nlocal;) {
    if (elmnts[i] <= mypicks[j+1]) {
      scounts[j]++;
      i++;
    }
    else {
      j++;
    }
  }
  BDMPI_Alltoall(scounts, 1, BDMPI_SIZET, rcounts, 1, BDMPI_SIZET, comm);
  if (mype == 0) printf("%zu Done with rcounts alltoall.\n", (size_t)time(NULL));

  sdispls[0] = rdispls[0] = 0;
  for (i=0; i<npes; i++) {
    sdispls[i+1] = sdispls[i] + scounts[i];
    rdispls[i+1] = rdispls[i] + rcounts[i];
  }

  /* allocate memory for sorted elements and receive them */
  nrecv   = rdispls[npes];  
  relmnts = gk_imalloc(nrecv, "relmnts");

  BDMPI_Alltoallv(elmnts, scounts, sdispls, BDMPI_INT,
                 relmnts, rcounts, rdispls, BDMPI_INT, 
                 comm);
  if (mype == 0) printf("%zu Done with data alltoall.\n", (size_t)time(NULL));

  /* do the local sort of the received elements */
  gk_isorti(nrecv, relmnts);

  //BDMPI_Barrier(comm);
  //if (mype == 0) printf("%zu Done with 2nd sort.\n", (size_t)time(NULL));

  /* hook the fields to the data that needs to be returned */
  gk_free((void **)r_elmnts, LTERM);
  *r_nlocal = nrecv;
  *r_elmnts = relmnts;

  BDMPI_Reduce(&nrecv, &nrecv_min, 1, BDMPI_SIZET, BDMPI_MIN, 0, comm);
  BDMPI_Reduce(&nrecv, &nrecv_max, 1, BDMPI_SIZET, BDMPI_MAX, 0, comm);

  if (mype == 0)
    printf("%zu min/max nrecv: %zu %zu\n", (size_t)time(NULL), nrecv_min, nrecv_max);

  //BDMPI_Barrier(comm);
  //if (mype == 0) printf("%zu Done.\n", (size_t)time(NULL));

  gk_free((void **)&mypicks, &rcounts, &scounts, &rdispls, &sdispls, LTERM);

  BDMPI_Comm_free(&comm);

  return 1;
}


/**************************************************************************/
/**************************************************************************/
int main(int argc, char **argv)
{
  int k, npes, mype, lastelmnt;
  size_t i, nlocal, nremote;
  int *elmnts;
  BDMPI_Status status;

  BDMPI_Init(&argc, &argv);
  BDMPI_Comm_size(BDMPI_COMM_WORLD, &npes);
  BDMPI_Comm_rank(BDMPI_COMM_WORLD, &mype);

  if (argc < 2) {
    if (mype == 0)
      fprintf(stderr, "Usage: %s nelems\n", argv[0]);

    BDMPI_Finalize();
    return EXIT_FAILURE;
  }

  nlocal = strtol(argv[1], NULL, 10);

  if (mype == 0) 
    printf("nlocal: %zd\n", nlocal);
  BDMPI_Barrier(BDMPI_COMM_WORLD);

  srand(5+mype);

  elmnts = gk_imalloc(nlocal, "elmnts");
  for (i=0; i<nlocal; i++) 
    elmnts[i] = rand();

  /* perform the sorting */
  if (samplesort(&nlocal, &elmnts, BDMPI_COMM_WORLD) == 0) {
    if (mype == 0)
      printf("Samplesort returned with an error: npes: %d\n", npes);
    gk_free((void **)&elmnts, LTERM);
    BDMPI_Finalize();
    return EXIT_FAILURE;
  }

  //BDMPI_Barrier(BDMPI_COMM_WORLD);
  //printf("[%3d] nlocal: %zu\n", mype, nlocal);

  /* the root will get all the data and write them to stdout */
  if (mype == 0) {
    //for (i=0; i<nlocal; i++)
    //  printf("%5d %d\n", mype, elmnts[i]);
    lastelmnt = elmnts[nlocal-1];
    gk_free((void **)&elmnts, LTERM);
  }


  for (k=1; k<npes; k++) {
    if (mype == k) {
      BDMPI_Send(&nlocal, 1, BDMPI_SIZET, 0, 1, BDMPI_COMM_WORLD);
      BDMPI_Send(elmnts, nlocal, BDMPI_INT, 0, 2, BDMPI_COMM_WORLD);
      gk_free((void **)&elmnts, LTERM);
    }
    else if (mype == 0) {
      BDMPI_Recv(&nremote, 1, BDMPI_SIZET, k, 1, BDMPI_COMM_WORLD, &status);
      elmnts = gk_imalloc(nremote, "elmnts");
      BDMPI_Recv(elmnts, nremote, BDMPI_INT, k, 2, BDMPI_COMM_WORLD, &status);
      if (elmnts[0] < lastelmnt)
        printf("Inversion: %d %d [k: %d]\n", lastelmnt, elmnts[0], k);
      //for (i=0; i<nremote; i++)
      //  printf("%5d %d\n", k, elmnts[i]);
      lastelmnt = elmnts[nremote-1];
      gk_free((void **)&elmnts, LTERM);
    }
    BDMPI_Barrier(BDMPI_COMM_WORLD);
  }

  BDMPI_Finalize();

  if (mype == 0) printf("%zu Done done!\n", (size_t)time(NULL));

  return EXIT_SUCCESS;
}



