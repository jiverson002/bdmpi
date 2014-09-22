/*!
\file
\brief A parallel samplesort program
\date Started 4/20/2013
\author George
*/


#include <GKlib.h>
#include <mpi.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>


/**************************************************************************/
/*! Sorts a distributed array of ints in increasing order. */
/**************************************************************************/
int samplesort(int *r_nlocal, int **r_elmnts, MPI_Comm icomm)
{
  int npes, mype, mystatus, status;
  int i, j, k, nlocal, ntotal, nrecv, nrecv_min, nrecv_max;
  int *scounts, *rcounts, *sdispls, *rdispls;
  int *elmnts, *relmnts, *mypicks, *allpicks;
  MPI_Comm comm;

  nlocal = *r_nlocal;
  elmnts = *r_elmnts;

  /* duplicate the communicator and get my coordinates */
  MPI_Comm_dup(icomm, &comm);
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  mystatus = (nlocal > npes ? 1 : 0);

  /* see if you can actually do the sort */
  MPI_Allreduce(&mystatus, &status, 1, MPI_INT, MPI_PROD, comm);

  MPI_Barrier(comm);
  if (mype == 0) printf("%zu Starting.\n", (size_t)time(NULL));

  if (status == 0) {
    MPI_Comm_free(&comm);
    return 0;
  }
  

  /* get memory for the counts and displacements */
  scounts = gk_imalloc(npes+1, "scounts");
  rcounts = gk_imalloc(npes+1, "rcounts");
  sdispls = gk_imalloc(npes+1, "sdispls");
  rdispls = gk_imalloc(npes+1, "rdispls");

  /* get memory for the local splitters */
  mypicks  = gk_imalloc(npes+1, "mypicks");

  /* sort the local elements in increasing order */
  gk_isorti(nlocal, elmnts);

//  MPI_Barrier(comm);
//  if (mype == 0) printf("%zu Done with local sort.\n", (size_t)time(NULL));

  /* select the local npes-1 equally spaced elements */
  for (i=1; i<npes; i++) 
    mypicks[i-1] = elmnts[i*(nlocal/npes)];

  /* get memory for all splitters */
  if (mype == 0) 
    allpicks = gk_imalloc(npes*npes, "allppicks");

  /* gather the picks to the root processor */
  MPI_Gather(mypicks, npes-1, MPI_INT, allpicks, npes-1, MPI_INT, 0, comm);

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
  MPI_Bcast(mypicks, npes+1, MPI_INT, 0, comm);
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
  MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm);
  if (mype == 0) printf("%zu Done with rcounts alltoall.\n", (size_t)time(NULL));

  sdispls[0] = rdispls[0] = 0;
  for (i=0; i<npes; i++) {
    sdispls[i+1] = sdispls[i] + scounts[i];
    rdispls[i+1] = rdispls[i] + rcounts[i];
  }

  /* allocate memory for sorted elements and receive them */
  nrecv   = rdispls[npes];  
  relmnts = gk_imalloc(nrecv, "relmnts");

  MPI_Alltoallv(elmnts, scounts, sdispls, MPI_INT,
                 relmnts, rcounts, rdispls, MPI_INT, 
                 comm);
  if (mype == 0) printf("%zu Done with data alltoall.\n", (size_t)time(NULL));

  /* do the local sort of the received elements */
  gk_isorti(nrecv, relmnts);

  //MPI_Barrier(comm);
  //if (mype == 0) printf("%zu Done with 2nd sort.\n", (size_t)time(NULL));

  /* hook the fields to the data that needs to be returned */
  gk_free((void **)r_elmnts, LTERM);
  *r_nlocal = nrecv;
  *r_elmnts = relmnts;

  nrecv_min = nrecv_max = nrecv;
  MPI_Reduce(&nrecv, &nrecv_min, 1, MPI_INT, MPI_MIN, 0, comm);
  MPI_Reduce(&nrecv, &nrecv_max, 1, MPI_INT, MPI_MAX, 0, comm);

  if (mype == 0)
    printf("%zu min/max nrecv: %d %d\n", (size_t)time(NULL), nrecv_min, nrecv_max);

  //MPI_Barrier(comm);
  //if (mype == 0) printf("%zu Done.\n", (size_t)time(NULL));

  gk_free((void **)&mypicks, &rcounts, &scounts, &rdispls, &sdispls, LTERM);

  MPI_Comm_free(&comm);

  return 1;
}


/**************************************************************************/
/**************************************************************************/
int main(int argc, char **argv)
{
  int k, npes, mype, lastelmnt = 0;
  int i, nlocal, nremote;
  int *elmnts;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);

  if (argc < 2) {
    if (mype == 0)
      fprintf(stderr, "Usage: %s nelems\n", argv[0]);

    MPI_Finalize();
    return EXIT_FAILURE;
  }

  nlocal = strtol(argv[1], NULL, 10);

  if (mype == 0) 
    printf("nlocal: %d\n", nlocal);
  MPI_Barrier(MPI_COMM_WORLD);

  srand(5+mype);

  elmnts = gk_imalloc(nlocal, "elmnts");
  for (i=0; i<nlocal; i++) 
    elmnts[i] = rand();

  /* perform the sorting */
  if (samplesort(&nlocal, &elmnts, MPI_COMM_WORLD) == 0) {
    if (mype == 0)
      printf("Samplesort returned with an error: npes: %d\n", npes);
    gk_free((void **)&elmnts, LTERM);
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  //MPI_Barrier(MPI_COMM_WORLD);
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
      MPI_Send(&nlocal, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
      MPI_Send(elmnts, nlocal, MPI_INT, 0, 2, MPI_COMM_WORLD);
      gk_free((void **)&elmnts, LTERM);
    }
    else if (mype == 0) {
      MPI_Recv(&nremote, 1, MPI_INT, k, 1, MPI_COMM_WORLD, &status);
      elmnts = gk_imalloc(nremote, "elmnts");
      MPI_Recv(elmnts, nremote, MPI_INT, k, 2, MPI_COMM_WORLD, &status);
      if (elmnts[0] < lastelmnt)
        printf("Inversion: %d %d [k: %d]\n", lastelmnt, elmnts[0], k);
      //for (i=0; i<nremote; i++)
      //  printf("%5d %d\n", k, elmnts[i]);
      lastelmnt = elmnts[nremote-1];
      gk_free((void **)&elmnts, LTERM);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Finalize();

  if (mype == 0) printf("%zu Done done!\n", (size_t)time(NULL));

  return EXIT_SUCCESS;
}



