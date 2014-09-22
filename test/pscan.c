/*!
\file
\brief A program to test prefix scan operations
\author George
*/

#include <GKlib.h>
#include <bdmpi.h>

int main(int argc, char **argv)
{
  int npes, mype;
  int i, j, k, nelems;
  int *a, *b, *c;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);

  if(argc < 2) {
    if (mype == 0)
      fprintf(stderr, "usage: %s nelems\n", argv[0]);

    MPI_Finalize();
    return EXIT_FAILURE;
  }

  nelems = strtol(argv[1], NULL, 10);

  a = (int*)gk_malloc(nelems*sizeof(int), "a");
  b = (int*)gk_malloc(nelems*sizeof(int), "b");
  c = (int*)gk_malloc(nelems*sizeof(int), "b");

  for (i=0; i<nelems; i++) 
    a[i] = i*i+i+mype;

  MPI_Barrier(MPI_COMM_WORLD);

  /* Inclusive scan */
  MPI_Scan(a, b, nelems, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  gk_iset(nelems, 0, c);
  for (i=0; i<nelems; i++) {
    for (j=0; j<=mype; j++)
      c[i] += i*i+i+j;
  }

  for (i=0; i<nelems; i++) {
    if (1 || b[i] != c[i])
      printf("INCL: [%3d] %4d %4d == %4d\n", mype, a[i], b[i], c[i]);
  }

  /* Exclusive scan */
  MPI_Exscan(a, b, nelems, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  gk_iset(nelems, 0, c);
  for (i=0; i<nelems; i++) {
    for (j=0; j<mype; j++)
      c[i] += i*i+i+j;
  }

  for (i=0; i<nelems; i++) {
    if (1 || b[i] != c[i])
      printf("EXCL: [%3d] %4d %4d == %4d\n", mype, a[i], b[i], c[i]);
  }

  MPI_Finalize();

  gk_free((void**)&a, &b, &c, LTERM);

  return EXIT_SUCCESS;
}

