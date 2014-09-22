#include "spmv.h"

int main(int argc, char** argv)
{
  int     m, n, nnz, vm;
  int     rank, np;
  double  tmr0, tmr1, tmr1m, tmr234, tmr5, tmr5p, tmr6, tmr7, tmr7m;
  csr_t   imat, mat;
  vec_t   ivec, bvec, yvec, ovec;
  cinfo_t rinfo, sinfo;

  /************
  * setup mpi *
  ************/
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  /***********/

  /***************************
  * initialize the variables *
  ***************************/
  m       = 0;
  n       = 0;
  nnz     = 0;
  vm      = 0;
  tmr0    = 0.0;
  tmr1    = 0.0;
  tmr1m   = 0.0;
  tmr234  = 0.0;
  tmr5    = 0.0;
  tmr5p   = 0.0;
  tmr6    = 0.0;
  tmr7    = 0.0;
  tmr7m   = 0.0;
  /**************************/

  /*******************
  * read input files *
  *******************/
  if (rank == 0) {
    start_timer(tmr0);
    start_timer(tmr1);
    start_timer(tmr1m);

    read_files(argv[1], argv[2], &imat, &ivec);

    stop_timer(tmr1m);

    m   = imat.m;
    n   = imat.n;
    nnz = imat.nnz;
    vm  = ivec.m;
  }
  /******************/

  /***********************************
  * broadcast the necessary counters *
  ***********************************/
  MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&vm, 1, MPI_INT, 0, MPI_COMM_WORLD);
  /**********************************/

  /*****************************************
  * adjust local values for number of rows *
  *****************************************/
  mat.m   = ceil(m*1.0/np);
  if ((rank+1)*mat.m > m)
    mat.m = m - (rank*mat.m);

  bvec.m  = ceil(vm*1.0/np);
  if ((rank+1)*bvec.m > vm)
    bvec.m = vm - (rank*bvec.m);
  /****************************************/

  /***************************
  * distibute the input data *
  ***************************/
  scatter_data(&imat, &mat, &ivec, &bvec, np, rank);

  if (rank == 0)
    stop_timer(tmr1);
  /**************************/

  /********************
  * free input memory *
  ********************/
  if (rank == 0) {
    free(imat.rnnz);
    free(imat.rowind);
    free(imat.colptr);
    free(imat.values);
    free(ivec.values);
  }
  /*******************/

  /**********************
  * setup exchange info *
  ***********************/
  if (rank == 0)
    start_timer(tmr234);

  setup_rinfo(&rinfo, m, mat.nnz, mat.colptr, np, rank);
  setup_sinfo(&sinfo, &rinfo, np, rank);

  if (rank == 0)
    stop_timer(tmr234);
  /**********************/

  /***************************
  * exchange nonlocal values *
  ***************************/
  if (rank == 0)
    start_timer(tmr5);

  exchange_nonlocal(m, bvec.local, &bvec.nonlocal, &sinfo, &rinfo, np, 
    rank);

  if (rank == 0)
    stop_timer(tmr5);
  /**************************/

  /*************************************************
  * reindex m.colptr to match the concatenation of *
  * local and nonlocal b elements                  *
  *************************************************/
  if (rank == 0)
    start_timer(tmr5p);

  g2l(m, mat.nnz, mat.colptr, &rinfo, &bvec.l2g, np, rank);
  /************************************************/

  /********************************************
  * concatenate local and nonlocal b elements *
  ********************************************/
  cat(bvec.m, rinfo.ptr[rinfo.n_nbrs], bvec.local, bvec.nonlocal, 
    &bvec.values, np, rank);
  /*******************************************/

  /*************************************************************
  * make sure that all nodes are completed to this point       *
  * since routines have run without any forced synchronization *
  * due to communication                                       *
  *************************************************************/
  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0)
    stop_timer(tmr5p);
  /************************************************************/

  /*******************************************************
  * perform the local sparse matrix-vector multiplcation *
  *******************************************************/
  if (rank == 0)
    start_timer(tmr6);

  sparsemv(&mat, &bvec, &yvec, np, rank);
  /******************************************************/

  /*************************************************************
  * make sure that all nodes are completed to this point       *
  * since routines have run without any forced synchronization *
  * due to communication                                       *
  *************************************************************/
  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0)
    stop_timer(tmr6);
  /************************************************************/

  /******************************
  * collect the data for output *
  ******************************/
  if (rank == 0)
    start_timer(tmr7);

  gather_data(m, &yvec, &ovec, np, rank);
  /*****************************/

  /*************
  * write file *
  *************/
  if (rank == 0) {
    start_timer(tmr7m);

    write_file(argv[2], &ovec);

    stop_timer(tmr7m);
    stop_timer(tmr7);
  }
  /************/

 /**********************
  * free program memory *
  **********************/
  free(mat.rnnz);
  free(mat.rowind);
  free(mat.colptr);
  free(mat.values);
  free(bvec.values);
  if (rank == 0)
    free(ovec.values);
  /*********************/

  /***************************************
  * output memory and timing information *
  ***************************************/
  if (rank == 0) {
    stop_timer(tmr0);

    printf("-----timing-----\n");
    printf("s1:   %.2fs\n", tmr1);
    printf("s1m:  %.2fs\n", tmr1m);
    printf("s234: %.2fs\n", tmr234);
    printf("s5:   %.2fs\n", tmr5);
    printf("s5p:  %.2fs\n", tmr5p);
    printf("s6:   %.2fs\n", tmr6);
    printf("s7:   %.2fs\n", tmr7);
    printf("s7m:  %.2fs\n", tmr7m);
    printf("s1-7: %.2fs\n", tmr0);
    printf("-----memory-----\n");
  }

  MPI_Barrier(MPI_COMM_WORLD);
  printf("m%d:   %.2fMB\n", rank, ji_memusage(1000));
  /**************************************/

  /*************
  * finish mpi *
  *************/
  MPI_Finalize();
  /************/

  return EXIT_SUCCESS;
}
