#include "spmv.h"

double gk_clockseconds()
{
#ifdef __GNUC__
  struct timeval ctime;
  gettimeofday(&ctime, NULL);
  return (double)ctime.tv_sec + (double).000001*ctime.tv_usec;
#else
  return (double)time(NULL);
#endif
}

double ji_memusage(int sz)
{
  struct rusage r;
			
	getrusage(RUSAGE_SELF, &r);
	return (r.ru_maxrss*1.0/sz);
}

int min(const int a, const int b)
{
  if (a < b)
    return a;
  return b;
}

int cmp0 (const void * a, const void * b)
{
  return (*(int *)a - *(int *)b);
}

int cmp1 (const void * a, const void * b)
{
  return ((*(kv_t *)a).k - (*(kv_t *)b).k);
}

void read_files(const char * f1, const char * f2, csr_t * mat, vec_t * vec)
{
  int     m, nnz, lnnz;
  int     c, r, r_;
  double  v;
  FILE    * fp;
  int     * rnnz, * rowind, * colptr;
  double  * values, * bvec;

  m   = 1;
  nnz = 0;
  r_  = 0;

  /***** READ A MATRIX *****/
  fp = fopen(f1, "r");

  while (fscanf(fp, "%d %d %lf\n", &r, &c, &v) != EOF) {
    nnz++;
    if (r != r_) {
      m++;
      r_ = r;
    }
  }

  rnnz   = mat->rnnz    = (int *) calloc((m), sizeof(int));
  rowind = mat->rowind  = (int *) malloc((m+1)*sizeof(int));
  colptr = mat->colptr  = (int *) malloc((nnz)*sizeof(int));
  values = mat->values  = (double *) malloc((nnz)*sizeof(double));

  rewind(fp);
  m     = 1;
  nnz   = 0;
  lnnz  = nnz;
  r_    = 0;

  rowind[0] = nnz;
  while (fscanf(fp, "%d %d %lf\n", &r, &c, &v) != EOF) {
    colptr[nnz]   = c;
    values[nnz++] = v;
    if (r != r_ && nnz != 1) {
        rnnz[m-1] = nnz-lnnz;
        lnnz      = nnz;

      rowind[m++] = nnz;
      r_ = r;
    }
  }
  rowind[m] = nnz;
  rnnz[m-1] = nnz-lnnz;

  fclose(fp);

  mat->m    = m;
  mat->n    = m;
  mat->nnz  = nnz;
  /*************************/

  /***** READ B VECTOR *****/
  fp = fopen(f2, "r");

  m = 0;

  while (fscanf(fp, "%lf\n", &v) != EOF)
    m++;

  bvec = vec->values = (double *) malloc((m)*sizeof(double));

  rewind(fp);
  m = 0;

  while (fscanf(fp, "%lf\n", &v) != EOF)
    bvec[m++] = v;

  fclose(fp);

  vec->m = m;
  /*************************/
}

void write_file(const char * f1, vec_t * vec)
{
  int     i, m;
  FILE    * fp;
  double  * values;
  char    fn[MAXLINE];

  m       = vec->m;
  values  = vec->values;

  sprintf(fn, "%s.sol", f1);

  fp = fopen(fn, "w");

  for (i=0; i<m; ++i)
    fprintf(fp, "%f\n", values[i]);

  fclose(fp);
}

void scatter_data(csr_t * imat, csr_t * mat, vec_t * ivec, vec_t * vec, int np, 
      int rank)
{
  int i, j, k, m, mm, nnz, tnnz;
  int * sendcounts, * displs;

  m   = imat->m;
  mm  = mat->m;

  if (rank == 0) {
    sendcounts  = (int *) malloc(np*sizeof(int));
    displs      = (int *) malloc(np*sizeof(int));

    for (i=0; i<np-1; ++i) {
      sendcounts[i] = mm;
      displs[i]     = i*mm;
    }
    sendcounts[np-1] = m-(np-1)*mm;
    displs[np-1]     = (np-1)*mm;
  }

  mat->rnnz    = (int *) malloc((mat->m)*sizeof(int));
  mat->rowind  = (int *) malloc((mat->m+1)*sizeof(int));
  vec->local   = (double *) malloc((vec->m)*sizeof(double));

  MPI_Scatterv(imat->rnnz, sendcounts, displs, MPI_INT, mat->rnnz, mm, MPI_INT, 
    0, MPI_COMM_WORLD);
  MPI_Scatterv(ivec->values, sendcounts, displs, MPI_DOUBLE, vec->local, mm, 
    MPI_DOUBLE, 0, MPI_COMM_WORLD);

  mat->rowind[0] = 0;
  for (i=1, mat->nnz=0; i<mat->m+1; ++i) {
    mat->nnz += mat->rnnz[i-1];
    mat->rowind[i] = mat->nnz;
  }

  mat->colptr  = (int *) malloc((mat->nnz)*sizeof(int));
  mat->values  = (double *) malloc((mat->nnz)*sizeof(double));

  if (rank == 0) {
    for (i=0, tnnz=0; i<np; ++i) {
      k = min((i+1)*mm, m);
      for (j=i*mm, nnz=0; j<k; ++j)
        nnz += imat->rnnz[j];
      sendcounts[i] = nnz;
      displs[i]     = tnnz;
      tnnz          += nnz;
    }
  }

  MPI_Scatterv(imat->colptr, sendcounts, displs, MPI_INT, mat->colptr, mat->nnz,
    MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatterv(imat->values, sendcounts, displs, MPI_DOUBLE, mat->values, 
    mat->nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#if DEBUG == 1
  printf("rank %d\n", rank);
  for (i=0; i<mat->nnz; ++i)
    printf(" %d", mat->colptr[i]);
  printf("\n");
  for (i=0; i<vec->m; ++i)
    printf(" %.2f", vec->local[i]);
  printf("\n");
#endif

  if (rank == 0) {
    free(sendcounts);
    free(displs);
  }
}

void setup_rinfo(cinfo_t * rinfo, int m, int nnz, int * ineed, int np, int rank)
{
  int i, j, l, mm, nunique;

  mm = ceil(m*1.0/np);

  rinfo->ind      = (int *) malloc(nnz*sizeof(int));
  rinfo->nbr_ids  = (int *) malloc(np*sizeof(int));
  rinfo->ptr      = (int *) malloc((np+1)*sizeof(int));

  for (i=0, j=0; i<nnz; ++i) {
    if (ineed[i] < rank*mm || ineed[i] >= (rank+1)*mm)
      rinfo->ind[j++] = ineed[i];
  }

  qsort(rinfo->ind, j, sizeof(int), cmp0);

  for (i=1, nunique=1; i<j; ++i) {
    if (rinfo->ind[i] != rinfo->ind[i-1])
      rinfo->ind[nunique++] = rinfo->ind[i];
  }

  rinfo->n_nbrs = 0;
  rinfo->ptr[0] = 0;
  for (i=0, j=0; i<np; ++i) {
    l = j;
    for (; j<nunique; ++j) {
      if (rinfo->ind[j] >= (i+1)*mm)
        break;
    }
    if (j>l) {
      rinfo->nbr_ids[rinfo->n_nbrs++] = i;
      rinfo->ptr[rinfo->n_nbrs]       = j;
    }
  }
}

void setup_sinfo(cinfo_t * sinfo, cinfo_t * rinfo, int np, int rank)
{ 
  int         i, nsend;
  int         * recv_row, * send_row;
  MPI_Request * send_req, * recv_req;
  MPI_Status  * status;

  recv_req  = (MPI_Request *) malloc(np*sizeof(MPI_Request));
  send_req  = (MPI_Request *) malloc(np*sizeof(MPI_Request));
  status    = (MPI_Status *) malloc(np*sizeof(MPI_Status));
  recv_row  = (int *) malloc(np*sizeof(int));
  send_row  = (int *) malloc(np*sizeof(int));

  for (i=0; i<np; ++i)
    recv_row[i] = 0;
  for (i=0; i<rinfo->n_nbrs; ++i)
    recv_row[rinfo->nbr_ids[i]] = rinfo->ptr[i+1]-rinfo->ptr[i];

  MPI_Alltoall(recv_row, 1, MPI_INT, send_row, 1, MPI_INT, MPI_COMM_WORLD);

  sinfo->n_nbrs = 0;
  for (i=0, nsend=0; i<np; ++i) {
    if (send_row[i] > 0) {
      sinfo->n_nbrs++;
      nsend += send_row[i];
    }
  }

  sinfo->nbr_ids  = (int *) malloc(sinfo->n_nbrs*sizeof(int));
  sinfo->ptr      = (int *) malloc((sinfo->n_nbrs+1)*sizeof(int));
  sinfo->ind      = (int *) malloc(nsend*sizeof(int));

  sinfo->n_nbrs = 0;
  sinfo->ptr[0] = 0;
  for (i=0, nsend=0; i<np; ++i) {
    if (send_row[i] > 0) {
      nsend                           += send_row[i];
      sinfo->nbr_ids[sinfo->n_nbrs++] = i;
      sinfo->ptr[sinfo->n_nbrs]       = nsend;
    }
  }

  for (i=0; i<sinfo->n_nbrs; ++i)
    MPI_Irecv(&sinfo->ind[sinfo->ptr[i]], sinfo->ptr[i+1]-sinfo->ptr[i], 
      MPI_INT, sinfo->nbr_ids[i], 1, MPI_COMM_WORLD, &recv_req[i]);

  for (i=0; i<rinfo->n_nbrs; ++i)
    MPI_Isend(&rinfo->ind[rinfo->ptr[i]], rinfo->ptr[i+1]-rinfo->ptr[i], 
      MPI_INT, rinfo->nbr_ids[i], 1, MPI_COMM_WORLD, &send_req[i]);

  MPI_Waitall(sinfo->n_nbrs, recv_req, status);
  MPI_Waitall(rinfo->n_nbrs, send_req, status);

  free(recv_row);
  free(send_row);
  free(recv_req);
  free(send_req);
  free(status);

  if (rank == -1) {} /* surpress unused parameter warning */
}

void exchange_nonlocal(int m, double * local_val, double ** recv_val_, 
      cinfo_t * sinfo, cinfo_t * rinfo, int np, int rank)
{
  int         i, firstval;
  double      * send_val, * recv_val;
  MPI_Request * send_req, * recv_req;
  MPI_Status  * status;

  recv_req  = (MPI_Request *) malloc(np*sizeof(MPI_Request));
  send_req  = (MPI_Request *) malloc(np*sizeof(MPI_Request));
  status    = (MPI_Status *) malloc(np*sizeof(MPI_Status));
  recv_val  = (*recv_val_)  = \
    (double *) malloc(rinfo->ptr[rinfo->n_nbrs]*sizeof(double));
  send_val  = (double *) malloc(sinfo->ptr[sinfo->n_nbrs]*sizeof(double));

  for (i=0; i<rinfo->n_nbrs; ++i)
    MPI_Irecv(&recv_val[rinfo->ptr[i]], rinfo->ptr[i+1]-rinfo->ptr[i], 
      MPI_DOUBLE, rinfo->nbr_ids[i], 1, MPI_COMM_WORLD, &recv_req[i]);

  firstval = ceil(m*1.0/np)*rank;
  for (i=0; i<sinfo->ptr[sinfo->n_nbrs]; ++i)
    send_val[i] = local_val[sinfo->ind[i]-firstval];

  for (i=0; i<sinfo->n_nbrs; ++i)
    MPI_Isend(&send_val[sinfo->ptr[i]], sinfo->ptr[i+1]-sinfo->ptr[i], 
      MPI_DOUBLE, sinfo->nbr_ids[i], 1, MPI_COMM_WORLD, &send_req[i]);

  MPI_Waitall(rinfo->n_nbrs, recv_req, status);
  MPI_Waitall(sinfo->n_nbrs, send_req, status);

#if DEBUG == 1
  for (i=0; i<rinfo->ptr[rinfo->n_nbrs]; ++i)
    printf(" %.2f", recv_val[i]);
  printf("\n");
#endif

  free(send_val);
  free(recv_req);
  free(send_req);
  free(status);
}

void cat(int mm, int nn, double * local, double * nonlocal, double ** values_,
      int np, int rank)
{
#if DEBUG == 1
  int i;
#endif
  double * values;

  values = (*values_) = (double *) malloc((mm+nn)*sizeof(double));

  memcpy(values, local, mm*sizeof(double));
  memcpy(&values[mm], nonlocal, nn*sizeof(double));

#if DEBUG == 1
  for (i=0; i<mm+nn; ++i)
    printf(" %.2f", values[i]);
  printf("\n");
#endif

  free(local);
  free(nonlocal);

  if (np == 0 || rank == 0) {} /* to surpress unused argument warning */
}

void g2l(int m, int nnz, int * ineed, cinfo_t * rinfo, int ** l2g_, int np, 
      int rank)
{
  int   i, j, k, mm, mmm;
  int   firstval;
  int   * l2g;
  kv_t  * kv;

  mm  = ceil(m*1.0/np);
  mmm = mm;

  if ((rank+1)*mm > m)
    mmm = m-(rank*mm);

  kv  = (kv_t *) malloc(nnz*sizeof(kv_t));
  l2g = (*l2g_) = (int *) malloc((mmm+rinfo->ptr[rinfo->n_nbrs])*sizeof(int));

  firstval = rank*mm;
  for (i=0, j=0; i<nnz; ++i) {
    if (ineed[i] < firstval || ineed[i] >= firstval+mm) {
      k++;
      kv[j].k = ineed[i];
      kv[j++].v = i;
    }else {
      ineed[i] -= firstval;
    }
  }

  qsort(kv, j, sizeof(kv_t), cmp1);

  kv[j].k = -1;
  for (i=0, k=mmm; i<j; ++i) {
    ineed[kv[i].v] = k;
    if (kv[i].k != kv[i+1].k)
      l2g[k++] = kv[i].k;
  }

  for (i=0; i<mmm; ++i)
    l2g[i] = firstval+i;

#if DEBUG == 1
  for (i=0; i<mmm+rinfo->ptr[rinfo->n_nbrs]; ++i)
    printf(" %d", l2g[i]);
  printf("\n");
#endif

  free(kv);
}

void l2g(int n, int * l2g, int * ineed, int np, int rank)
{
  int i;

  for (i=0; i<n; ++i)
    ineed[i] = l2g[ineed[i]];

#if DEBUG == 1
  for (i=0; i<n; ++i)
    printf(" %d", ineed[i]);
  printf("\n");
  printf("-----\n");
#endif

  if (np == 0 || rank == 0) {} /* to surpress unused argument warning */
}

void sparsemv(csr_t * mat, vec_t * bvec, vec_t * yvec, int np, int rank)
{
  int     i, j, m, nnz;
  int     * rnnz, * colptr;
  double  * y, * mat_val, * vec_val;

  m       = mat->m;
  rnnz    = mat->rnnz;
  colptr  = mat->colptr;
  mat_val = mat->values;
  vec_val = bvec->values;
  y       = yvec->values = (double *) malloc(m*sizeof(double));

  for (i=0, nnz=0; i<m; ++i) {
    for (j=0, y[i]=0.0; j<rnnz[i]; ++j) {
      y[i] += (mat_val[nnz] * vec_val[colptr[nnz]]);
      nnz++;
    }
  }
  yvec->m = m;

  if (np == 0 || rank == 0) {} /* to surpress unused argument warning */
}

void gather_data(int m, vec_t * yvec, vec_t * ovec, int np, int rank)
{
  int i, mm;
  int * recvcounts, * displs;

  mm          = yvec->m;
  recvcounts  = NULL;
  displs      = NULL;

  if (rank == 0) {
    ovec->m       = m;
    ovec->values  = (double *) malloc(m*sizeof(double));
    recvcounts    = (int *) malloc(np*sizeof(int));
    displs        = (int *) malloc(np*sizeof(int));

    for (i=0; i<np-1; ++i) {
      recvcounts[i] = mm;
      displs[i]     = i*mm;
    }
    recvcounts[np-1] = m-(np-1)*mm;
    displs[np-1]     = (np-1)*mm;
  }

  MPI_Gatherv(yvec->values, mm, MPI_DOUBLE, ovec->values, recvcounts, displs, 
    MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    free(recvcounts);
    free(displs);
  }
}

