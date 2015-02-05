#include <GKlib.h>
#include <bdmpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int lV;         /* number of local vertices */
int lE;         /* number of local edges */
int * ia=NULL;  /* edge index array */
int * ja=NULL;  /* edge array */

int * weights;  /* global vertex weights */
int * colors;   /* global vertex colors */

int root=0; /* root process */
int rank;   /* rank in MPI communicator */
int npes;   /* number of processes */
int * off;  /* offsets for each process */
int * cnt;  /* counts for each process */

int V;  /* number of global vertices */
int E;  /* number of global edges */

// chromaticity= minimum number of colors required to color the graph
int chromaticity_upper = -1;  // upper bound on the chromaticity of the graph


/* todo: this code should be fixed so that it does not require the whole V
 * elements of colors or weights to be on each process, rather, the processes
 * should communicate the necessary interface values ... think sparse mat-vec.
 * */


int compare (const void * a, const void * b)
{
  return (*(int*)a-*(int*)b);
}


void distribute(int to)
{
  int i, num_v_per_p, remainder_v_per_p, first_v, last_v;
  MPI_Status status;
  int meta[3];

  if (rank == root) {
    meta[0] = V;
    meta[1] = E;
    meta[2] = chromaticity_upper;
    MPI_Send(meta, sizeof(meta), MPI_INT, to, 0, MPI_COMM_WORLD);
    MPI_Send(ia, cnt[to]+1, MPI_INT, to, 1, MPI_COMM_WORLD);
    MPI_Send(ja, ia[cnt[to]]-ia[0], MPI_INT, to, 2, MPI_COMM_WORLD);
  }
  else {
    MPI_Recv(meta, sizeof(meta), MPI_INT, root, 0, MPI_COMM_WORLD, &status);

    V  = meta[0];
    E  = meta[1];
    chromaticity_upper = meta[2];

    if (npes > 1) {
      num_v_per_p = V/npes;
      for (i=0; i<npes; ++i) {
        remainder_v_per_p = (V+i)%npes;
        first_v = num_v_per_p*i+remainder_v_per_p*(remainder_v_per_p<i);
        last_v = (i+1)*num_v_per_p+(remainder_v_per_p+1)*(remainder_v_per_p<i)-1;
        off[i] = first_v;
        cnt[i] = last_v-first_v+1;
      }
    }
    lV = cnt[rank];

    if (NULL == (ia=(int*)malloc((lV+1)*sizeof(int))))
      abort();

    MPI_Recv(ia, lV+1, MPI_INT, root, 1, MPI_COMM_WORLD, &status);

    lE = ia[lV]-ia[0];
    if (NULL == (ja=(int*)malloc(lE*sizeof(int))))
      abort();

    MPI_Recv(ja, lE, MPI_INT, root, 2, MPI_COMM_WORLD, &status);
  }
}


void read_graph(char * filename)
{
  int i, j, jj=0, to=0, n=0;
  int num_v_per_p, remainder_v_per_p, first_v, last_v;
  int prev_row_idx=0, row_idx=0, col_idx, off_idx=0;
  char * token;
  FILE * file;
  char line[1000];

  if (NULL == (file=fopen(filename, "r"))) {
    printf("Unable to open file %s\n", filename);
    return;
  }

  //read file line by line
  i = 0;
  j = 0;
  while (NULL != fgets(line, 1000, file)) {
    //tokenize line
    strtok(line, " ");

    //read number of vertices and edges
    if (0 == strcmp(line, "p")) {
      token = (char *)strtok(NULL," ");
      token = (char *)strtok(NULL," ");
      V = atoi(token);
      token = (char *)strtok(NULL," ");
      E = atoi(token);

      /*if (V%2 == 0)
        chromaticity_upper = V-1;
      else
        chromaticity_upper = V;*/
      chromaticity_upper = 256;

      if (npes > 1) {
        num_v_per_p = V/npes;
        for (i=0; i<npes; ++i) {
          remainder_v_per_p = (V+i)%npes;
          first_v = num_v_per_p*i+remainder_v_per_p*(remainder_v_per_p<i);
          last_v = (i+1)*num_v_per_p+(remainder_v_per_p+1)*(remainder_v_per_p<i)-1;
          off[i] = first_v;
          cnt[i] = last_v-first_v+1;
        }
      }
    }

    //read edges into graph matrix
    //1.Allocate memory for graph
    if (NULL == ia && V > 0) {
      if (NULL == (ia=(int *) malloc((V+1)*sizeof(int))))
        abort();
      if (NULL == (ja=(int *) malloc(E*sizeof(int))))
        abort();
      ia[0] = 0;
    }

    //2.then load edges into it
    if (0 == strcmp(line, "e")) {
      token = (char *)strtok(NULL," ");
      row_idx = atoi(token)-1;//0 based index
      token = (char *)strtok(NULL," ");
      col_idx = atoi(token)-1;//0 based index

      if (row_idx != prev_row_idx) {
        ia[row_idx-off_idx] = jj;

        if (to < npes-1 && row_idx == off[to+1]) {
          ia[row_idx+1-off_idx] = jj;

          distribute(to++);

          off_idx = row_idx;
          jj = 0;
        }

        prev_row_idx = row_idx;
      }
      ja[jj++] = col_idx;
      j++;
    }
  }
  ia[row_idx+1-off_idx] = jj;

  lV = cnt[rank];
  lE = ia[lV]-ia[0];

  if (NULL == (ia=(int*)realloc(ia, (lV+1)*sizeof(int))))
    abort();
  if (NULL == (ja=(int*)realloc(ja, lE*sizeof(int))))
    abort();

  if (row_idx+1 != V)
    printf("Incorrect vertex reading: There are %d vertices but read %d.\n",
      V, row_idx);
  if (j != E)
    printf("Incorrect edge reading: There are %d edges but read %d.\n", E, j);

  fclose(file);
}


void jones_plassmann(void)
{
  int i, j, k, jj;
  int i_weight, num_colors, min_color;
  int i_weight_is_max;
  int * i_colors,* neighbor_colors;

  // Go through the vertices of this process and compare their weights with
  // neighboring vertices to find which of them are local maxima. Those form
  // the independent set in each iteration
  // The minimum number of colors, i.e. the chromaticity of the graph
  // can not be larger than chromaticity_upper so only iterate that many times

  if (NULL == (i_colors=(int *) calloc(lV, sizeof(int))))
    abort();

  for (k=0; k<chromaticity_upper; ++k) {
    for (i=0; i<lV; ++i) {
      //get the vertex weight
      i_weight        = weights[off[rank]+i];
      i_weight_is_max = 1;
      num_colors      = 0;

      if (NULL == (neighbor_colors=(int *) calloc(V, sizeof(int))))
        abort();

      // compare vertex weight to weights of its non-colored neighbors to see
      // if it is a maximum. Also gather the colors of all neighbors of the
      // vertex i that have been colored
      for (j=ia[i]; j<ia[i+1]; ++j) {
        jj = ja[j];
        if (0 != colors[jj]) {
          //if neighbor is colored just add its color to the neighbor_colors
          neighbor_colors[num_colors++] = colors[jj];
        }
        else if (i_weight < weights[jj] ||
                 (i_weight == weights[jj] && jj > i)) {
          //if the weights match, solve conflict by looking at the vertices
          //ids and taking the vertex with higher id as the max
          i_weight_is_max = 0;
          break;
        }
      }

      // if the vertex weight is a max and vertex hasnt been colored, color it
      // with the smallest color possible that is not one of neighbor_colors
      if (1 == i_weight_is_max && 0 == colors[off[rank]+i]) {
        /* find smallest color to assign to the j vertex that color is either
            a)  1 if none of the neighbors is colored or the smallest color of
                a neighbor is >1
            b)  In between a color in the array of neighbors colors if there
                is a gap between two of the (sorted) neighbors colors
            c)  1 more than the last color in the sorted array of neighbors
                colors sort neighbors colors. */
        qsort(neighbor_colors, num_colors, sizeof(int), compare);

        if (0 == num_colors || 1 < neighbor_colors[0]) {
          min_color = 1;
        }
        else {
          for (j=0; j<num_colors-1; ++j) {
            if (1 < neighbor_colors[j+1]-neighbor_colors[j]) {
              min_color = neighbor_colors[j]+1;
              break;
            }
          }
          if (j == num_colors-1)
            min_color = neighbor_colors[num_colors-1]+1;
        }

        i_colors[i] = min_color;
      }

      free(neighbor_colors);
    }

    MPI_Allgatherv(i_colors, lV, MPI_INT, colors, cnt, off, MPI_INT,
      MPI_COMM_WORLD);
  }
  free(i_colors);
}


int main(int argc, char * argv[])
{
  int i, j, k, jj;
  char * input_filename;

  //Initialize
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);

  root = npes-1;

  if (argc != 2) {
    if (rank == root)
      printf("Usage: graphcoloring input_filename\n");

    MPI_Finalize();
    return -1;
  }

  if (NULL == (input_filename=malloc(900)))
    abort();
  strncpy(input_filename, argv[1], 900);

  /* allocate offsets */
  if (NULL == (off=(int *)malloc(npes*sizeof(int))))
    abort();
  if (NULL == (cnt=(int *)malloc(npes*sizeof(int))))
    abort();

  /* only root reads the file and loads the full graph */
  if (rank == root) {
    read_graph(input_filename);
    printf("V=%d E=%d Chromaticity Upper Bound=%d\n", V, E, chromaticity_upper);
    fflush(stdout);
  }
  else {
    distribute(rank);
  }

#if 0
  for (i=0; i<npes; ++i) {
    if (rank == i) {
      for (j=0; j<lV; ++j) {
        for (k=ia[j]; k<ia[j+1]; ++k)
          printf("e %d %d\n", off[rank]+j+1, ja[k]+1);
      }
    }
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif

  printf("p[%d] gets %d,%d vertices and %d edges starting at vertex "
    "%d\n", rank, lV, cnt[rank], lE, off[rank]);

  /* randomly initialize weights */
  if (NULL == (weights=(int *)malloc(V*sizeof(int))))
    abort();
  srand(5);
  for (i=0; i<V; ++i)
    weights[i] = rand()%(V*1000);

  /* initialize colors to 0 */
  if (NULL == (colors=(int *)calloc(V, sizeof(int))))
    abort();

  /* Jones-Plassman algorithm */
  jones_plassmann();

  /* find out how many vertices are uncolered */
  if (rank == root) {
    int num_uncolored=0;
    int * map = (int *) calloc(V, sizeof(int));
    if (NULL == map)
      abort();
    for (i=0; i<V; ++i) {
      if (0 == colors[i])
        num_uncolored++;

      map[colors[i]] = 1;
    }
    if (0 != num_uncolored)
      printf("Not all vertices have been colored\n");
    for (i=1; i<V; ++i) {
      if (0 == map[i])
        break;
    }
    printf("num colors=%d\n", i);
    for (; i<V; ++i) {
      if (0 != map[i])
        abort();
    }
    free(map);
  }
  /*for (i=0; i<lV; ++i) {
    for (j=ia[i]; j<ia[i+1]; ++j) {
      if (colors[off[rank]+i] == colors[ja[j]])
        abort();
    }
  }*/

  free(input_filename);
  free(ia);
  free(ja);
  free(weights);
  free(colors);
  free(off);

  MPI_Finalize();

  return 0;
}
