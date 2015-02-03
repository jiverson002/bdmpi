#include <GKlib.h>
#include <bdmpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef int bool;
#define true 1
#define false 0

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

int V;  /* number of global vertices */
int E;  /* number of global edges */

// chromaticity= minimum number of colors required to color the graph
int chromaticity_upper = -1;  // upper bound on the chromaticity of the graph
int *graph=NULL;  // the symmetric adjacency graph matrix. The index of
                  // element at col,row is idx = row*V+col. Similarly the
                  // element at index i is at row = i/V,col=i%V


int compare (const void * a, const void * b)
{
  return (*(int*)a-*(int*)b);
}


void distribute(int to, int n)
{
  int i;
  MPI_Status status;
  int meta[4];

  if (rank == root) {
    meta[0] = V;
    meta[1] = E;
    meta[2] = chromaticity_upper;
    meta[3] = n;
    MPI_Send(meta, sizeof(meta), MPI_INT, to, 0, MPI_COMM_WORLD);
    MPI_Send(ia, n+1, MPI_INT, to, 1, MPI_COMM_WORLD);
    MPI_Send(ja, ia[n]-ia[0], MPI_INT, to, 2, MPI_COMM_WORLD);
  }
  else {
    MPI_Recv(meta, sizeof(meta), MPI_INT, root, 0, MPI_COMM_WORLD, &status);

    V   = meta[0];
    E   = meta[1];
    chromaticity_upper = meta[2];
    lV  = meta[3];

    if (NULL == (off=(int*)malloc(npes*sizeof(int))))
      abort();

    if (npes > 1) {
      off[0] = V/npes;
      for (i=1; i<npes; ++i)
        off[i] = off[i-1]+V/npes;
    }

    if (NULL == (ia=(int*)malloc((lV+1)*sizeof(int))))
      abort();

    MPI_Recv(ia, lV+1, MPI_INT, root, 1, MPI_COMM_WORLD, &status);

    lE = ia[lV]-ia[0];
    if (NULL == (ja=(int*)malloc(lE*sizeof(int))))
      abort();

    MPI_Recv(ja, lE, MPI_INT, root, 2, MPI_COMM_WORLD, &status);

    printf("p[%d] gets %d vertices and %d edges starting at vertex "
      "%d\n", to, lV, lE, V/npes);
  }
}


void read_graph(char * filename)
{
  int i, j, jj=0, to=0;
  int prev_row_idx=0, row_idx=0, col_idx, off_idx=0;
  char * token;
  FILE * file;
  char line[1000];

  if (NULL == (file=fopen(filename,"r"))) {
    printf("Unable to open file %s\n",filename);
    return;
  }

  //read file line by line
  i = 0;
  j = 0;
  while (fgets(line,1000,file) != NULL) {
    //tokenize line
    strtok(line," ");
    //read number of vertices and edges
    if (strcmp(line,"p")==0) {
      token = (char *)strtok(NULL," ");
      token = (char *)strtok(NULL," ");
      V = atoi(token);
      token = (char *)strtok(NULL," ");
      E = atoi(token);

      lV = V/npes;
    }

    //read edges into graph matrix
    //1.Allocate memory for graph
    if (NULL == graph && V>0) {
      if (NULL == (ia=(int *) malloc((V+1)*sizeof(int))))
        abort();
      if (NULL == (ja=(int *) malloc(E*sizeof(int))))
        abort();
      ia[0] = 0;

      if (NULL == (graph=(int *) malloc(V*V*sizeof(int))))
        abort();
      memset(graph,0,V*V*sizeof(int));
    }

    //2.then load edges into it
    if (strcmp(line,"e") == 0) {
      token = (char *)strtok(NULL," ");
      row_idx = atoi(token)-1;//0 based index
      token = (char *)strtok(NULL," ");
      col_idx = atoi(token)-1;//0 based index

      if (row_idx != prev_row_idx) {
        ia[row_idx-off_idx] = jj;

        if (0 == row_idx%lV && to < npes) {
          ia[row_idx+1-off_idx] = jj;

          distribute(to++, row_idx-off_idx);

          off_idx = row_idx;
          jj = 0;
        }

        prev_row_idx = row_idx;
      }
      ja[jj++] = col_idx;

      graph[row_idx*V+col_idx]=1;
      graph[col_idx*V+row_idx]=1;//symmetric matrix, unidirectional graph
      j++;
    }
  }
  ia[row_idx+1-off_idx] = jj;

  lV = row_idx+1-off_idx;
  lE = ia[lV];

  printf("p[%d] gets %d vertices and %d edges starting at vertex %d\n", root,
    lV, lE, off_idx);
  if (NULL == (ia=(int*)realloc(ia, (lV+1)*sizeof(int))))
    abort();
  if (NULL == (ja=(int*)realloc(ja, lE*sizeof(int))))
    abort();

  if (row_idx+1 != V)
    printf("Incorrect vertex reading: There are %d vertices but read %d.\n",
      V, row_idx);
  if (j != E && j != E/2)
    printf("Incorrect edge reading: There are %d edges but read %d.\n", E, j);

  if (V%2 == 0)
    chromaticity_upper = V-1;
  else
    chromaticity_upper = V;

  fclose(file);
}


void distribute_graph(int ** range, int ** vertex_offsets, int ** p_graph)
{
  int i, first_v, last_v, num_v_per_p, remainder_v_per_p;
  int * p_graph_size, * offsets;

  //root broadcasts E,V and chromaticity_upper to all other processes
  MPI_Bcast(&V, 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(&E, 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(&chromaticity_upper, 1, MPI_INT, root, MPI_COMM_WORLD);

  //Distribute the number of vertices as uniformly as possible among the
  //processors: that means (V+rank)/npes vertices per process (integer
  //division) with npes being the processor id 0..np-1
  if (NULL == (p_graph_size = (int *)malloc(npes*sizeof(int))))
    abort();
  if (NULL == (offsets = (int *)malloc(npes*sizeof(int))))
    abort();
  if (NULL == (*vertex_offsets = (int *)malloc(npes*sizeof(int))))
    abort();
  if (NULL == (*range = (int *)malloc(npes*sizeof(int))))
    abort();
  num_v_per_p     = V/npes;

  for (i=0; i<npes; ++i) {
    remainder_v_per_p = (V+i)%npes;
    //index of the first vertex for process i
    first_v = num_v_per_p*i+remainder_v_per_p*(remainder_v_per_p<i);
    //index of the last vertex for process i
    last_v = (i+1)*num_v_per_p+(remainder_v_per_p+1)*(remainder_v_per_p<i)-1;
    (*range)[i] = last_v-first_v + 1;
    p_graph_size[i] = (*range)[i]*V;

    offsets[0] = 0;
    (*vertex_offsets)[0]=0;
    if (i>0) {
      offsets[i] = offsets[i-1] + p_graph_size[i-1];
      (*vertex_offsets)[i] = (*vertex_offsets)[i-1] + (*range)[i-1];
    }
  }

  //root sends portion of the graph array corresponding to vertices of the
  //process to each process
  if (NULL == (*p_graph=(int *) malloc(p_graph_size[rank]*sizeof(int))))
    abort();

  MPI_Scatterv(graph, p_graph_size, offsets, MPI_INT, *p_graph,
    p_graph_size[rank], MPI_INT, root, MPI_COMM_WORLD);

  free(p_graph_size);
  free(offsets);
}


void jones_plassmann(int* range, int * vertex_offsets, int * p_graph)
{
  int i,j, k;
  int num_v_per_p,remainder_v_per_p,first_v;
  int * j_colors,*neighbor_colors;  // i_colors holds the colors of this
                                    // vertices process, neighbor_colors the
                                    // colors of one of these
  int j_weight, num_colors, min_color;
  bool j_weight_is_max;
#if 0
  int * i_colors,*neighbor_colors;  // i_colors holds the colors of this
                                    // vertices process, neighbor_colors the
                                    // colors of one of these
  int i_weight, num_colors, min_color;
  bool i_weight_is_max;
#endif

  //Go through the vertices of this process and compare their weights with
  //neighboring vertices to find which of them are local maxima. Those form
  //the independent set in each iteration
  //The minimum number of colors, i.e. the chromaticity of the graph
  //can not be larger than chromaticity_upper so only iterate that many times
  num_v_per_p = V/npes;
  if (NULL == (j_colors=(int *) calloc(range[rank], sizeof(int))))
    abort();

#if 0
  for (k=0; k<chromaticity_upper; ++k) {
    for (i=0; i<lV; ++i) {
      //get the vertex weight
      i_weight        = weights[off+i];
      i_weight_is_max = 1;
      num_colors      = 0;

      if (NULL == (neighbor_colors=(int *) calloc(V*sizeof(int))))
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
        else if (i_weight < weights[jj] || (i_weight == weights[jj] && jj>i)) {
          //if the weights match, solve conflict by looking at the vertices
          //ids and taking the vertex with higher id as the max
          i_weight_is_max = 0;
          break;
        }
      }

      // if the vertex weight is a max and vertex hasnt been colored, color it
      // with the smallest color possible that is not one of neighbor_colors
      if (1 == i_weight_is_max && 0 == colors[off+i]) {
        // find smallest color to assign to the j vertex that color is either
        //  a)  1 if none of the neighbors is colored or the smallest color of
        //      a neighbor is >1
        //  b)  In between a color in the array of neighbors colors if there
        //      is a gap between two of the (sorted) neighbors colors
        //  c)  1 more than the last color in the sorted array of neighbors
        //      colors sort neighbors colors.
        qsort(neighbor_colors, num_colors, sizeof(int), compare);

        if (0 == num_colors || 1 < neighbor_colors[0]) {
          min_color =1;
        }
        else {
          for (j=0; j<num_colors; ++j) {
            if (j != num_colors-1 &&
                neighbor_colors[j+1] != neighbor_colors[j]+1) {
              min_color = neighbor_colors[j]+1;
              break;
            }
          }
          if (j == num_colors)
            min_color = neighbor_colors[num_colors-1]+1;
        }
        i_colors[i] = min_color;
      }

      free(neighbor_colors);
    }

    // todo: make this Allgatherv
    //each process sends the colors of its vertices to root
    MPI_Gatherv(i_colors, lV, MPI_INT, colors, range, vertex_offsets,
      MPI_INT, root, MPI_COMM_WORLD);

    //root synchronizes colors on all processes
    MPI_Bcast(colors, V, MPI_INT, root, MPI_COMM_WORLD);
  }
  free(i_colors);
#endif

  for (i=0;i<chromaticity_upper;i++) {
    //for each vertex in this process
    remainder_v_per_p = (V + rank) % npes;
    //index of the first vertex for process i
    first_v = num_v_per_p*rank+remainder_v_per_p*(remainder_v_per_p< rank);

     //for each vertex of this process
    for (j=0;j<range[rank];j++) {
      //get the vertex weight
      j_weight = weights[first_v+j];
      j_weight_is_max = true;
      if (NULL == (neighbor_colors=(int *) malloc(V*sizeof(int))))
        abort();
      memset(neighbor_colors, 0, V*sizeof(int));
      num_colors=0;

      //compare vertex weight to weights of its non-colored neighbors to see
      //if it is a maximum. Also gather the colors of all neighbors of the
      //vertex j that have been colored
      for (k=0; k<V; ++k) {
        //if there is an edge between j vertex and neighbor k vertex
        if (1 == p_graph[j*V+k]) {
          //if neighbor is colored just add its color to the neighbor_colors
          if (0 != colors[k]) {
            neighbor_colors[num_colors++]=colors[k];
          }
          //if the weights match, solve conflict by looking at the vertices
          //ids and taking the vertex with higher id as the max
          else if (j_weight< weights[k] || (j_weight==weights[k] && k>j)) {
            j_weight_is_max = false;
            break;
          }
        }
      }

      //if the vertex weight is a max and vertex hasnt been colored,
      //color it with the smallest color possible that is not one of
      //neighbor_colors
      if (j_weight_is_max==true && colors[first_v+j]==0) {
        /* find smallest color to assign to the j vertex that color is either
            a)  1 if none of the neighbors is colored or the smallest color of
                a neighbor is >1
            b)  In between a color in the array of neighbors colors if there
                is a gap between two of the (sorted) neighbors colors
            c)  1 more than the last color in the sorted array of neighbors
                colors sort neighbors colors. */
        qsort(neighbor_colors,num_colors,sizeof(int),compare);
        if (num_colors==0 || neighbor_colors[0]>1) {
          min_color =1;
        }
        else {
          for (k=0;k<V;k++) {
            if (k<V-1 && (neighbor_colors[k+1]-neighbor_colors[k]>1)) {
              min_color = neighbor_colors[k]+1;
              break;
            }
            else {
              min_color = neighbor_colors[num_colors-1]+1;
            }
          }
        }
        j_colors[j] = min_color;
      }
      free(neighbor_colors);
    }
    //each process sends the colors of its vertices to root
    MPI_Gatherv(j_colors, range[rank], MPI_INT, colors, range, vertex_offsets,
                MPI_INT, root, MPI_COMM_WORLD);
    //root synchronizes colors on all processes
    MPI_Bcast(colors, V, MPI_INT, root, MPI_COMM_WORLD);
  }
  free(j_colors);
}


int main(int argc, char * argv[])
{
  int i;
  //graph with edges corresponding to the
  //vertices in this process and their size. Note that if say process p
  //has vertices 0 and 1, then p_graph will be a 2xV matrix, so it has
  //all the edges between 1,2 and all vertices
  //offsets is the index in graph where to start copying to p_graph
  int * range, * vertex_offsets, * p_graph;
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

  /* only root reads the file and loads the full graph */
  if (rank == root) {
    read_graph(input_filename);
    printf("V=%d E=%d Chromaticity Upper Bound=%d\n",V,E,chromaticity_upper);
    fflush(stdout);
  }
  else {
    distribute(rank, 0);
  }

  distribute_graph(&range, &vertex_offsets, &p_graph);

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
  jones_plassmann(range, vertex_offsets, p_graph);

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

  free(input_filename);
  free(ia);
  free(ja);
  free(weights);
  free(colors);
  free(off);

  MPI_Finalize();

  return 0;
}
