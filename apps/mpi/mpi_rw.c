/*!
\file
\brief A parallel random walk program
\date Started 5/6/2013
\author Huzefa Rangwala
*/


#include <GKlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>


/*** RW.C in GKLib *******************/

/*************************************************************************/
/*! Data structures for the code */
/*************************************************************************/
typedef struct {
  int niter;
  int ntvs;
  int ppr;
  float eps;
  float lamda;
  char *infile;
  char *outfile;
} params_t;

/*************************************************************************/
/*! Constants */
/*************************************************************************/
#define CMD_NITER       1
#define CMD_EPS         2
#define CMD_LAMDA       3
#define CMD_PPR         4
#define CMD_NTVS        5
#define CMD_HELP        10


/*************************************************************************/
/*! Local variables */
/*************************************************************************/
static struct gk_option long_options[] = {
  {"niter",      1,      0,      CMD_NITER},
  {"lamda",      1,      0,      CMD_LAMDA},
  {"eps",        1,      0,      CMD_EPS},
  {"ppr",        1,      0,      CMD_PPR},
  {"ntvs",       1,      0,      CMD_NTVS},
  {"help",       0,      0,      CMD_HELP},
  {0,            0,      0,      0}
};


/*-------------------------------------------------------------------*/
/* Mini help  */
/*-------------------------------------------------------------------*/
static char helpstr[][100] = {
" ",
"Usage: rw [options] <graph-file> <out-file>",
" ",
" Required parameters",
"  graph-file",
"     The name of the file storing the transactions. The file is in ",
"     Metis' graph format.",
" ",
" Optional parameters",
"  -niter=int",
"     Specifies the maximum number of iterations. [default: 100]",
" ",
"  -lamda=float",
"     Specifies the follow-the-adjacent-links probability. [default: 0.80]",
" ",
"  -eps=float",
"     Specifies the error tollerance. [default: 1e-10]",
" ",
"  -ppr=int",
"     Specifies the source of the personalized PR. [default: -1]",
" ",
"  -ntvs=int",
"     Specifies the number of test-vectors to compute. [default: -1]",
" ",
"  -help",
"     Prints this message.",
""
};

static char shorthelpstr[][100] = {
" ",
"   Usage: rw [options] <graph-file> <out-file>",
"          use 'rw -help' for a summary of the options.",
""
};
 


/*************************************************************************/
/*! Function prototypes */
/*************************************************************************/
void print_init_info(params_t *params, gk_csr_t *mat);
void print_final_info(params_t *params);
params_t *parse_cmdline(int argc, char *argv[]);

int mpi_rw_PageRank (params_t *params, float lamda, float eps, int max_niter, MPI_Comm icomm)
{
	int npes, mype, mystatus, status;
	int npp;
	double *rscale, *prnew, *prold, *prtmp;
	int iter;

	int i, j, k, ntotal, nrecv, nrecv_min, nrecv_max, nlocal;
	MPI_Comm comm;
	int nodes_per_processor; //avg. nodes per processor
	int *edge_elem_per_processor;
	int *tmp_nodes_pp, *nodes_per_p;
	int *start_outgoing;
	int new_processor_flag;
	int niter;	
	float *pr;
	int ee;
	int *mymatrowptr;
	int *mymatind;
	int *displ, *displ_per_processor;
	int *node_owner; // which processor own which node;
	gk_csr_t *mat;
	int total_nodes;
	int *outgoing_linked_nodes;
	
	ssize_t *rowptr;
	int *rowind;
	float *rowval;
	
	
	rowptr = NULL;
	rowind = NULL;

	MPI_Comm_dup (icomm, &comm);

	// Split the matrix by npes - number of processor.

	MPI_Comm_size (comm, &npes);
	MPI_Comm_rank (comm, &mype);
	
	
	edge_elem_per_processor = gk_ismalloc (npes, 0,  "Setting the number of elements per processor");
	nodes_per_p = gk_ismalloc (npes, 0, "nodes per processor");
	tmp_nodes_pp = gk_ismalloc (npes, 0, "nodes per processor");
	displ_per_processor = gk_ismalloc (npes, 0, "displ_processor");

	start_outgoing = gk_imalloc (npes,  "mpi_rw_PageRank");
	
	printf ("Debug: Processor %d: %d\n",mype, npes);
	if (mype == 0){
	//Send matrix to other processors;
	
	int processor_index; 
	
 	mat = gk_csr_Read(params->infile, GK_CSR_FMT_METIS, 1, 1);
	
	node_owner = gk_ismalloc (mat->nrows,0,  "mpi_rw");
	total_nodes = mat->nrows;
	//nodes_per_processor =  (int) ((1.0*mat->nrows)/npes+0.5);
        nodes_per_processor = (int) ceil ((mat->nrows *1.0)/npes); 
	printf ("Debug %d %d NPP: %d\n", mat->nrows,npes,nodes_per_processor);	

	displ = gk_ismalloc (npes, 0, "mpi:pagerank");

	print_init_info (params, mat);
	
	rowptr = mat->rowptr;
	rowind = mat->rowind;
	rowval = mat->rowval;

	// mat->rowind has all the outgoing edges;
	// mat->rowptr has all info. on the 
	// mat->rowval has the values
	for (i=0;i<=mat->nrows;i++)
		printf ("%zd ", rowptr[i]);
	printf ("\n");
	for (i=0;i<=mat->nrows;i++){
	// printf ("%d %d\n", rowptr[i], rowptr[i+1]);
	 for (j=rowptr[i]; j<rowptr[i+1];j++){
	 	printf ("%d %d\n", i+1, rowind[j]+1);
	 }
	 
	}	


	new_processor_flag = -1;
	displ_per_processor[0] = 0;
	for(i=0;i<mat->nrows;i++){
	 processor_index = (int) i / nodes_per_processor;
	 edge_elem_per_processor[processor_index] +=  (rowptr[i+1] - rowptr[i]); // This determines the size of the outgoing edges
	 node_owner [i] = processor_index;
	 nodes_per_p [processor_index] ++;
	}
	
	displ[0] = 0;
	for (i=1;i<npes;i++){
		displ[i] = displ[i-1]+ edge_elem_per_processor[i-1];
		displ_per_processor[i] = displ_per_processor[i-1] + nodes_per_p[i-1];
	}
	for (i=0;i<npes;i++)	
		printf ("NPP: %d %d %d\n", i, nodes_per_p[i], displ_per_processor[i]);
	//Broadcast information to all processors about size
	
	// * Compute the outgoing edges for every processor
	// Store this in an array
	// Broadcast this array
	// All_to_all_v for the elements
	
	
	}	
	
	MPI_Bcast (&nodes_per_processor, 1, MPI_INT, 0, comm);

	MPI_Barrier (comm);
	printf ("Check Processor %d\n", mype);
	MPI_Scatter (edge_elem_per_processor, 1, MPI_INT, &ee, 1, MPI_INT,  0, comm); 
	MPI_Bcast (&total_nodes, 1, MPI_INT, 0, comm);
	MPI_Bcast (nodes_per_p, npes, MPI_INT, 0, comm);
	for (i=0; i<npes; i++)
		tmp_nodes_pp[i] = 1 + nodes_per_p[i];
	
	MPI_Barrier (comm);	
	
	if (mype!=0){
		node_owner = gk_imalloc (total_nodes, "mpi_rw");
		displ_per_processor = gk_imalloc (npes,"mpi_rw");
		displ = gk_imalloc (npes, "mpi_rw");
		
	}
	printf ("Total Nodes: %d, %d\n", mype, total_nodes);

	MPI_Bcast (node_owner, total_nodes, MPI_INT, 0, comm);
	//MPI_Bcast (displ_per_processor, npes, MPI_INT, 0, comm);
	//MPI_Bcast (edge_elem_per_processor, npes, MPI_INT, 0, comm);
	//MPI_Bcast (displ, npes, MPI_INT, 0, comm);
	MPI_Barrier (comm);	

	printf ("PR#%d: EE%d\n", mype, ee);	
	mymatind = gk_imalloc (ee, "mymatind");
	mymatrowptr = gk_imalloc (nodes_per_p[mype]+1, "mymatrowptr");

	//Source for bug.
	//for (i=0;i<npes;i++)
	//	printf("DBG: %d %d\n",tmp_nodes_pp[i], displ_per_processor[i]);
	MPI_Barrier (comm);
	MPI_Scatterv(rowptr, tmp_nodes_pp, displ_per_processor, MPI_INT, mymatrowptr, tmp_nodes_pp[mype], MPI_INT, 0, comm);
	MPI_Scatterv(rowind, edge_elem_per_processor, displ, MPI_INT, mymatind, ee, MPI_INT, 0, comm); 

	MPI_Barrier(comm);

	//printf ("Scatter Done \n");

	//mymatrowptr[nodes_per_p[mype]] = ee;
		
	for (i=0;i<ee;i++)
	 printf ("P%d: %d\n", mype,mymatind[i]+1);

//	for (i=0;i<nodes_per_p[mype];i++){
//	printf ("Node: %d\n",i);
//	for (j=mymatrowptr[i];j<mymatrowptr[i+1];j++)
//	printf ("PrXXX%d : Node;%d :: %d\n",mype, i+1, mymatind[j-mymatrowptr[i]]+1);
//
//		printf ("Processor %d: RowPtr %d\n", mype, mymatrowptr[i]); 
//	}	

	
	//Initialize the probabilities;
 	
	npp = nodes_per_p[mype];
	printf ("Debgugging for Processor %d, Max_iter %d, NPP %d\n", mype, max_niter, npp);

	prold  = gk_dsmalloc (npp, 0.0, "prnew");
	prnew  = gk_dsmalloc (npp, 0.0, "ornew");
	rscale = gk_dsmalloc (npp, 0.0, "rscale");
	/* compute the scalling factors to get adjacency weights into
		transition matrix */
	for (i=0;i<npp;i++)
	  rscale[i] = 1.0/npes;
	// Assuming no rowval;
	
	// Do the priteration //
	for(iter=0; iter<max_niter&& iter<10;iter++){
	  MPI_Barrier (comm); 
	  gk_SWAP (prnew, prold, prtmp);
	  gk_dset (npp, 0.0, prnew);

	  // push random-walk scores to outlinks;
	  
	 outgoing_linked_nodes = gk_ismalloc (edge_elem_per_processor[mype], 0, "outgoint");
	 
 
	for (i=0;i<npp;i++){
	//   printf ("Check PROCESSOR#%d RWPTR[%d]=%d\n", mype,i, mymatrowptr[i]);
	
// mymatrowptr[i+1]);
	   for (j=mymatrowptr[i];j<mymatrowptr[i+1];j++){
		// determine index of node to send to
		// determine which processor has it
		// determine what needs to be sent; prold[i], rscale[i]
	//	outgoing_linked_nodes[j-mymatrowptr[i]] = node_owner[mymatind[j-mymatrowptr[i]]];
		printf ("Processor[%d] Node-From: %d Node-To %d in Processor[%d]j=%d mymatrowptr[i]=%d\n",mype,i+mype*nodes_per_processor, mymatind[j-mymatrowptr[0]], node_owner[mymatind[j-mymatrowptr[0]]],j,mymatrowptr[0]);		
				 
	
            
	    }
	}
		




	gk_free ((void **) &outgoing_linked_nodes, LTERM);


	}
	

	niter = 10;
	MPI_Barrier (comm);
	return niter;




}


/**************************************************************************/
/**************************************************************************/
int main(int argc, char **argv)
{
  int k, npes, mype, lastelmnt;
  int nlocal, nremote;
  int *elmnts;
  MPI_Status status;

  ssize_t i, j, niter;
  params_t *params;
  FILE *fpout;
  float *probs;


  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);


  MPI_Barrier (MPI_COMM_WORLD);


  
  /* Read in all the parameters */
 
  /* Read in the graph file */

 //if (mype == 0){

  /* get command-line options */

  params = parse_cmdline(argc, argv);

  /* read the data */
 // mat = gk_csr_Read(params->infile, GK_CSR_FMT_METIS, 1, 1);

  /* display some basic stats */
  //print_init_info(params, mat);
 //} 

 //MPI_Barrier (MPI_COMM_WORLD);

 printf ("Beginning PageRank: %d\n", mype);
 //probs = gk_fsmalloc(mat->nrows, 0.0, "main:probs");
 
 niter = mpi_rw_PageRank(params, params->lamda,params->eps,params->niter,MPI_COMM_WORLD); 
  
  /* Output the probability vector */
  
  //fpout = gk_fopen (params->outfile, "w", "main: outfile");
//  for (i=0;i<mat->nrows;i++)
//	fprintf (fpout, "%.4e\n", probs[i]);
 // gk_fclose (fpout);
  
  //gk_free ((void *) &probs, LTERM);

   


 // if (mype == 0)	
 // 	gk_csr_Free(&mat);
  MPI_Finalize ();

  return 0;
}


/*************************************************************************/
/*! This function prints run parameters */
/*************************************************************************/
void print_init_info(params_t *params, gk_csr_t *mat)
{
  printf("*******************************************************************************\n");
  printf(" fis\n\n");
  printf("Matrix Information ---------------------------------------------------------\n");
  printf(" input file=%s, [%d, %d, %zd]\n", 
      params->infile, mat->nrows, mat->ncols, mat->rowptr[mat->nrows]);

  printf("\n");
  printf("Options --------------------------------------------------------------------\n");
  printf(" niter=%d, ntvs=%d, ppr=%d, lamda=%f, eps=%e\n",
      params->niter, params->ntvs, params->ppr, params->lamda, params->eps);

  printf("\n");
  printf("Performing random walks... ----------------------------------------------\n");
}


/*************************************************************************/
/*! This function prints final statistics */
/*************************************************************************/
void print_final_info(params_t *params)
{
  printf("\n");
  printf("Memory Usage Information -----------------------------------------------------\n");
  printf("   Maximum memory used:              %10zd bytes\n", (ssize_t) gk_GetMaxMemoryUsed());
  printf("   Current memory used:              %10zd bytes\n", (ssize_t) gk_GetCurMemoryUsed());
  printf("********************************************************************************\n");
}


/*************************************************************************/
/*! This is the entry point of the command-line argument parser */
/*************************************************************************/
params_t *parse_cmdline(int argc, char *argv[])
{
  int i;
  int c, option_index;
  params_t *params;

  params = (params_t *)gk_malloc(sizeof(params_t), "parse_cmdline: params");

  /* initialize the params data structure */
  params->niter     = 100;
  params->ppr       = -1;
  params->ntvs      = -1;
  params->eps       = 1e-10;
  params->lamda     = 0.80;
  params->infile    = NULL;
  params->outfile   = NULL;


  /* Parse the command line arguments  */
  while ((c = gk_getopt_long_only(argc, argv, "", long_options, &option_index)) != -1) {
    switch (c) {
      case CMD_NITER:
        if (gk_optarg) params->niter = atoi(gk_optarg);
        break;
      case CMD_NTVS:
        if (gk_optarg) params->ntvs = atoi(gk_optarg);
        break;
      case CMD_PPR:
        if (gk_optarg) params->ppr = atoi(gk_optarg);
        break;
      case CMD_EPS:
        if (gk_optarg) params->eps = atof(gk_optarg);
        break;
      case CMD_LAMDA:
        if (gk_optarg) params->lamda = atof(gk_optarg);
        break;

      case CMD_HELP:
        for (i=0; strlen(helpstr[i]) > 0; i++)
          printf("%s\n", helpstr[i]);
        exit(0);
        break;
      case '?':
      default:
        printf("Illegal command-line option(s)\nUse %s -help for a summary of the options.\n", argv[0]);
        exit(0);
    }
  }

  if (argc-gk_optind != 2) {
    printf("Unrecognized parameters.");
    for (i=0; strlen(shorthelpstr[i]) > 0; i++)
      printf("%s\n", shorthelpstr[i]);
    exit(0);
  }

  params->infile  = gk_strdup(argv[gk_optind++]);
  params->outfile = gk_strdup(argv[gk_optind++]);

  if (!gk_fexists(params->infile))
    errexit("input file %s does not exist.\n", params->infile);

  if (params->ppr != -1 && params->ntvs != -1)
    errexit("Only one of the -ppr and -ntvs options can be specified.\n");

  return params;
}


