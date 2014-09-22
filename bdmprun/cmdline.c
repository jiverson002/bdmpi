/*!
\file cmdline.c
\brief Command-line argument parsing 

\date 3/31/2013
\author George
\version\verbatim $Id: cmdline_gpmetis.c 13901 2013-03-24 16:17:03Z karypis $\endverbatim
*/

#include "bdmprun.h"


/*-------------------------------------------------------------------
 * Command-line options 
 *-------------------------------------------------------------------*/
static struct gk_option long_options[] = {
  {"ns",    1,      0,      BDMPRUN_CMD_NS},
  {"nr",    1,      0,      BDMPRUN_CMD_NR}, 
  {"wd",    1,      0,      BDMPRUN_CMD_WDIR}, 
  {"sm",    1,      0,      BDMPRUN_CMD_SMSIZE}, 
  {"im",    1,      0,      BDMPRUN_CMD_IMSIZE}, 
  {"mm",    1,      0,      BDMPRUN_CMD_MMSIZE}, 
  {"sb",    1,      0,      BDMPRUN_CMD_SBSIZE}, 
  {"nlm",   0,      0,      BDMPRUN_CMD_NOLOCKMEM}, 

  {"dl",    1,      0,      BDMPRUN_CMD_DBGLVL},
  {"h",     0,      0,      BDMPRUN_CMD_HELP},
  {0,       0,      0,      0}
};


/*-------------------------------------------------------------------
 * Mini help
 *-------------------------------------------------------------------*/
static char helpstr[][100] =
{
" ",
"Usage: bdmprun [options] exefile [options for the exe-file]",
" ",
" Required parameters",
"    exefile     The program to be executed.",
" ",
" Optional parameters",
"  -ns=int [Default: 1]",
"     Specifies the number of slave processes on each node.",
" ",
"  -nr=int [Default: 1]",
"     Specifies the maximum number of concurrently running slaves.",
" ",
"  -nc=int [Default: 1]",
"     Specifies the maximum number of slaves in a critical section.",
" ",
"  -sm=int [Default: 20]",
"     Specifies the number of shared memory pages allocated for each slave.",
" ",
"  -im=int [Default: 4]",
"     Specifies the maximum size of a message that will be buffered",
"     in the memory of the master. Messages longer than that are buffered",
"     on disk. The size is in terms of memory pages.",
" ",
"  -mm=int [Default: 32]",
"     Specifies the maximum size of the buffer to be used by MPI for ",
"     inter-node communication. The size is in terms of memory pages.",
/*
" ",
"  -dm=int [Default: 8]",
"     Specifies the minimum size of the message that will lead to disk-based",
"     buffering. Disk-based buffering can happen for smaller messages due to",
"     the -im limit.",
" ",
"  -nlm",
"     Specifies that the master will not lock in memory data used for ",
"     bcast/reduce operations.",
*/
" ",
"  -sb=int [Default: 32]",
"     Specifies the size of allocations for which the explicit storage backed",
"     subsystem should be used. The size is in terms of number of pages and a",
"     value of 0 turns it off.",
" ",
"  -wd=string [Default: "BDMPRUN_DEFAULT_WDIR"]",
"     Specifies where working files will be stored.",
" ",
"  -dl=int [Default: 0]",
"     Selects the dbglvl.",
" ",
"  -h",
"     Prints this message.",
""
};


/*************************************************************************
* This is the entry point of the command-line argument parser
**************************************************************************/
mjob_t *parse_cmdline(int argc, char *argv[])
{
  int i, nargs;
  int c, option_index;
  mjob_t *bdmp;

  bdmp = (mjob_t *)gk_malloc(sizeof(mjob_t), "parse_cmdline");
  memset((void *)bdmp, 0, sizeof(mjob_t));


  /* initialize with default values */
  bdmp->ns       = BDMPRUN_DEFAULT_NS;
  bdmp->nr_input = BDMPRUN_DEFAULT_NR;
  bdmp->nc       = BDMPRUN_DEFAULT_NC;
  bdmp->smsize   = BDMPRUN_DEFAULT_SMSIZE;
  bdmp->imsize   = BDMPRUN_DEFAULT_IMSIZE;
  bdmp->mmsize   = BDMPRUN_DEFAULT_MMSIZE;
  bdmp->sbsize   = BDMPRUN_DEFAULT_SBSIZE;
  bdmp->lockmem  = BDMPRUN_DEFAULT_LOCKMEM;
  bdmp->dbglvl   = BDMPRUN_DEFAULT_DBGLVL;
  bdmp->iwdir    = NULL;
  bdmp->exeargv  = NULL;



  /* Parse the command line arguments  */
  while ((c = gk_getopt_long_only(argc, argv, "+", long_options, &option_index)) != -1) {
    switch (c) {
      case BDMPRUN_CMD_NS:
        if (gk_optarg) bdmp->ns = atoi(gk_optarg);
        break;

      case BDMPRUN_CMD_NR:
        if (gk_optarg) bdmp->nr_input = atoi(gk_optarg);
        break;

      case BDMPRUN_CMD_SMSIZE:
        if (gk_optarg) bdmp->smsize = (size_t)atoi(gk_optarg);
        break;

      case BDMPRUN_CMD_IMSIZE:
        if (gk_optarg) bdmp->imsize = (size_t)atoi(gk_optarg);
        break;

      case BDMPRUN_CMD_MMSIZE:
        if (gk_optarg) bdmp->mmsize = (size_t)atoi(gk_optarg);
        break;

      case BDMPRUN_CMD_SBSIZE:
        if (gk_optarg) bdmp->sbsize = (size_t)atoi(gk_optarg);
        break;

      case BDMPRUN_CMD_DBGLVL:
        if (gk_optarg) bdmp->dbglvl = atoi(gk_optarg);
        break;

      case BDMPRUN_CMD_WDIR:
        if (gk_optarg) bdmp->iwdir = gk_strdup(gk_optarg);
        break;

      case BDMPRUN_CMD_NOLOCKMEM:
        bdmp->lockmem = 0;
        break;

      case BDMPRUN_CMD_HELP:
        for (i=0; strlen(helpstr[i]) > 0; i++)
          printf("%s\n", helpstr[i]);
        exit(0);
        break;

      case '?':
      default:
        break;
    }
  }

  if (bdmp->iwdir == NULL)
    bdmp->iwdir = gk_strdup(BDMPRUN_DEFAULT_WDIR);


  nargs = argc-gk_optind;
  if (nargs == 0) {
    printf("Missing parameters.");
    for (i=0; strlen(helpstr[i]) > 0; i++)
      printf("%s\n", helpstr[i]);
    exit(0);
  }

  /* allocate memory for bdmp->exeargv[] */
  bdmp->exeargv = (char **)gk_malloc(sizeof(char *)*(nargs+1), "exeargv");

  for (i=0; i<nargs; i++) 
    bdmp->exeargv[i] = gk_strdup(argv[gk_optind++]);
  bdmp->exeargv[nargs] = NULL;

  bdmp->exefile = gk_strdup(bdmp->exeargv[0]);

  /* check for reasonable values */
  if (bdmp->ns < 1) 
    errexit("The value for -ns should be greater than 0.");
  if (bdmp->nr_input < 1) 
    errexit("The value for -nr should be greater than 0.");
  if (bdmp->nr_input > bdmp->ns) 
    errexit("The value for -nr should be less than or equal to that of -ns.");
  if (bdmp->nc < 1) 
    errexit("The value for -nc should be greater than 0.");
  if (bdmp->nc > SEM_VALUE_MAX) 
    errexit("The value for -nc should be less than %d.", (int)SEM_VALUE_MAX);
  if (bdmp->smsize < 1) 
    errexit("The value for -sm should be greater than 0.");
  if (bdmp->imsize < 1) 
    errexit("The value for -im should be greater than 0.");
  if (bdmp->mmsize < 1) 
    errexit("The value for -mm should be greater than 0.");
  if (bdmp->sbsize < 0) 
    errexit("The value for -sb should be non-negative.");

  return bdmp;
}


