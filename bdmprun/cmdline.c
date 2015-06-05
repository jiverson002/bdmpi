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
  {"nlm",   0,      0,      BDMPRUN_CMD_NOLOCKMEM},

  {"pg",    1,      0,      BDMPRUN_CMD_PGSIZE},
  {"rm",    1,      0,      BDMPRUN_CMD_RMSIZE},
  {"sbma",  1,      0,      BDMPRUN_CMD_SBMA},

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
"  -pg=int [Default: 4]",
"     Specifies the number of system pages which make a single sbpage.",
" ",
"  -rm=int [Default: 917504]",
"     Specifies the maximum resident set size for the slave processes on",
"     each node. The size is in terms of number of system pages.",
" ",
"  -wd=string [Default: "BDMPRUN_DEFAULT_WDIR"]",
"     Specifies where working files will be stored.",
" ",
"  -sbma=string [Default: "BDMPRUN_DEFAULT_SBMA"]",
"     Speicifies the operating mode of the SBMA library.  Valid options are:",
"       araw  Aggressive read / aggressive write",
"       arlw  Aggressive read / lazy write",
"       lraw  Lazy read       / aggressive write",
"       lrlw  Lazy read       / lazy write",
" ",
"     The `laziness' of the lr* methods is controlled by the system page ",
"     multiplier command line parameter `-pg='.",
"     The `laziness' of the *lw methods is controlled by the resident memory ",
"     command line parameter `-rm='.",
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
  bdmp->sbopts   = BDMPRUN_DEFAULT_SBOPTS;
  bdmp->smsize   = BDMPRUN_DEFAULT_SMSIZE;
  bdmp->imsize   = BDMPRUN_DEFAULT_IMSIZE;
  bdmp->mmsize   = BDMPRUN_DEFAULT_MMSIZE;
  bdmp->rmsize   = BDMPRUN_DEFAULT_RMSIZE;
  bdmp->pgsize   = BDMPRUN_DEFAULT_PGSIZE;
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

      case BDMPRUN_CMD_PGSIZE:
        if (gk_optarg) bdmp->pgsize = atoi(gk_optarg);
        break;

      case BDMPRUN_CMD_RMSIZE:
        if (gk_optarg) bdmp->rmsize = atoi(gk_optarg);
        break;

      case BDMPRUN_CMD_SBMA:
        if (0 == strncmp(gk_optarg, "araw", 5)) {
        }
        else if (0 == strncmp(gk_optarg, "arlw", 5)) {
          bdmp->sbopts &= ~BDMPI_SB_SAVEALL;
          bdmp->sbopts |= BDMPI_SB_LAZYWRITE;
        }
        else if (0 == strncmp(gk_optarg, "lraw", 5)) {
          bdmp->sbopts |= BDMPI_SB_LAZYREAD;
        }
        else if (0 == strncmp(gk_optarg, "lrlw", 5)) {
          bdmp->sbopts &= ~BDMPI_SB_SAVEALL;
          bdmp->sbopts |= BDMPI_SB_LAZYREAD|BDMPI_SB_LAZYWRITE;
        }
        else {
          printf("invalid sbma option `%s'\n", gk_optarg);
        }
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
  if (bdmp->pgsize < 0)
    errexit("The value for -pg should be non-negative.");
  if (bdmp->rmsize < 0)
    errexit("The value for -rm should be non-negative.");

  return bdmp;
}
