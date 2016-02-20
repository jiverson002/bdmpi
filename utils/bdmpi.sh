#! /bin/bash


die() {
#{{{1
  #RET=$? ; echo "`date '+%s'` ${1} returned ${RET}" 1>&2 ; exit ${RET}
  RET=$? ; echo "`date '+%s'` ${1} returned ${RET}" 1>&2
#}}}1
}


#------------------------------------------------------------------------------
# CONFIGURATION VARIABLES
#------------------------------------------------------------------------------
#{{{1
# Output directory
OUT_DIR="${HOME}/local/results/bdmpi/hdd"

# Error file
ERR_FILE="${OUT_DIR}/ERR"

# Working directory for bdmpi runtime
WRK_DIR="/scratch/bdmp"

# Bitmask of experiments to run
BITMASK=127

# Program used to run experiment
PROG=pr1d

# Number of simultaneously running slaves
NR=1

# Number of trials for each experiment
NT=1

# Memory threshold for sbma runtime
MEM=910000

# Page size for sbma runtime
PAGE=4

# Debugging indicator variable
DEBUG=0
#}}}1


#------------------------------------------------------------------------------
# COMMAND-LINE PARSING
#------------------------------------------------------------------------------
#{{{1
while [[ $# > 0 ]] ; do
#{{{2
  key=$1
  case $key in
  -b)
    shift # past argument
    key="$1"
    ;&
  --bitmask=*)
    BITMASK="${key#*=}"
    ;;
  -m)
    shift # past argument
    key="$1"
    ;&
  --memory=*)
    MEM="${key#*=}"
    ;;
  -s)
    shift # past argument
    key="$1"
    ;&
  --pagesize=*)
    PAGE="${key#*=}"
    ;;
  -r)
    shift # past argument
    key="$1"
    ;&
  --nr=*)
    NR="${key#*=}"
    ;;
  -t)
    shift # past argument
    key="$1"
    ;&
  --trials=*)
    NT="${key#*=}"
    ;;
  -d|--debug)
    DEBUG=1
    ;;
  -p)
    shift # past argument
    key="$1"
    ;&
  --prog=*)
    PROG="${key#*=}"
    ;;
  -h|--help)
    echo "usage: bdmpi.sh"
    echo "  -b | --bitmask="
    echo "  -m | --memory="
    echo "  -s | --pagesize="
    echo "  -r | --nr="
    echo "  -t | --trials="
    echo "  -p | --prog="
    echo "  -d | --debug"
    echo "  -h | --help"
    exit
    ;;
  *)
    # unknown option
    echo "unknown option '$key'"
  ;;
  esac
  shift # past value or argument=value or argument with no value
#}}}2
done
#}}}1


#------------------------------------------------------------------------------
# SBMA CONFIG INFO
#------------------------------------------------------------------------------
#{{{1
SBMA_DIR="/scratch/jeremy/src/sbma"
SBMA_HASH="`git -C ${SBMA_DIR} log -n 1 | head -n 1 | awk '{print substr($2,0,7)}'`"
#}}}1


#------------------------------------------------------------------------------
# BDMPI CONFIG INFO
#------------------------------------------------------------------------------
#{{{1
# Path to bdmprun executable
BDMP_CMD=/scratch/jeremy/bin/bdmprun

# Arguments passed to bdmprun executable
BDMP_ARG="-ns=\${NS} -nr=\${NR} -pg=\${PAGE} -rm=\${MEM} -sbma=\${SBMA[\${i}]} -wd=${WRK_DIR}"

# Source directory for bdmpi source code
BDMPI_DIR="/scratch/jeremy/src/bdmpi"

# Git commit of bdmpi source code
BDMPI_HASH="`git -C ${BDMPI_DIR} log -n 1 | head -n 1 | awk '{print substr($2,0,7)}'`"
#}}}1


#------------------------------------------------------------------------------
# PROGRAM VARIABLES
#------------------------------------------------------------------------------
#{{{1
  declare -A PROG_NS   # number of slaves
  declare -A PROG_CMD  # path to program executable
  declare -A PROG_ARG  # arguments to be passed to program executable
  declare -A SBMA      # sbma operating modes

  SBMA=([1]="osvmm" [2]="aggrd -sa" [4]="aggrd " [8]="lzyrd -sa" [16]="lzyrd")

  #============================================================================
  # MATRIX-FACTORIZATION
  #============================================================================
  #{{{2
  #BDMP_CMD="mpirun -n 4 -f /home/jeremy/local/host ${BDMP_CMD}"
  PROG_NS[mf1d2]=16
  PROG_CMD[mf1d2]="/scratch/jeremy/bin/bdmp_mf1d2"
  PROG_ARG[mf1d2]="/data/graphs/nlpkkt200.bcsr 5 150"
  #}}}2

  #============================================================================
  # PAGERANK
  #============================================================================
  #{{{2
  PROG_NS[pr1d]=14
  PROG_CMD[pr1d]="/scratch/jeremy/bin/bdmp_pr1d"
  PROG_ARG[pr1d]="/scratch-ssd/uk2007.bcsr-shuff 5 1 0"
  #}}}2

  #============================================================================
  # KMEANS
  #============================================================================
  #{{{2
  PROG_NS[sphkmeans]=16
  PROG_CMD[sphkmeans]="/scratch/jeremy/bin/bdmp_sphkmeans"
  PROG_ARG[sphkmeans]="/data/kmeans/new1600k.bcsr 150 1 5 1 0"
  #}}}2

  #============================================================================
  # PARMETIS
  #============================================================================
  #{{{2
  PROG_NS[parmetis]=6
  PROG_CMD[parmetis]="/scratch/jeremy/bin/parmetis"
  PROG_ARG[parmetis]="/data/graphs/nlpkkt240.graph 1 \${NS} 0 0 7 0"
  #}}}2

  #============================================================================
  # SPLATT
  #============================================================================
  #{{{2
  PROG_NS[splatt]=12
  PROG_CMD[splatt]="/scratch/jeremy/bin/splatt"
  PROG_ARG[splatt]="cpd --tol=0 -r 16 -d 4x3x2 --nowrite -i 5 \
                    /scratch-ssd/tensors/nell_large_rand.fixed.tns"
  #}}}2
#}}}1


#------------------------------------------------------------------------------
# INPUT VALIDATION
#------------------------------------------------------------------------------
#{{{1
NS=${PROG_NS[$PROG]}
CMD=${PROG_CMD[$PROG]}
ARG=${PROG_ARG[$PROG]}

case $NR in
#{{{2
[[:digit:]]*)
  if [[ $NR -gt $NS ]] ; then
    NR=$NS
  fi
  ;;
ns)
  NR=$NS
  ;;
*)
  ;;
#}}}2
esac

# Template log name, parameterized on input
TEMPLATE="${OUT_DIR}/${SBMA_HASH}-${BDMPI_HASH}-${PROG}-${NS}-${NR}"
#}}}1


#------------------------------------------------------------------------------
# ERROR LOG CONFIGURATION
#------------------------------------------------------------------------------
#{{{1
if [[ $DEBUG -ne 1 ]] ; then
  # Touch error file
  touch "${ERR_FILE}"
  # Close STDERR FD
  exec 2<&-
  # Open STDERR FD as $ERR_FILE file for write.
  exec 2>${ERR_FILE}
fi
#}}}1


#------------------------------------------------------------------------------
# RUN EXPERIMENTS
#------------------------------------------------------------------------------
#{{{1
if [[ $DEBUG -eq 1 ]] ; then
  echo "BITMASK  = $BITMASK"
  echo "TEMPLATE = $TEMPLATE"
  echo ""
fi

for i in 1 2 4 8 16 ; do
#{{{2
  for t in `seq 1 ${NT}` ; do
  #{{{3
    if (( ${i} == (${BITMASK} & ${i}) )) ; then
    #{{{4
      # Increment the identifier in the filename until the filename is unique
      LOG="${TEMPLATE}-${SBMA[${i}]// /}-0"
      while [ -e ${LOG} ] ; do
      #{{{5
        ID=${LOG##*-}
        ID=$((${ID}+1))
        LOG="${TEMPLATE}-${SBMA[${i}]// /}-${ID}"
      #}}}5
      done

      if [[ $DEBUG -eq 1 ]] ; then
      #{{{5
        eval "echo ${BDMP_CMD} ${BDMP_ARG} ${CMD} ${ARG} \&\>\> ${LOG}"
      #}}}5
      else
      #{{{5
        rm -rf ${WRK_DIR}/* /dev/mqueue/* /dev/shm/*
        sudo fstrim -v /scratch-ssd

        echo "${LOG}"
        eval "echo ${BDMP_CMD} ${BDMP_ARG} ${CMD} ${ARG} > ${LOG}"
        (eval "${BDMP_CMD} ${BDMP_ARG} ${CMD} ${ARG} &>> ${LOG}") \
          || die ${LOG}
      #}}}5
      fi
    #}}}4
    fi
  #}}}3
  done
#}}}2
done
#}}}1
