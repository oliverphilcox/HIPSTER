#!/bin/bash
# Bash wrapper for the HIPSTER code developed by Oliver Philcox & Daniel Eisenstein

# Define Inputs
SHORT=h
LONG=dat:,l_max:,R0:,k_bin:,nthreads:,string:,subsample:,f_rand:

# read the options
OPTS=$(getopt --options $SHORT --long $LONG --name "$0" -- "$@")

if [ $? != 0 ] ; then echo "Failed to parse options...exiting." >&2 ;
exit 1 ; fi

eval set -- "$OPTS"

 # set initial values
PARAM_COUNT=0
NTHREADS=10
SUBSAMPLE=1
STRING=hipster

# Help dialogue
function usageText ()
{
    echo "USAGE"
    echo "-----"
    echo "-h: Display this message"
    echo "--dat: Data file (from simulation) in (x,y,z) co-ordinates"
    echo "--l_max: Maximum Legendre multipole"
    echo "--R0: Pair count truncation radius"
    echo "--k_bin: k-space binning file"
    echo "--f_rand: Ratio of random particles to data (used for DDR) counts. This should be order a few."
    echo "--string (Optional): Identification string for output file names."
    echo "--nthreads (Optional): Number of CPU threads on which to run. Default: 10"
    echo "--subsample (Optional): Factor by which to sub-sample the data. Default: 1 (no subsampling)"
    echo
}

# extract options and their arguments into variables.

while true ; do
  case "$1" in
    -h )
      usageText;
      exit 1;
      ;;
    --dat )
      DATA="$2"
      shift 2
      let PARAM_COUNT++
      ;;
    --l_max )
      MAX_L="$2"
      shift 2
      let PARAM_COUNT++
      ;;
    --R0 )
      R0="$2"
      shift 2
      let PARAM_COUNT++
      ;;
    --k_bin )
      BINFILE="$2"
      shift 2
      let PARAM_COUNT++
      ;;
    --f_rand )
      FRAND="$2"
      shift 2
      let PARAM_COUNT++
      ;;
    --string )
      STRING="$2"
      shift 2
      ;;
    --nthreads )
      NTHREADS="$2"
      shift 2
      ;;
    --subsample )
      SUBSAMPLE="$2"
      shift 2
      ;;
    -- )
      shift
      break
      ;;
    *)
      echo "Internal error!"
      usageText;
      exit 1
      ;;
  esac
done

# Check all parameters are specified
if [ "$PARAM_COUNT" -ne 5 ]; then echo; echo 'Not all command line parameters specified!'; echo; usageText; exit 1; fi


# Print the variables
echo
echo "INPUT PARAMETERS"
echo "----------------"
echo "Data file: $DATA"
echo "Maximum Legendre multipole: $MAX_L"
echo "Pair count truncation radius: $R0"
echo "k-space binning file: $BINFILE"
echo "Random-to-data ratio: $FRAND"
echo "Output string: $STRING"
echo "CPU-threads: $NTHREADS"
echo "Subsampling: $SUBSAMPLE"
echo

# Check some variables:

if [ "$MAX_L" -ge 7 ]; then echo "Only multipoles up to ell = 6 currently implemented. Exiting;"; exit 1; fi;
if [ $((MAX_L%2)) -eq 1 ]; then echo "Maximum Legendre multipole must be even. Exiting;"; exit 1; fi;

if [ "$FRAND" -ge 40 ]; then echo "Sampling with $FRAND times more randoms than data will be very slow. Exiting;"; exit 1; fi;
if [ "$FRAND" -lt 1 ]; then echo "Should have at least as many randoms as data points. Exiting;"; exit 1; fi;

if ! ( test -f "$DATA" ); then
    echo "Data file: $DATA does not exist. Exiting;"; exit 1;
fi
if ! ( test -f "$BINFILE" ); then
    echo "Binning file: $BINFILE does not exist. Exiting;"; exit 1;
fi
if [[ $(bc <<< "$SUBSAMPLE < 1.") -eq 1 ]]; then echo "Subsampling parameter must be greater than or equal to unity. Exiting;"; exit 1; fi;

# Work out where the code is installed
CODE_DIR=`dirname "$0"`

# Compute number of k bins in file
K_BINS=`wc -l < $BINFILE`
let K_BINS++

# Subsample the data if necessary
if [[ $(bc <<< "$SUBSAMPLE > 1.") -eq 1 ]]; then
  # Count number of particles in file
  N_GAL=`wc -l < $DATA`
  let N_GAL++

  # Compute number after subsampling and ensure it's an integer
  N_SUB=`echo "$N_GAL/$SUBSAMPLE" | bc`
  N_SUB=`echo "($N_SUB+0.5)/1" | bc`

  # Run subsampling script and create new file name
  NEW_DATA=$DATA.sub
  echo "Subsampling data with subsampling ratio $SUBSAMPLE"
  python $CODE_DIR/python/take_subset_of_particles.py $DATA $NEW_DATA $N_SUB
  DATA=$NEW_DATA

  echo "Using $N_SUB particles in $DATA, from $N_GAL particles originally"
fi

# Define file names
OUTPUT_FILE=$CODE_DIR/output/${STRING}_bispectrum_n${K_BINS}_l${MAX_L}_R0${R0}.txt

# Compile code
echo "COMPILING C++ CODE"
echo
bash $CODE_DIR/clean
make Periodic="-DPERIODIC" Bispectrum="-DBISPECTRUM" --directory $CODE_DIR

# Compute bispectrum counts
echo
echo "COMPUTING BISPECTRUM COUNTS"

$CODE_DIR/power -in $DATA -binfile $BINFILE -output $CODE_DIR/output -out_string ${STRING} -max_l $MAX_L -R0 $R0 -nthread $NTHREADS -perbox -f_rand $FRAND

# Ensure that the files have actually been created
if ! (test -f "$OUTPUT_FILE"); then
    echo
    echo "Bispectrum has not been computed. This indicates an error. Exiting."
    echo
    exit 1;
fi

echo
echo "COMPUTATIONS COMPLETE AND SAVED TO $OUTPUT_FILE"
echo
