#!/bin/bash
# Bash wrapper for the HIPSTER code developed by Oliver Philcox & Daniel Eisenstein

# Define Inputs
SHORT=h
LONG=dat:,l_max:,R0:,k_bin:,nthreads:,string:

# read the options
OPTS=$(getopt --options $SHORT --long $LONG --name "$0" -- "$@")

if [ $? != 0 ] ; then echo "Failed to parse options...exiting." >&2 ;
exit 1 ; fi

eval set -- "$OPTS"

 # set initial values
PARAM_COUNT=0
NTHREADS=10

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
    echo "--string: (Optional): Identification string for output file names."
    echo "--nthreads: (Optional): Number of CPU threads on which to run. Default: 10"
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
    --string )
      STRING="$2"
      shift 2
      ;;
    --nthreads)
      NTHREADS="$2"
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
if [ "$PARAM_COUNT" -ne 4 ]; then echo 'Not all command line parameters specified! Exiting;'; exit 1; fi


# Print the variables
echo
echo "INPUT PARAMETERS"
echo "----------------"
echo "Data file: $DATA"
echo "Maximum Legendre multipole: $MAX_L"
echo "Pair count truncation radius: $R0"
echo "k-space binning file: $BINFILE"
echo "Output string: $STRING"
echo "CPU-threads: $NTHREADS"
echo

# Check some variables:

if [ "$MAX_L" -ge 5 ]; then echo "Only multipoles up to ell = 4 currently implemented. Exiting;"; exit 1; fi;
if [ $((MAX_L%2)) -eq 1 ]; then echo "Maximum Legendre multipole must be even. Exiting;"; exit 1; fi;

if ! ( test -f "$DATA" ); then
    echo "Data file: $DATA does not exist. Exiting;"; exit 1;
fi
if ! ( test -f "$BINFILE" ); then
    echo "Binning file: $BINFILE does not exist. Exiting;"; exit 1;
fi

# Work out where the code is installed
CODE_DIR=`dirname "$0"`

# Compute number of k bins in file
K_BINS=`wc -l < $BINFILE`
let K_BINS++

# Define file names
CORRECTION_FILE=$CODE_DIR/output/correction_function_${STRING}_R0${R0}_lmax${MAX_L}.txt
DD_FILE=$CODE_DIR/output/${STRING}_DD_power_counts_n${K_BINS}_l${MAX_L}_full.txt
OUTPUT_FILE=$CODE_DIR/output/${STRING}_power_spectrum_n{$K_BINS}_l${MAX_L}.txt

# If correction file does not exist re-create it!
if $PRELOAD; then
    if ! (test -f "$CORRECTION_FILE"); then
        echo "Survey correction file: $CORRECTION_FILE has not been previously computed. The correction function and RR counts will now be recomputed from scratch."
        PRELOAD=false
    fi
fi

# Compile code
echo
echo "COMPILING C++ CODE"
echo
bash $CODE_DIR/clean
make --directory $CODE_DIR Periodic="-DPERIODIC"

# Check that the preloaded RR counts actually exist!
if $PRELOAD; then
    if ! (test -f "$RR_FILE"); then
        echo "Weighted RR count file $RR_FILE has not been previously computed. This will be computed from scratch."
        PRELOAD=false
    fi
fi

# Compute DD pair counts (always need to be computed)
echo
echo "COMPUTING DD PAIR COUNTS"
$CODE_DIR/power -in $DATA -in2 $DATA -binfile $BINFILE -output $CODE_DIR/output/ -out_string ${STRING}_DD -max_l $MAX_L -R0 $R0 -nthread $NTHREADS $PERIODIC_FLAG
# Ensure that the files have actually been created
if ! (test -f "$DD_FILE"); then
    echo
    echo "Weighted DD counts have not been computed. This indicates an error. Exiting."
    exit 1;
fi

# Combine files and output
echo
echo "COMBINING PAIR COUNTS TO FORM POWER SPECTRUM"
echo
python $CODE_DIR/python/reconstruct_power_periodic.py $DD_FILE $DATA $OUTPUT_FILE

echo "COMPUTATIONS COMPLETE"
