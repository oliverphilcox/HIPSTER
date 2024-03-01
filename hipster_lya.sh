#!/bin/bash
# Bash wrapper for the HIPSTER-Lya code developed by Oliver Philcox & Roger de Belsunce

# Define Inputs
SHORT=h
LONG=dat:,l_max:,R0:,k_bin:,nthreads:,string:

# read the options
OPTS=$(getopt --options $SHORT --long $LONG --name "$0" -- "$@")

if [ $? != 0 ] ; then echo "Failed to parse options...exiting." >&2 ;
exit 1 ; fi

eval set -- "$OPTS"

 # set initial values
PERIODIC=false
PARAM_COUNT=0
NTHREADS=10
STRING=hipster

# Help dialogue
function usageText ()
{
    echo "USAGE"
    echo "-----"
    echo "-h: Display this message"
    echo "--dat: Data file in (x,y,z,weight*density,skewer_id) co-ordinates"
    echo "--l_max: Maximum Legendre multipole"
    echo "--R0: Pair count truncation radius"
    echo "--k_bin: k-space binning file"
    echo "--string (Optional): Identification string for output file names."
    echo "--nthreads (Optional): Number of CPU threads on which to run. Default: 10"
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
if [ "$PARAM_COUNT" -ne 4 ]; then echo; echo "Not all command line parameters specified!"; echo; usageText; exit 1; fi

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
echo "Periodic: $PERIODIC"
echo

# Check some variables:

if [ "$MAX_L" -ge 7 ]; then echo "Only multipoles up to ell = 6 currently implemented. Exiting;"; exit 1; fi;
if [ $((MAX_L%2)) -eq 1 ]; then echo "Maximum Legendre multipole must be even. Exiting;"; exit 1; fi;

if $PERIODIC; then echo 'Assuming periodic boundary conditions';
else echo 'Assuming non-periodic boundary conditions'; fi;

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
#let K_BINS++

# Define file names
R0int=$( printf "%.0f" $R0 )
POWER_FILE=$CODE_DIR/output/${STRING}_power_spectrum_n${K_BINS}_l${MAX_L}_R0${R0int}.txt

# Compile code
echo
echo "COMPILING C++ CODE"
pushd $CODE_DIR
echo
bash clean
popd
make Lya="-DLYA" --directory $CODE_DIR

# Compute power spectrum
echo
echo "RUNNING HIPSTER"
$CODE_DIR/power -in $DATA -in2 $DATA -binfile $BINFILE -output $CODE_DIR/output -out_string ${STRING} -max_l $MAX_L -R0 $R0 -nthread $NTHREADS $PERIODIC_FLAG
# Ensure that the files have actually been created
if ! (test -f "$POWER_FILE"); then
    echo
    echo "Power spectrum has not been computed. This indicates an error. Exiting."
    exit 1;
fi

echo
echo "COMPUTATIONS COMPLETE"
echo
