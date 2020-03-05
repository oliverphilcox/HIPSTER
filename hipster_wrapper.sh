#!/bin/bash
# Bash wrapper for the HIPSTER code developed by Oliver Philcox & Daniel Eisenstein

# Define Inputs
SHORT=h
LONG=dat:,ran_DR:,ran_RR:,l_max:,R0:,k_bin:,nthreads:,string:,subsample:,load_RR

# read the options
OPTS=$(getopt --options $SHORT --long $LONG --name "$0" -- "$@")

if [ $? != 0 ] ; then echo "Failed to parse options...exiting." >&2 ;
exit 1 ; fi

eval set -- "$OPTS"

 # set initial values
PRELOAD=false
PERIODIC=false
PERIODIC_MAKEFLAG=
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
    echo "--dat: Data file in (x,y,z,weight) co-ordinates"
    echo "--ran_DR: Random file for DR pair counting"
    echo "--ran_RR: Random file for RR pair counting (and survey correction function estimation)"
    echo "--l_max: Maximum Legendre multipole"
    echo "--R0: Pair count truncation radius"
    echo "--k_bin: k-space binning file"
    echo "--string (Optional): Identification string for output file names."
    echo "--nthreads (Optional): Number of CPU threads on which to run. Default: 10"
    echo "--load_RR: If set, load previously computed RR pair counts and survey correction functions for a large speed boost."
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
    --ran_DR )
      RAN_DR="$2"
      shift 2
      let PARAM_COUNT++
      ;;
    --ran_RR )
      RAN_RR="$2"
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
    --subsample )
      SUBSAMPLE="$2"
      shift 2
      ;;
    --load_RR )
      PRELOAD=true
      shift
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
if [ "$PARAM_COUNT" -ne 6 ]; then echo; echo "Not all command line parameters specified!"; echo; usageText; exit 1; fi


# Print the variables
echo
echo "INPUT PARAMETERS"
echo "----------------"
echo "Data file: $DATA"
echo "Random file for RR counts: $RAN_RR"
echo "Random file for DR counts: $RAN_DR"
echo "Maximum Legendre multipole: $MAX_L"
echo "Pair count truncation radius: $R0"
echo "k-space binning file: $BINFILE"
echo "Output string: $STRING"
echo "CPU-threads: $NTHREADS"
echo "Periodic: $PERIODIC"
echo "Subsampling: $SUBSAMPLE"
echo

# Check some variables:

if [ "$MAX_L" -ge 7 ]; then echo "Only multipoles up to ell = 6 currently implemented. Exiting;"; exit 1; fi;
if [ $((MAX_L%2)) -eq 1 ]; then echo "Maximum Legendre multipole must be even. Exiting;"; exit 1; fi;

if $PRELOAD; then echo 'Using pre-computed survey correction function and (weighted) random pair counts';
else echo 'Survey correction function and (weighted) random pair counts will be computed from scratch'; fi

if $PERIODIC; then echo 'Assuming periodic boundary conditions';
else echo 'Assuming non-periodic boundary conditions'; fi;

if ! ( test -f "$DATA" ); then
    echo "Data file: $DATA does not exist. Exiting;"; exit 1;
fi
if ! ( test -f "$RAN_RR" ); then
    echo "RR random file: $RAN_RR does not exist. Exiting;"; exit 1;
fi
if ! ( test -f "$RAN_DR" ); then
    echo "RR random file: $RAN_DR does not exist. Exiting;"; exit 1;
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

# Subsample files if necessary
if [[ $(bc <<< "$SUBSAMPLE > 1.") -eq 1 ]]; then
  # Count number of particles in file
  N_GAL=`wc -l < $DATA`
  let N_GAL++
  if ! $PRELOAD;  then
    LEN_RAND_RR=`wc -l < $RAN_RR`
    let LEN_RAND_RR++
  fi
  LEN_RAND_DR=`wc -l < $RAN_DR`
  let LEN_RAND_DR++

  # Compute number after subsampling and ensure it's an integer
  N_SUB=`echo "$N_GAL/$SUBSAMPLE" | bc`
  N_SUB=`echo "($N_SUB+0.5)/1" | bc`
  if ! $PRELOAD;  then
    N_SUB_RR=`echo "$LEN_RAND_RR/$SUBSAMPLE" | bc`
    N_SUB_RR=`echo "($N_SUB_RR+0.5)/1" | bc`
  fi
  N_SUB_DR=`echo "$LEN_RAND_DR/$SUBSAMPLE" | bc`
  N_SUB_DR=`echo "($N_SUB_DR+0.5)/1" | bc`

  # Run subsampling script and create new file names
  NEW_DATA=$DATA.sub
  echo "Subsampling data with subsampling ratio $SUBSAMPLE"
  python $CODE_DIR/python/take_subset_of_particles.py $DATA $NEW_DATA $N_SUB
  DATA=$NEW_DATA

  if ! $PRELOAD;  then
    NEW_RR=$RAN_RR.sub
    echo "Subsampling RR file with subsampling ratio $SUBSAMPLE"
    python $CODE_DIR/python/take_subset_of_particles.py $RAN_RR $NEW_RR $N_SUB_RR
    RAN_RR=$NEW_RR
  fi

  NEW_DR=$RAN_DR.sub
  echo "Subsampling DR file with subsampling ratio $SUBSAMPLE"
  python $CODE_DIR/python/take_subset_of_particles.py $RAN_DR $NEW_DR $N_SUB_DR
  RAN_DR=$NEW_DR

  echo "Using $N_SUB particles in $DATA, from $N_GAL particles originally"
fi


# Define file names
R0int=$( printf "%.0f" $R0 )
CORRECTION_FILE=$CODE_DIR/output/correction_function_${STRING}_R0${R0int}_lmax${MAX_L}.txt
RR_FILE=$CODE_DIR/output/${STRING}_RR_power_counts_n${K_BINS}_l${MAX_L}_R0${R0int}.txt
DR_FILE=$CODE_DIR/output/${STRING}_DR_power_counts_n${K_BINS}_l${MAX_L}_R0${R0int}.txt
DD_FILE=$CODE_DIR/output/${STRING}_DD_power_counts_n${K_BINS}_l${MAX_L}_R0${R0int}.txt
OUTPUT_FILE=$CODE_DIR/output/${STRING}_power_spectrum_n${K_BINS}_l${MAX_L}_R0${R0int}.txt

# If correction file does not exist re-create it!
if $PRELOAD; then
    if ! (test -f "$CORRECTION_FILE"); then
        echo "Survey correction file: $CORRECTION_FILE has not been previously computed. The correction function and RR counts will now be recomputed from scratch."
        PRELOAD=false
    fi
fi

# Compute survey correction function if necessary
if ! $PRELOAD;  then
    echo
    echo "COMPUTING SURVEY CORRECTION FUNCTION"
    echo
    python $CODE_DIR/python/compute_correction_function.py $RAN_RR $CORRECTION_FILE $R0 $R0 100 $NTHREADS
fi

# Count number of randoms in each file
N_RAND_DR=`wc -l < $RAN_DR`
let N_RAND_DR++
N_RAND_RR=`wc -l < $RAN_RR`
let N_RAND_RR++

# Compile code
echo
echo "COMPILING C++ CODE"
echo
pushd $CODE_DIR
bash clean
popd
make --directory $CODE_DIR
# compile without periodic behavior

# Check that the preloaded RR counts actually exist!
if $PRELOAD; then
    if ! (test -f "$RR_FILE"); then
        echo "Weighted RR count file $RR_FILE has not been previously computed. This will be computed from scratch."
        PRELOAD=false
    fi
fi

# Compute RR counts if necessary
if ! $PRELOAD; then
      echo
      echo "COMPUTING RR PAIR COUNTS"
      $CODE_DIR/power -in $RAN_RR -in2 $RAN_RR -binfile $BINFILE -output $CODE_DIR/output -out_string ${STRING}_RR -max_l $MAX_L -R0 $R0 -inv_phi_file $CORRECTION_FILE -nthread $NTHREADS $PERIODIC_FLAG
      # Ensure that the files have actually been created
      if ! (test -f "$RR_FILE"); then
          echo
          echo "Weighted RR counts have not been computed. This indicates an error. Exiting."
          exit 1;
      fi
fi

# Compute DR pair counts (always need to be computed)
echo
echo "COMPUTING DR PAIR COUNTS"
$CODE_DIR/power -in $DATA -in2 $RAN_RR -binfile $BINFILE -output $CODE_DIR/output -out_string ${STRING}_DR -max_l $MAX_L -R0 $R0 -inv_phi_file $CORRECTION_FILE -nthread $NTHREADS $PERIODIC_FLAG
# Ensure that the files have actually been created
if ! (test -f "$DR_FILE"); then
    echo
    echo "Weighted DR counts have not been computed. This indicates an error. Exiting."
    exit 1;
fi


# Compute DD pair counts (always need to be computed)
echo
echo "COMPUTING DD PAIR COUNTS"
$CODE_DIR/power -in $DATA -in2 $DATA -binfile $BINFILE -output $CODE_DIR/output -out_string ${STRING}_DD -max_l $MAX_L -R0 $R0 -inv_phi_file $CORRECTION_FILE -nthread $NTHREADS $PERIODIC_FLAG
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
python $CODE_DIR/python/reconstruct_power.py $DD_FILE $DR_FILE $RR_FILE $DATA $N_RAND_RR $N_RAND_DR $OUTPUT_FILE
# Ensure this ran as expected
if ! (test -f "$OUTPUT_FILE"); then
    echo
    echo "Output power file has not been computed. This indicates an error. Exiting."
    exit 1;
fi

echo
echo "COMPUTATIONS COMPLETE"
echo
