#!/bin/bash
# Bash wrapper for the HIPSTER-Lya code developed by Oliver Philcox & Roger de Belsunce

########## computation env 
# Load any necessary modules or activate virtual environment

########## home dir
HOMEDIR=</enter/home/directory/>
cd $HOMEDIR

########## input parameters for code 
LMAX=4 #multipoles for power spectrum
RMAX=200 # pair truncation parameter
LMAX_RR=8 # max multipole for window matrix
NBIN_RR=1000 # number of bins for window matrix
NTHREADS=256 # number of parallel threads for computation 
K_DIR=${HOMEDIR}'/binning.csv'

########## data 
DIR= </enter/data/directory/>
DATA=${DIR}'lya_data.xyzwrwgth'
STRING='pk_hipster'

########## compile code

# Compile code once before running it
# Work out where the code is installed
echo
echo "COMPILING C++ CODE"
pushd $HOMEDIR
echo
bash clean
popd
make Lya="-DLYA" --directory $HOMEDIR

########## define .EXE
EXE="./hipster_wrapper_lya.sh \
        --dat ${DATA} \
        --l_max ${LMAX} \
        --R0 ${RMAX} \
        --k_bin ${K_DIR} \
        --nthread ${NTHREADS} \
        --string ${STRING} \
        --l_max_RR ${LMAX_RR} \
        --nbin_RR ${NBIN_RR} "

# print command
echo ${EXE}

echo " =============== Job started at $(date)  =============== "
# run command to submit job to queue
$EXE
echo " =============== Job ended at $(date)  =============== "