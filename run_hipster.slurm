#!/bin/bash
#
#SBATCH --job-name=quijote_hipster2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ophilcox@princeton.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=08:00:01
#SBATCH --mem-per-cpu=1000
#SBATCH --output=/home/ophilcox/out/hipster_run.slurm

# Load modules
module load rh/devtoolset/8

#sim=$SLURM_ARRAY_TASK_ID

bash run_hipster
#dat=/tigress/ophilcox/hipster_data/quijote_halos_1.txt
# ./power -in $dat -in2 $dat -in3 $dat -binfile /tigress/ophilcox/hipster_data/k_bin_small3.csv  -output /tigress/ophilcox/hipster_data -out_string sim1_run -perbox -max_l 4 -R0 50 -nthread 20
