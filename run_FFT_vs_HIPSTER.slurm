#!/bin/bash
#
#SBATCH --job-name=quijote_hipster
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ophilcox@princeton.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=04:00:01
#SBATCH --mem-per-cpu=1000
#SBATCH --output=/home/ophilcox/out/hipster_quijote_%a.log

# Load modules
module load rh/devtoolset/8

sim=$SLURM_ARRAY_TASK_ID

bash run_hipster
