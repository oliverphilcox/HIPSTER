#!/bin/bash
#
#SBATCH --job-name=quijote_pylians
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ophilcox@princeton.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=02:00:01
#SBATCH --mem-per-cpu=500
#SBATCH --array=0-100
#SBATCH --output=/home/ophilcox/out/pylians_%a

# Load modules
module load anaconda3
source activate py2

base=$SLURM_ARRAY_TASK_ID
let min=base*16
let max=base*16+16

for i in $(eval echo "{$min..$max}")
do
echo "Using iteration $i"
python pylians_bk.py $i
done

#bash run_hipster
#dat=/tigress/ophilcox/hipster_data/quijote_halos_1.txt
# ./power -in $dat -in2 $dat -in3 $dat -binfile /tigress/ophilcox/hipster_data/k_bin_small3.csv  -output /tigress/ophilcox/hipster_data -out_string sim1_run -perbox -max_l 4 -R0 50 -nthread 20
