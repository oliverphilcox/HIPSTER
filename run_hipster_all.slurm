#!/bin/bash
#
#SBATCH --job-name=quijote_hipster_all
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ophilcox@princeton.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=0-02:00:01
#SBATCH --mem-per-cpu=3000
#SBATCH --output=/home/ophilcox/out/hipster_run_quijote_%a
#SBATCH --array=0-100

sim=$SLURM_ARRAY_TASK_ID

binfile=/tigress/ophilcox/hipster_data/k_bin_big.csv
dat=/tigress/ophilcox/hipster_data/quijote_mass_cut_halo_$sim.txt
out_str=quijote.$sim

./power -in $dat -binfile $binfile  -output /tigress/ophilcox/hipster_data/quijote_all_bk -out_string $out_str -perbox -max_l 6 -R0 100 -nthread 20  -f_rand 3.

echo "Sim $sim complete"
