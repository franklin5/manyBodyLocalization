#!/bin/bash
#SBATCH --partition=commons
#SBATCH --nodes=1
##SBATCH --exclusive
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=30000m
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=shangshun89@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --export=ALL

echo "My job ran on:"
echo $SLURM_NODELIST
cd $SLURM_SUBMIT_DIR
matlab -r "fock; decimation; quit" > log.txt
