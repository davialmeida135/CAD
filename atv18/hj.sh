#!/bin/bash
#SBATCH --partition gpu-8-v100 
#SBATCH --nodes 1
#SBATCH --time 00:02:00
#SBATCH --job-name vadd
#SBATCH --output saida-%j.out

cd $SLURM_SUBMIT_DIR

ulimit -s unlimited

./heat
