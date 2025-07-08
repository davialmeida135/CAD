#!/bin/bash
#SBATCH --partition gpu-4-a100
#SBATCH --gpus-per-node=1   # Número GPUs por nó
#SBATCH --nodes 1
#SBATCH --time 00:02:00
#SBATCH --job-name vadd
#SBATCH --output saida-%j.out

cd $SLURM_SUBMIT_DIR

ulimit -s unlimited

./heat