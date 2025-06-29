#!/bin/bash
#SBATCH --partition=gpu-4-a100
#SBATCH --gpus-per-node=1  # Número GPUs por nó
#SBATCH --cpus-per-task=6
#SBATCH --time=0-03:00
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

./vadd
