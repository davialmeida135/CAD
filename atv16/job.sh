#!/bin/bash
#SBATCH --job-name=MPI_hello  
#SBATCH --output=saida%j.out
#SBATCH --error=erro%j.err
#SBATCH --partition=amd-512
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-0:5

mpirun -n 20 ./matriz