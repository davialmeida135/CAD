# Computação de Alto Desempenho
Repositório das atividades da disciplina de IMD1116 - Computação de Alto Desempenho 2025.1

## Compilação com openmp
gcc -fopenmp -o programa programa.c -lm

## Compilação com MPI
mpicc -o programa programa.c -lm

## Compilação para programação em GPU no NPAD
```bash
#GPU
module load compilers/nvidia/nvhpc/24.11

sinfo | grep gpu

nvc gputest.c -mp=gpu

squeue -u $USER
```
