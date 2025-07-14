# Computação de Alto Desempenho
Repositório das atividades da disciplina de IMD1116 - Computação de Alto Desempenho 2025.1

[Primeiros passos NPAD](https://npad.ufrn.br/npad/primeirospassos)

## Compilação com openmp
gcc -fopenmp -o programa programa.c -lm

## Compilação com MPI
mpicc -o programa programa.c -lm

## Compilação para programação em GPU no NPAD
[Repositório tutorial de GPU do NPAD](https://github.com/NPAD-UFRN/openmp-tutorial)
```bash
#GPU
module load compilers/nvidia/nvhpc/24.11

sinfo | grep gpu

nvc gputest.c -mp=gpu

squeue -u $USER
```
