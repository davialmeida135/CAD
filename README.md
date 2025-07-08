Repositório das atividades da disciplina de Computação de Alto Desempenho 2025.1

gcc -fopenmp -o kkkk kkkk.c


```bash
#GPU
module load compilers/nvidia/nvhpc/24.11

sinfo | grep gpu

nvc gputest.c -mp=gpu

squeue -u $USER
```