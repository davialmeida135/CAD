# Relatório Atividade 16

Aluno: Davi Diógenes Ferreira de Almeida
Disciplina: IMD1116 Computação de Alto Desempenho/2025.1

Para esta atividade, implementamos um programa em MPI que multiplica uma matriz NxN por um vetor de tamanho N. As linhas da matriz são divididas entre processos por meio de MPI_Scatter e depois agregadas em um vetor com MPI_Gather, já o vetor é passado por completo para todos os processos com MPI_Bcast.

O código desenvolvido foi o seguinte:

- Inicialmente, fazemos o broadcast das dimensões da matriz para todos os processos
```c
MPI_Bcast(&M_dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&N_dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
``` 
- Após isso, alocamos a matriz, o vetor que será multiplicado e o vetor resultante, além de dividir o número de linhas da matriz e alocar o espaço da fração da matriz que será atribuida para cada processo
```c
rows_per_proc = M_dim / size;
if (rank == 0) {
    A_glob = (double*)malloc((long)M_dim * N_dim * sizeof(double));
    x_glob = (double*)malloc((long)N_dim * sizeof(double));
    y_glob = (double*)malloc((long)M_dim * sizeof(double));
    // Criação de elementos        
}
local_A = (double*)malloc((long)rows_per_proc * N_dim * sizeof(double));
local_y = (double*)malloc((long)rows_per_proc * sizeof(double));
```
- Em seguida, é feito o broadcast do vetor para todos os processos e o espalhamento de partes de mesmo tamanho da matriz.
```c
// Broadcast vetor B
MPI_Bcast(x_glob, N_dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
// Scatter matriz A
// |A global|Tamanho da parte|tipo|A local| Tamanho da parte| Tipo|Origem|Cmm world|
MPI_Scatter(A_glob, rows_per_proc * N_dim, MPI_DOUBLE, local_A, rows_per_proc * N_dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
```
- Por fim, a computação é feita e o vetor resultante é agregado
```c
// Local computation: local_y = local_A * x_glob
for (int i = 0; i < rows_per_proc; i++) {
    local_y[i] = 0.0;
    for (int j = 0; j < N_dim; j++) {
        local_y[i] += local_A[i * N_dim + j] * x_glob[j];
    }
}

MPI_Gather(local_y, rows_per_proc, MPI_DOUBLE, y_glob, rows_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

```

Após executar com diferentes tamanhos de matriz e quantidades de processo, os resultados foram:
- Tempo Total (ms)

| Matrix Dim | 2 Procs | 4 Procs | 5 Procs | 8 Procs | 10 Procs |
|------------|---------|---------|---------|---------|----------|
| 1000       | 7.390   | 8.553   | 10.485  | 9.490   | 9.070    |
| 2000       | 20.093  | 23.914  | 28.544  | 25.070  | 20.465   |
| 4000       | 61.432  | 52.839  | 68.223  | 52.367  | 56.026   |
| 8000       | 219.335 | 177.636 | 245.484 | 141.496 | 172.138  |
| 10000      | 332.762 | 252.090 | 296.066 | 229.939 | 313.371  |
