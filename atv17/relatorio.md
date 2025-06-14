# Relatório Atividade 17

Aluno: Davi Diógenes Ferreira de Almeida
Disciplina: IMD1116 Computação de Alto Desempenho/2025.1

Para esta atividade, implementamos um programa em MPI que multiplica uma matriz NxN por um vetor de tamanho N. As colunas da matriz são divididas entre proces$

O código desenvolvido foi o seguinte:

- Inicialmente, fazemos o broadcast das dimensões da matriz para todos os processos
```c
MPI_Bcast(&M_dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&N_dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
```
- Após isso, alocamos a matriz, o vetor que será multiplicado e o vetor resultante, além de dividir o número de colunas da matriz e alocar o espaço da fração $
```c
rows_per_proc = N_dim / size;
if (rank == 0) {
    A_glob = (double*)malloc((long)M_dim * N_dim * sizeof(double));
    x_glob = (double*)malloc((long)N_dim * sizeof(double));
    y_glob = (double*)malloc((long)M_dim * sizeof(double));
    // Criação de elementos
}
local_A = (double*)malloc((long)rows_per_proc * N_dim * sizeof(double));
local_y = (double*)malloc((long)rows_per_proc * sizeof(double));
```

- Em seguida, usamos MPI_Type_vector para descrever um padrão para armazenar

```c
MPI_Datatype column_type, resized_column_type;
//N de blocos(N de linhas)|Blocklength(colunas/process)|Stride(Tamanho da linha)|Novo Tipo
MPI_Type_vector(M_dim, cols_per_proc, N_dim, MPI_DOUBLE, &column_type);
//Extende o datatype (0, cols_per_proc * sizeof(double))
//É tipo um stride para o scatter
MPI_Type_create_resized(column_type, 0, cols_per_proc * sizeof(double), &resized_column_type);
MPI_Type_commit(&resized_column_type);
```

Após executar com diferentes tamanhos de matriz e quantidades de processo, os resultados foram:
- Tempo Total (ms)
