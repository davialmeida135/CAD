# Relatório Atividade 17

Aluno: Davi Diógenes Ferreira de Almeida
Disciplina: IMD1116 Computação de Alto Desempenho/2025.1

Para esta atividade, implementamos um programa em MPI que multiplica uma matriz NxN por um vetor de tamanho N. As colunas da matriz são divididas entre processos.

O código desenvolvido foi o seguinte:

- Inicialmente, fazemos o broadcast das dimensões da matriz para todos os processos
```c
MPI_Bcast(&M_dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&N_dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
```
- Após isso, alocamos a matriz, o vetor que será multiplicado e o vetor resultante, além de dividir o número de colunas da matriz e alocar o espaço da fração
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

- Em seguida, usamos MPI_Type_vector para descrever um padrão para armazenar as colunas

```c
MPI_Datatype column_type, resized_column_type;
//N de blocos(N de linhas)|Blocklength(colunas/process)|Stride(Tamanho da linha)|Novo Tipo
MPI_Type_vector(M_dim, cols_per_proc, N_dim, MPI_DOUBLE, &column_type);

//Reduz o tamanho do bloco para o tamanho que é capturado por cada bloco (0, cols_per_proc * sizeof(double))
//Segmentation fault se não usar isso
//Originalmente o 2o objeto começa apenas após o fim do ultimo bloco
MPI_Type_create_resized(column_type, 0, cols_per_proc * sizeof(double), &resized_column_type);
MPI_Type_commit(&resized_column_type);
```

- Espalhamos as colunas e os pedaços do vetor x para os processos

```c
MPI_Scatter(A_glob, 1, resized_column_type, local_A, M_dim * cols_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Scatter(x_glob, cols_per_proc, MPI_DOUBLE, local_x, cols_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
```

- Executamos os cálculos parciais
```c
// Cada processo calcula y parcial
for (int i = 0; i < M_dim; i++) {
    for (int j = 0; j < cols_per_proc; j++) {
        local_y[i] += local_A[i * cols_per_proc + j] * local_x[j];
    }
}
```

- Somamos todas as contas parciais por meio de MPI_Reduce

```c
//Vetor local| Vetor soma| Tamanho| Tipo| Operação| Destino| Comm_world
MPI_Reduce(local_y, y_glob, M_dim, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
```

## Resultados

Testing on 2 MPI processes (column-wise distribution):

| N     | Col major time | Row major time |
|-------|---------------|----------------|
| 1000  | 0.006251      | 0.007390       |
| 2000  | 0.016031      | 0.020093       |
| 4000  | 0.058866      | 0.061432       |
| 8000  | 0.224108      | 0.219335       |
| 10000 | 0.347464      | 0.332762       |

Testing on 4 MPI processes (column-wise distribution):

| N     | Col major time | Row major time |
|-------|---------------|----------------|
| 1000  | 0.012371      | 0.008553       |
| 2000  | 0.020498      | 0.023914       |
| 4000  | 0.063258      | 0.052839       |
| 8000  | 0.188353      | 0.177636       |
| 10000 | 0.281764      | 0.252090       |

Os resultados mostram que os tempos foram muito parecidos e, com tamanhos de matriz maiores, o formato Row-major trouxe tempos menores.
O formato de multiplicação por colunas aumenta muito o número