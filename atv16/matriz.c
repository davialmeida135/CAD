#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

int main(int argc, char *argv[]) {
    int rank, size;
    int M, N;
    double *A = NULL, *x = NULL, *y = NULL;
    double *local_A, *local_y;
    int rows_per_proc;
    double start_time, end_time, local_start_time, local_end_time;
    double computation_time, communication_time, total_time;
    double local_computation_start, local_computation_end, local_communication_start, local_communication_end;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        // Usando valores fixos para M e N para simplicidade
        M = 10000;
        N = 10000;
        printf("Executando com Matriz A(%d x %d) e Vetor x(%d)\n", M, N);

        if (M % size != 0) {
            fprintf(stderr, "Erro: O número de linhas (M) deve ser divisível pelo número de processos.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        A = (double*)malloc(M * N * sizeof(double));
        x = (double*)malloc(N * sizeof(double));
        y = (double*)malloc(M * sizeof(double));

        // Inicializando matriz A e vetor x com valores aleatórios
        srand(time(NULL));
        for (int i = 0; i < M * N; i++) A[i] = (double)rand() / RAND_MAX;
        for (int i = 0; i < N; i++) x[i] = (double)rand() / RAND_MAX;
    }

    // Barreira para garantir que todos os processos iniciem a medição de tempo juntos
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    // Envia M e N para todos os processos
    local_communication_start = MPI_Wtime();
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    local_communication_end = MPI_Wtime();
    communication_time = local_communication_end - local_communication_start;


    rows_per_proc = M / size;
    local_A = (double*)malloc(rows_per_proc * N * sizeof(double));

    if(rank != 0) {
        x = (double*)malloc(N * sizeof(double));
    }

    // Distribui o vetor x para todos os processos
    local_communication_start = MPI_Wtime();
    MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    local_communication_end = MPI_Wtime();
    communication_time += local_communication_end - local_communication_start;

    // Distribui as linhas da matriz A entre os processos
    local_communication_start = MPI_Wtime();
    MPI_Scatter(A, rows_per_proc * N, MPI_DOUBLE, local_A, rows_per_proc * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    local_communication_end = MPI_Wtime();
    communication_time += local_communication_end - local_communication_start;


    local_y = (double*)malloc(rows_per_proc * sizeof(double));

    // Cálculo local do produto matriz-vetor
    local_computation_start = MPI_Wtime();
    for (int i = 0; i < rows_per_proc; i++) {
        local_y[i] = 0.0;
        for (int j = 0; j < N; j++) {
            local_y[i] += local_A[i * N + j] * x[j];
        }
    }
    local_computation_end = MPI_Wtime();
    computation_time = local_computation_end - local_computation_start;

    // Reúne os resultados parciais de y no processo 0
    local_communication_start = MPI_Wtime();
    MPI_Gather(local_y, rows_per_proc, MPI_DOUBLE, y, rows_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    local_communication_end = MPI_Wtime();
    communication_time += local_communication_end - local_communication_start;

    end_time = MPI_Wtime();
    total_time = end_time - start_time;

    // Imprime os tempos de cada processo
    printf("Processo %d: Tempo de Computação = %f s, Tempo de Comunicação = %f s, Tempo Total = %f s\n",
           rank, computation_time, communication_time, total_time);


    // O processo 0 pode verificar os resultados (opcional)
    if (rank == 0) {
        printf("\nResultados agregados no processo 0.\n");
        printf("Tempo total de execução (Processo 0): %f segundos\n", total_time);
        // Exemplo para verificar o primeiro elemento de y
        // double y0_check = 0.0;
        // for(int j=0; j<N; j++) y0_check += A[j] * x[j];
        // printf("y[0] (calculado) = %f, y[0] (verificação) = %f\n", y[0], y0_check);

        free(A);
        free(y);
    }

    free(x);
    free(local_A);
    free(local_y);

    MPI_Finalize();
    return 0;
}