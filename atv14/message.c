#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define NUM_EXCHANGES 5
#define MESSAGE_CONTENT 42

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    printf("Numero do processo: %d\n", world_rank);
    printf("Tamanho da lista: %d\n", world_size);

    if (world_size < 2) {
        fprintf(stderr, "Este programa deve ser executado com pelo menos 2 processos.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int message = 55;
    double start_time, end_time, total_time;

    // Sincroniza os processos antes de iniciar o temporizador
    MPI_Barrier(MPI_COMM_WORLD);

    // Inicia o temporizador no processo 0
    if (world_rank == 0) {
        start_time = MPI_Wtime();
    }

    for (int i = 0; i < NUM_EXCHANGES; ++i) {
        if (world_rank == 0) {
            // O processo 0 envia a mensagem para o processo 1
            printf("Mensagem de tamanho %ld enviado pelo processo %d para o processo %d",sizeof(message),world_rank,1);
            MPI_Send(&message, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);

            // O processo 0 recebe a mensagem de volta do processo 1
            MPI_Recv(&message, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Mensagem de tamanho %ld recebido pelo processo %d do processo %d",sizeof(message),world_rank,1);

        } else if (world_rank == 1) {
            // O processo 1 recebe a mensagem do processo 0
            MPI_Recv(&message, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // O processo 1 envia a mesma mensagem de volta para o processo 0
            MPI_Send(&message, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }

    // Sincroniza os processos antes de parar o temporizador
    MPI_Barrier(MPI_COMM_WORLD);

    // Para o temporizador e calcula o tempo total no processo 0
    if (world_rank == 0) {
        end_time = MPI_Wtime();
        total_time = end_time - start_time;
        printf("Tempo total para %d trocas de mensagens: %f segundos\n", NUM_EXCHANGES, total_time);
        printf("Tempo médio por troca: %f microssegundos\n", (total_time / NUM_EXCHANGES) * 1e6);
    }

    // Finaliza o ambiente MPI
    MPI_Finalize();

    return 0;
}
//mpicc message.c -o message
//mpirun ./message --allow-run-as-root