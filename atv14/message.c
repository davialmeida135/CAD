#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 

#define MESSAGE_SIZE 1024 * 1024// Exemplo: 1 KB = 1024


int main(int argc, char *argv[]) {
    int rank, size;

    char message_sent[MESSAGE_SIZE];
    char message_received[MESSAGE_SIZE];
    MPI_Status status;
    double start_time, end_time, total_time;

    // Define numero de trocas como 1000 se não for declarado
    int num_exchanges = 1000;
    if (argc > 1) {
        num_exchanges = atoi(argv[1]);
        if (num_exchanges <= 0) {
            num_exchanges = 1000;
        }
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != 2) {
        if (rank == 0) {
            fprintf(stderr, "Erro: Este programa deve ser executado com exatamente 2 processos.\n");
            fprintf(stderr, "Uso: mpirun -np 2 %s [numero_de_trocas]\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    // Fill mensagem com letras
    for (int i = 0; i < MESSAGE_SIZE; i++) {
        message_sent[i] = 'A' + (i % 26); // Preenche com A-Z repetidamente
    }

    int message_length = MESSAGE_SIZE;
    MPI_Barrier(MPI_COMM_WORLD); // Sincroniza os processos antes de iniciar o timer

    if (rank == 0) {
        printf("Iniciando %d trocas de mensagens (ping-pong).\n", num_exchanges);
        printf("Tamanho da mensagem: %d bytes\n", message_length);
        start_time = MPI_Wtime();

        for (int i = 0; i < num_exchanges; i++) {
            // 1. Processo 0 envia para o processo 1
            MPI_Send(message_sent, message_length, MPI_CHAR, 1, 0, MPI_COMM_WORLD);

            // 4. Processo 0 recebe do processo 1
            MPI_Recv(message_received, message_length, MPI_CHAR, 1, 0, MPI_COMM_WORLD, &status);
        }

        end_time = MPI_Wtime();
        total_time = end_time - start_time;

        printf("--------------------------------------------------\n");
        printf("Tempo total para %d trocas: %f segundos\n", num_exchanges, total_time);
        printf("Tempo médio por troca (RTT): %e segundos\n", total_time / num_exchanges);
        if (message_length > 0 && total_time > 0) {
            double total_data_moved_one_way_gb = (double)num_exchanges * message_length / (1024.0*1024.0*1024.0);
            double effective_bandwidth_gbps = (total_data_moved_one_way_gb * 2) / total_time; // GB/s, *2 para ping-pong
            printf("Largura de banda efetiva (ping-pong): %.3f GB/s\n", effective_bandwidth_gbps * 8); // Gbps
        }
        printf("--------------------------------------------------\n");

    } else if (rank == 1) {
        for (int i = 0; i < num_exchanges; i++) {
            // 2. Processo 1 recebe do processo 0
            MPI_Recv(message_received, message_length, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);

            // 3. Processo 1 envia a mesma mensagem de volta para o processo 0
            MPI_Send(message_received, message_length, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();
    return 0;
}