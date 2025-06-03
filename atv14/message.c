#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 

// Exemplo: 1 KB = 1024
//#define MESSAGE_SIZE 1024*1024*24

int main(int argc, char *argv[]) {
    int rank, size;
    long MESSAGE_SIZE;
    char *message_sent = NULL;     // Ponteiro para o buffer de envio
    char *message_received = NULL; // Ponteiro para o buffer de recebimento
    MPI_Status status;
    double start_time, end_time, total_time;

    // Define numero de trocas como 1000 se não for declarado
    int num_exchanges = 1000;
        // Verifica se o tamanho da mensagem foi passado como argumento
    if (argc > 1) {
        MESSAGE_SIZE = atol(argv[1]); // Usa atol para números maiores
        if (MESSAGE_SIZE <= 0) {
            fprintf(stderr, "Tamanho da mensagem inválido. Usando o padrão: 1MB\n");
            MESSAGE_SIZE = 1024 * 1024; // Define um padrão, por exemplo 1MB
        }
    } else {
        // Se nenhum tamanho de mensagem for fornecido, usa um padrão ou sai com erro
        fprintf(stderr, "Uso: %s [num_exchanges] <message_size_bytes>\n", argv[0]);
        fprintf(stderr, "Definindo tamanho da mensagem padrão para 1MB.\n");
        MESSAGE_SIZE = 1024 * 1024; // Define um padrão, por exemplo 1MB
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

    // Alocação dinâmica dos buffers de mensagem
    message_sent = (char *)malloc(MESSAGE_SIZE);
    message_received = (char *)malloc(MESSAGE_SIZE);

    if (message_sent == NULL || message_received == NULL) {
        fprintf(stderr, "Rank %d: Falha ao alocar buffers de %d bytes.\n", rank, MESSAGE_SIZE);
        perror("Detalhe do erro de malloc");
        MPI_Abort(MPI_COMM_WORLD, 1); // Aborta todos os processos MPI
        // Não é necessário 'return 1;' aqui pois MPI_Abort já termina.
    }

    // O comprimento da mensagem é o MESSAGE_SIZE completo
    int message_length = MESSAGE_SIZE;

    // Fill mensagem com letras (APÓS alocação bem-sucedida)
    if (rank == 0) { // Apenas o processo 0 precisa preencher o buffer que ele envia
        for (int i = 0; i < MESSAGE_SIZE; i++) {
            message_sent[i] = 'A' + (i % 26); // Preenche com A-Z repetidamente
        }
    }


    MPI_Barrier(MPI_COMM_WORLD); // Sincroniza os processos antes de iniciar o timer

    if (rank == 0) {
        printf("Iniciando %d trocas de mensagens (ping-pong).\n", num_exchanges);
        printf("Tamanho da mensagem: %d bytes\n", message_length);
        start_time = MPI_Wtime();

        for (int i = 0; i < num_exchanges; i++) {
            // 1. Processo 0 envia para o processo 1
            // |mensagem|tamanho|tipo|pra quem|id|comm_world|
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
            double effective_bandwidth_gbps = (total_data_moved_one_way_gb * 2) / total_time; // GB/s, *2 para ida e volta
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