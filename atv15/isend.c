#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h> // Para memcpy, se necessário (usaremos troca de ponteiros)

// Parâmetros da Simulação
#define L_DOMAIN 1.0                // Comprimento da barra
#define T_FINAL_TIME 0.01           // Tempo total de simulação
#define ALPHA_DIFFUSIVITY 0.01      // Difusividade térmica
#define N_GLOBAL_POINTS 33333       // Número total de pontos espaciais na barra global
#define CFL_STABILITY_FACTOR 0.4    // Fator CFL para estabilidade (<= 0.5)

// Condições de Contorno Globais
#define BC_GLOBAL_LEFT 0.0
#define BC_GLOBAL_RIGHT 0.0

// Tags MPI para comunicação não bloqueante
#define TAG_DATA_GOES_LEFT 0  // Mensagem contém ponto de borda que vai para o processo à esquerda
#define TAG_DATA_GOES_RIGHT 1 // Mensagem contém ponto de borda que vai para o processo à direita

// Função para calcular os parâmetros locais de cada processo (divisão de N_GLOBAL_POINTS)
void calculate_local_params(int rank, int mpi_size, int* local_n_ptr, int* global_offset_ptr) {
    int base_chunk = N_GLOBAL_POINTS / mpi_size;
    int remainder = N_GLOBAL_POINTS % mpi_size;
    if (rank < remainder) {
        *local_n_ptr = base_chunk + 1;
        *global_offset_ptr = rank * (base_chunk + 1);
    } else {
        *local_n_ptr = base_chunk;
        *global_offset_ptr = remainder * (base_chunk + 1) + (rank - remainder) * base_chunk;
    }
}

// Função para inicializar os dados (condição inicial u(x,0) = sin(pi*x/L))
void initialize_data(double* local_u_current, int local_n, int global_offset, double dx) {
    for (int i = 0; i < local_n; ++i) {
        double x = (global_offset + i) * dx;
        local_u_current[i + 1] = sin(M_PI * x / L_DOMAIN);
    }
    local_u_current[0] = 0.0; // Fantasma esquerdo inicializado
    if (local_n > 0) {
         local_u_current[local_n + 1] = 0.0; // Fantasma direito inicializado
    }
}

// Versão MODIFICADA: Simulação com MPI_Isend/MPI_Irecv e um loop de MPI_Wait
void run_simulation_nonblocking_individual_waits(double* u_curr, double* u_next,
                                                 int local_n, int rank, int mpi_size,
                                                 int num_steps, double factor) {
    MPI_Request requests[4]; // Máximo 2 envios e 2 recebimentos por processo
    int num_reqs;

    for (int t = 0; t < num_steps; ++t) {
        num_reqs = 0; // Reseta o contador de requisições para cada passo de tempo

        // 1. Iniciar todos os MPI_Irecv
        // Receber do vizinho da esquerda (se existir) -> preenche u_curr[0]
        if (rank > 0) {
            MPI_Irecv(&u_curr[0], 1, MPI_DOUBLE, rank - 1, TAG_DATA_GOES_RIGHT,
                      MPI_COMM_WORLD, &requests[num_reqs++]);
        }
        // Receber do vizinho da direita (se existir) -> preenche u_curr[local_n + 1]
        if (rank < mpi_size - 1) {
            MPI_Irecv(&u_curr[local_n + 1], 1, MPI_DOUBLE, rank + 1, TAG_DATA_GOES_LEFT,
                      MPI_COMM_WORLD, &requests[num_reqs++]);
        }

        // 2. Iniciar todos os MPI_Isend
        // Enviar para o vizinho da direita (se existir) o ponto u_curr[local_n]
        if (rank < mpi_size - 1) {
            MPI_Isend(&u_curr[local_n], 1, MPI_DOUBLE, rank + 1, TAG_DATA_GOES_RIGHT,
                      MPI_COMM_WORLD, &requests[num_reqs++]);
        }
        // Enviar para o vizinho da esquerda (se existir) o ponto u_curr[1]
        if (rank > 0) {
            MPI_Isend(&u_curr[1], 1, MPI_DOUBLE, rank - 1, TAG_DATA_GOES_LEFT,
                      MPI_COMM_WORLD, &requests[num_reqs++]);
        }

        // 3. Esperar que todas as comunicações não bloqueantes completem usando MPI_Wait individualmente
        //    (Em vez de MPI_Waitall)
        if (num_reqs > 0) {
            for (int i = 0; i < num_reqs; i++) {
                // Para cada requisição no array 'requests', esperamos sua conclusão.
                MPI_Wait(&requests[i], MPI_STATUS_IGNORE); 
            }
        }

        // 4. Computar novos valores de u_next usando os dados de u_curr (agora com bordas atualizadas)
        for (int i = 1; i <= local_n; ++i) {
            if (rank == 0 && i == 1) { // Condição de contorno global à esquerda
                u_next[i] = BC_GLOBAL_LEFT;
            } else if (rank == mpi_size - 1 && i == local_n) { // Condição de contorno global à direita
                u_next[i] = BC_GLOBAL_RIGHT;
            } else { // Pontos internos (usam células fantasmas já recebidas)
                u_next[i] = u_curr[i] + factor * (u_curr[i - 1] - 2.0 * u_curr[i] + u_curr[i + 1]);
            }
        }

        // 5. Trocar os ponteiros de u_curr e u_next para a próxima iteração
        double* temp = u_curr;
        u_curr = u_next;
        u_next = temp;
    }
}


int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    if (mpi_size != 3) {
        if (rank == 0) {
            fprintf(stderr, "Erro: Este programa foi projetado para ser executado com exatamente 3 processos MPI.\n");
            fprintf(stderr, "Por favor, execute com: mpirun -np 3 %s\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    if (N_GLOBAL_POINTS < 2) {
        if (rank == 0) fprintf(stderr, "Erro: N_GLOBAL_POINTS deve ser ao menos 2.\n");
        MPI_Finalize();
        return 1;
    }
     if (N_GLOBAL_POINTS < mpi_size && rank == 0) { // mpi_size é 3
        printf("Aviso: N_GLOBAL_POINTS (%d) é menor que o número de processos (3).\n", N_GLOBAL_POINTS);
        printf("         Alguns processos podem não ter pontos para calcular (local_n = 0).\n");
    }

    int local_n, global_offset;
    calculate_local_params(rank, mpi_size, &local_n, &global_offset);

    if (local_n == 0) {
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        return 0;
    }

    double dx = L_DOMAIN / (N_GLOBAL_POINTS - 1);
    double dt = CFL_STABILITY_FACTOR * dx * dx / ALPHA_DIFFUSIVITY;
    int num_steps = (int)(T_FINAL_TIME / dt);
    double stability_param = ALPHA_DIFFUSIVITY * dt / (dx * dx);

    if (rank == 0) {
        printf("## Simulação de Difusão de Calor 1D com MPI (para 3 processos) ##\n");
        printf("L: %.2f, T_final: %.3f, Alpha: %.4f, N_global: %d\n", L_DOMAIN, T_FINAL_TIME, ALPHA_DIFFUSIVITY, N_GLOBAL_POINTS);
        printf("dx: %e, dt: %e, Num_steps: %d\n", dx, dt, num_steps);
        printf("Parâmetro de estabilidade (Alpha*dt/dx^2): %f\n", stability_param);
        printf("Processos MPI: %d (fixo)\n", mpi_size);
        printf("---------------------------------------------------\n");
    }

    double* u_curr_nb_loop_wait = (double*)malloc((local_n + 2) * sizeof(double));
    double* u_next_nb_loop_wait = (double*)malloc((local_n + 2) * sizeof(double));

    if (!u_curr_nb_loop_wait || !u_next_nb_loop_wait) {
        fprintf(stderr, "Rank %d: Falha na alocação de memória.\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    double time_start, time_end, exec_time;

    initialize_data(u_curr_nb_loop_wait, local_n, global_offset, dx);
    MPI_Barrier(MPI_COMM_WORLD);
    time_start = MPI_Wtime();

    // Chamando a função com loop de MPI_Wait
    run_simulation_nonblocking_individual_waits(u_curr_nb_loop_wait, u_next_nb_loop_wait, local_n, rank, mpi_size, num_steps, stability_param);

    MPI_Barrier(MPI_COMM_WORLD);
    time_end = MPI_Wtime();
    exec_time = time_end - time_start;

    if (rank == 0) {
        printf("Tempo de execução (Versão MPI_Isend/Irecv + loop de MPI_Wait): %f segundos\n", exec_time);
    }

    free(u_curr_nb_loop_wait);
    free(u_next_nb_loop_wait);

    MPI_Finalize();
    return 0;
}