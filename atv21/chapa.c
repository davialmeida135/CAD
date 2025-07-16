#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <string.h>

// Parâmetros da Simulação 2D
#define L_DOMAIN_X 1.0               // Comprimento da chapa em X
#define L_DOMAIN_Y 1.0               // Comprimento da chapa em Y
#define T_FINAL_TIME 0.01            // Tempo total de simulação
#define ALPHA_DIFFUSIVITY 0.01       // Difusividade térmica
#define N_GLOBAL_X 10000               // Número total de pontos em X
#define N_GLOBAL_Y 10000              // Número total de pontos em Y
#define CFL_STABILITY_FACTOR 0.2     // Fator CFL para estabilidade 2D (<= 0.25)

// Condições de Contorno Globais
#define BC_GLOBAL 0.0

// Tags MPI para comunicação não bloqueante
#define TAG_DATA_GOES_UP 0    // Mensagem vai para processo acima
#define TAG_DATA_GOES_DOWN 1  // Mensagem vai para processo abaixo

// Função para calcular os parâmetros locais de cada processo (divisão de linhas)
void calculate_local_params(int rank, int mpi_size, int* local_rows_ptr, int* row_offset_ptr) {
    int base_chunk = N_GLOBAL_Y / mpi_size;
    int remainder = N_GLOBAL_Y % mpi_size;
    if (rank < remainder) {
        *local_rows_ptr = base_chunk + 1;
        *row_offset_ptr = rank * (base_chunk + 1);
    } else {
        *local_rows_ptr = base_chunk;
        *row_offset_ptr = remainder * (base_chunk + 1) + (rank - remainder) * base_chunk;
    }
}

// Adicionar no topo do arquivo
static inline int map_2d_to_1d(int i, int j, int cols) {
    return i * cols + j;
}

// Função para inicializar os dados (condição inicial u(x,y,0) = sin(pi*x/Lx)*sin(pi*y/Ly))
void initialize_data(double* local_u, int local_rows, int row_offset, double dx, double dy) {
    // Array incluindo células fantasmas: [linha_fantasma_superior][linhas_locais][linha_fantasma_inferior]
    // Cada linha tem N_GLOBAL_X pontos
    
    // Inicializar células fantasmas
    for (int j = 0; j < N_GLOBAL_X; j++) {
        local_u[map_2d_to_1d(0, j, N_GLOBAL_X)] = 0.0; // Fantasma superior
        local_u[map_2d_to_1d(local_rows + 1, j, N_GLOBAL_X)] = 0.0; // Fantasma inferior
    }
    
    // Inicializar pontos locais
    for (int i = 0; i < local_rows; i++) {
        for (int j = 0; j < N_GLOBAL_X; j++) {
            double x = j * dx;
            double y = (row_offset + i) * dy;
            local_u[map_2d_to_1d(i + 1, j, N_GLOBAL_X)] = sin(M_PI * x / L_DOMAIN_X) * sin(M_PI * y / L_DOMAIN_Y);
        }
    }
}

// Versão com MPI_Isend/MPI_Irecv e MPI_Wait para matriz 2D
void run_simulation_2d_heat(double* u_curr, double* u_next, int local_rows, int rank, int mpi_size, int num_steps, double factor_x, double factor_y) {
    MPI_Request requests[4]; // 2 envios e 2 recebimentos por processo
    int num_reqs;

    //Não pode ser paralelizado pois próximo passo depende do resultado do passo anterior
    
    for (int t = 0; t < num_steps; ++t) {
        num_reqs = 0;
        
        // 1. Iniciar comunicações não-bloqueantes
        if (rank > 0) {
            MPI_Irecv(&u_curr[map_2d_to_1d(0, 0, N_GLOBAL_X)], N_GLOBAL_X, MPI_DOUBLE, 
                     rank - 1, TAG_DATA_GOES_DOWN, MPI_COMM_WORLD, &requests[num_reqs++]);
        }
        if (rank < mpi_size - 1) {
            MPI_Irecv(&u_curr[map_2d_to_1d(local_rows + 1, 0, N_GLOBAL_X)], N_GLOBAL_X, MPI_DOUBLE, 
                     rank + 1, TAG_DATA_GOES_UP, MPI_COMM_WORLD, &requests[num_reqs++]);
        }
        if (rank < mpi_size - 1) {
            MPI_Isend(&u_curr[map_2d_to_1d(local_rows, 0, N_GLOBAL_X)], N_GLOBAL_X, MPI_DOUBLE, 
                     rank + 1, TAG_DATA_GOES_DOWN, MPI_COMM_WORLD, &requests[num_reqs++]);
        }
        if (rank > 0) {
            MPI_Isend(&u_curr[map_2d_to_1d(1, 0, N_GLOBAL_X)], N_GLOBAL_X, MPI_DOUBLE, 
                     rank - 1, TAG_DATA_GOES_UP, MPI_COMM_WORLD, &requests[num_reqs++]);
        }
        
        // 2. Calcular pontos internos (que não dependem das células fantasmas)
        #pragma omp parallel for collapse(2) schedule(guided)
        for (int i = 2; i <= local_rows - 1; ++i) {
            for (int j = 1; j < N_GLOBAL_X - 1; ++j) {
                int idx = map_2d_to_1d(i, j, N_GLOBAL_X);
                int idx_up = map_2d_to_1d(i - 1, j, N_GLOBAL_X);
                int idx_down = map_2d_to_1d(i + 1, j, N_GLOBAL_X);
                int idx_left = map_2d_to_1d(i, j - 1, N_GLOBAL_X);
                int idx_right = map_2d_to_1d(i, j + 1, N_GLOBAL_X);
                
                u_next[idx] = u_curr[idx] + 
                             factor_y * (u_curr[idx_up] - 2.0 * u_curr[idx] + u_curr[idx_down]) +
                             factor_x * (u_curr[idx_left] - 2.0 * u_curr[idx] + u_curr[idx_right]);
            }
        }

        // 3. Esperar comunicações
        if (num_reqs > 0) {
            MPI_Waitall(num_reqs, requests, MPI_STATUSES_IGNORE);  // MUDANÇA: usar Waitall
        }

        // 5. Calcular pontos das bordas que dependem das células fantasmas
        // Primeira linha local (i = 1)
        if (local_rows >= 1) {
            for (int j = 1; j < N_GLOBAL_X - 1; ++j) {
                int idx = map_2d_to_1d(1, j, N_GLOBAL_X);
                
                if (rank == 0) { // Contorno global superior
                    int idx_down = map_2d_to_1d(2, j, N_GLOBAL_X);
                    int idx_left = map_2d_to_1d(1, j - 1, N_GLOBAL_X);
                    int idx_right = map_2d_to_1d(1, j + 1, N_GLOBAL_X);
                    
                    u_next[idx] = u_curr[idx] + 
                                 factor_y * (BC_GLOBAL - 2.0 * u_curr[idx] + u_curr[idx_down]) +
                                 factor_x * (u_curr[idx_left] - 2.0 * u_curr[idx] + u_curr[idx_right]);
                } else { // Usa célula fantasma superior
                    int idx_up = map_2d_to_1d(0, j, N_GLOBAL_X);
                    int idx_down = map_2d_to_1d(2, j, N_GLOBAL_X);
                    int idx_left = map_2d_to_1d(1, j - 1, N_GLOBAL_X);
                    int idx_right = map_2d_to_1d(1, j + 1, N_GLOBAL_X);
                    
                    u_next[idx] = u_curr[idx] + 
                                 factor_y * (u_curr[idx_up] - 2.0 * u_curr[idx] + u_curr[idx_down]) +
                                 factor_x * (u_curr[idx_left] - 2.0 * u_curr[idx] + u_curr[idx_right]);
                }
            }
        }

        // Última linha local (i = local_rows)
        if (local_rows > 1) {
            for (int j = 1; j < N_GLOBAL_X - 1; ++j) {
                int idx = map_2d_to_1d(local_rows, j, N_GLOBAL_X);
                
                if (rank == mpi_size - 1) { // Contorno global inferior
                    int idx_up = map_2d_to_1d(local_rows - 1, j, N_GLOBAL_X);
                    int idx_left = map_2d_to_1d(local_rows, j - 1, N_GLOBAL_X);
                    int idx_right = map_2d_to_1d(local_rows, j + 1, N_GLOBAL_X);
                    
                    u_next[idx] = u_curr[idx] + 
                                 factor_y * (u_curr[idx_up] - 2.0 * u_curr[idx] + BC_GLOBAL) +
                                 factor_x * (u_curr[idx_left] - 2.0 * u_curr[idx] + u_curr[idx_right]);
                } else { // Usa célula fantasma inferior
                    int idx_up = map_2d_to_1d(local_rows - 1, j, N_GLOBAL_X);
                    int idx_down = map_2d_to_1d(local_rows + 1, j, N_GLOBAL_X);
                    int idx_left = map_2d_to_1d(local_rows, j - 1, N_GLOBAL_X);
                    int idx_right = map_2d_to_1d(local_rows, j + 1, N_GLOBAL_X);
                    
                    u_next[idx] = u_curr[idx] + 
                                 factor_y * (u_curr[idx_up] - 2.0 * u_curr[idx] + u_curr[idx_down]) +
                                 factor_x * (u_curr[idx_left] - 2.0 * u_curr[idx] + u_curr[idx_right]);
                }
            }
        }

        // 6. Aplicar condições de contorno nas bordas laterais (x = 0 e x = L_DOMAIN_X)
        for (int i = 1; i <= local_rows; ++i) {
            u_next[map_2d_to_1d(i, 0, N_GLOBAL_X)] = BC_GLOBAL; // Borda esquerda
            u_next[map_2d_to_1d(i, N_GLOBAL_X - 1, N_GLOBAL_X)] = BC_GLOBAL; // Borda direita
        }

        // 7. Trocar os ponteiros de u_curr e u_next para a próxima iteração
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

    if (N_GLOBAL_Y < mpi_size) {
        if (rank == 0) {
            fprintf(stderr, "Erro: Número de linhas (%d) deve ser >= número de processos (%d).\n", 
                   N_GLOBAL_Y, mpi_size);
        }
        MPI_Finalize();
        return 1;
    }

    int local_rows, row_offset;
    calculate_local_params(rank, mpi_size, &local_rows, &row_offset);

    if (local_rows == 0) {
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        return 0;
    }

    double dx = L_DOMAIN_X / (N_GLOBAL_X - 1);
    double dy = L_DOMAIN_Y / (N_GLOBAL_Y - 1);
    double dt = CFL_STABILITY_FACTOR * fmin(dx * dx, dy * dy) / ALPHA_DIFFUSIVITY;
    int num_steps = 200;
    double factor_x = ALPHA_DIFFUSIVITY * dt / (dx * dx);
    double factor_y = ALPHA_DIFFUSIVITY * dt / (dy * dy);

    if (rank == 0) {
        printf("## Simulação de Difusão de Calor 2D com MPI ##\n");
        printf("Domínio: %.2f x %.2f, T_final: %.3f, Alpha: %.4f\n", 
               L_DOMAIN_X, L_DOMAIN_Y, T_FINAL_TIME, ALPHA_DIFFUSIVITY);
        printf("Grid: %d x %d, dx: %e, dy: %e, dt: %e\n", 
               N_GLOBAL_X, N_GLOBAL_Y, dx, dy, dt);
        printf("Num_steps: %d, Factor_x: %f, Factor_y: %f\n", 
               num_steps, factor_x, factor_y);
        printf("Processos MPI: %d\n", mpi_size);
        printf("---------------------------------------------------\n");
    }

    // Alocar memória: (local_rows + 2) linhas x N_GLOBAL_X colunas
    // +2 para as células fantasmas (superior e inferior)
    int total_size = (local_rows + 2) * N_GLOBAL_X;
    double* u_curr = (double*)malloc(total_size * sizeof(double));
    double* u_next = (double*)malloc(total_size * sizeof(double));

    if (!u_curr || !u_next) {
        fprintf(stderr, "Rank %d: Falha na alocação de memória.\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    double time_start, time_end, exec_time;

    initialize_data(u_curr, local_rows, row_offset, dx, dy);
    MPI_Barrier(MPI_COMM_WORLD);
    time_start = MPI_Wtime();

    run_simulation_2d_heat(u_curr, u_next, local_rows, rank, mpi_size, num_steps, factor_x, factor_y);

    MPI_Barrier(MPI_COMM_WORLD);
    time_end = MPI_Wtime();
    exec_time = time_end - time_start;

    if (rank == 0) {
        printf("Tempo de execução (Simulação 2D Heat): %f segundos\n", exec_time);
    }

    free(u_curr);
    free(u_next);

    MPI_Finalize();
    return 0;
}