#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <cuda_runtime.h>

// --- Configuration Constants ---
#define GRID_SIZE_X 128        // Number of grid points in X dimension
#define GRID_SIZE_Y 128        // Number of grid points in Y dimension
#define GRID_SIZE_Z 128        // Number of grid points in Z dimension
#define NUM_TIME_STEPS 5      // Total number of simulation time steps
#define CELL_SPACING_X 1.0f     // Spatial step size (delta X)
#define CELL_SPACING_Y 1.0f     // Spatial step size (delta Y)
#define CELL_SPACING_Z 1.0f     // Spatial step size (delta Z)
#define TIME_STEP_INCREMENT 0.1f // Time step size (delta T)
#define DIFFUSION_COEFFICIENT 0.1f // Diffusion coefficient (Nu)
#define PERTURBATION_MAGNITUDE 10.0f // Magnitude of perturbation
#define RADIUS 1
#define BLOCK_SIZE_X 8 // Number of threads in X dimension per block
#define BLOCK_SIZE_Y 8 // Number of threads in Y dimension per block
#define BLOCK_SIZE_Z 8 // Number of threads in Z dimension per block

// Macro to convert 3D grid coordinates to a 1D array index
#define MAP_3D_TO_1D_INDEX(x_coord, y_coord, z_coord) ((z_coord) * GRID_SIZE_Y * GRID_SIZE_X + (y_coord) * GRID_SIZE_X + (x_coord))

__global__ void update_temperature_field(double* temperature_field_new, double* temperature_field_old, 
    int grid_points_x, int grid_points_y, int grid_points_z, double diffusion_alpha) {

    __shared__ float tile[BLOCK_SIZE_X+3*RADIUS][BLOCK_SIZE_Y+3*RADIUS][BLOCK_SIZE_Z+3*RADIUS]; // Shared memory tile for block
    // Calcula as coordenadas 3D da thread atual dentro do grid
    // Cada thread processará um ponto específico da grade 3D
    // Block Idx = posição do bloco no grid
    // Block Dim = número de threads por bloco
    // Thread Idx = posição da thread dentro do bloco
    int thread_x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x;  // Posição X da thread no grid global
    int thread_y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y;  // Posição Y da thread no grid global
    int thread_z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;  // Posição Z da thread no grid global

    int tile_x = threadIdx.x + RADIUS; // Posição X da thread dentro do bloco
    int tile_y = threadIdx.y + RADIUS; // Posição Y da thread dentro do bloco
    int tile_z = threadIdx.z + RADIUS; // Posição Z da thread dentro do bloco

    // Verifica se a thread está dentro dos limites válidos da grade
    // Exclui as bordas (pontos nas extremidades) pois eles não têm todos os vizinhos necessários
    // Para aplicar o operador de diferenças finitas, precisamos de 6 vizinhos (±X, ±Y, ±Z)
    if (thread_x > 0 && thread_x < grid_points_x-1 &&     // Não está na borda esquerda/direita
        thread_y > 0 && thread_y < grid_points_y-1 &&     // Não está na borda frente/trás  
        thread_z > 0 && thread_z < grid_points_z-1) {     // Não está na borda superior/inferior
        
        
        int current_index = thread_z * grid_points_y * grid_points_x + thread_y * grid_points_x + thread_x;
        
        // Copia pixel central
        if (tile_x < BLOCK_SIZE_X && tile_y < BLOCK_SIZE_Y && tile_z < BLOCK_SIZE_Z)
            tile[tile_z][tile_y][tile_x] = temperature_field_old[current_index];

        // Copia vizinhos para a memória compartilhada
        if (tile_x - 1 >= 0 && tile_x - 1 < BLOCK_SIZE_X && tile_y < BLOCK_SIZE_Y && tile_z < BLOCK_SIZE_Z)
            tile[tile_z][tile_y][tile_x - 1] = temperature_field_old[current_index - 1]; // Vizinho à esquerda
        if (tile_x + 1 < BLOCK_SIZE_X && tile_y < BLOCK_SIZE_Y && tile_z < BLOCK_SIZE_Z)
            tile[tile_z][tile_y][tile_x + 1] = temperature_field_old[current_index + 1]; // Vizinho à direita
        if (tile_y - 1 >= 0 && tile_x < BLOCK_SIZE_X && tile_z < BLOCK_SIZE_Z)
            tile[tile_z][tile_y - 1][tile_x] = temperature_field_old[current_index - grid_points_x]; // Vizinho acima
        if (tile_y + 1 < BLOCK_SIZE_Y && tile_x < BLOCK_SIZE_X && tile_z < BLOCK_SIZE_Z)
            tile[tile_z][tile_y + 1][tile_x] = temperature_field_old[current_index + grid_points_x]; // Vizinho abaixo
        if (tile_z - 1 >= 0 && tile_x < BLOCK_SIZE_X && tile_y < BLOCK_SIZE_Y)
            tile[tile_z - 1][tile_y][tile_x] = temperature_field_old[current_index - grid_points_x * grid_points_y]; // Vizinho atrás
        if (tile_z + 1 < BLOCK_SIZE_Z && tile_x < BLOCK_SIZE_X && tile_y < BLOCK_SIZE_Y)
            tile[tile_z + 1][tile_y][tile_x] = temperature_field_old[current_index + grid_points_x * grid_points_y]; // Vizinho à frente

        __syncthreads();
        
        // Aplica a equação no tile shared
        temperature_field_new[current_index] = tile[tile_z][tile_y][tile_x] + diffusion_alpha * (
            tile[tile_z][tile_y][tile_x - 1] + tile[tile_z][tile_y][tile_x + 1] +   
            tile[tile_z][tile_y - 1][tile_x] + tile[tile_z][tile_y + 1][tile_x] +    
            tile[tile_z - 1][tile_y][tile_x] + tile[tile_z + 1][tile_y][tile_x] -   
            6 * tile[tile_z][tile_y][tile_x]);                                         
    }
}

float calculate_max_deviation_from_one(double* temperature_field, int grid_points_x, int grid_points_y, int grid_points_z) {
    float max_deviation = 0.0f;
    int total_points = grid_points_x * grid_points_y * grid_points_z;
    
    for (int point_index = 0; point_index < total_points; point_index++) {
        float current_deviation = fabs(temperature_field[point_index] - 1.0);
        if (current_deviation > max_deviation) {
            max_deviation = current_deviation;
        }
    }
    return max_deviation;
}

float run_cuda_simulation(int grid_points_x, int grid_points_y, int grid_points_z, int number_of_time_steps, double diffusion_alpha) {
    size_t memory_size = grid_points_x * grid_points_y * grid_points_z * sizeof(double);
    
    // Allocate host memory
    double *host_temperature_initial = (double*)malloc(memory_size);
    double *host_temperature_result = (double*)malloc(memory_size);
    
    if (!host_temperature_initial || !host_temperature_result) {
        printf("Failed to allocate host memory\n");
        return -1.0f;
    }
    
    // Initialize host array: all points to 1.0
    int total_grid_points = grid_points_x * grid_points_y * grid_points_z;
    for (int point_index = 0; point_index < total_grid_points; point_index++) {
        host_temperature_initial[point_index] = 1.0;
    }
    
    // Add perturbation at center
    int center_x_coordinate = grid_points_x / 2;
    int center_y_coordinate = grid_points_y / 2;
    int center_z_coordinate = grid_points_z / 2;
    int center_point_index = center_z_coordinate * grid_points_y * grid_points_x + 
                            center_y_coordinate * grid_points_x + center_x_coordinate;
    host_temperature_initial[center_point_index] += PERTURBATION_MAGNITUDE;
    
    // Allocate device memory
    double *device_temperature_old, *device_temperature_new;
    cudaMalloc(&device_temperature_old, memory_size);
    cudaMalloc(&device_temperature_new, memory_size);
    
    // Copy initial data to device
    cudaMemcpy(device_temperature_old, host_temperature_initial, memory_size, cudaMemcpyHostToDevice);
    
    // Define 3D block and grid
    dim3 threads_per_block(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z);
    dim3 blocks_per_grid((grid_points_x + BLOCK_SIZE_X - 1) / BLOCK_SIZE_X, 
                        (grid_points_y + BLOCK_SIZE_Y - 1) / BLOCK_SIZE_Y, 
                        (grid_points_z + BLOCK_SIZE_Z - 1) / BLOCK_SIZE_Z);
    
    printf("Starting CUDA simulation...\n");
    printf("Grid size: %dx%dx%d\n", grid_points_x, grid_points_y, grid_points_z);
    printf("Block size: %dx%dx%d\n", BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z);
    printf("Grid dimensions: %dx%dx%d\n", blocks_per_grid.x, blocks_per_grid.y, blocks_per_grid.z);
    
    // Create CUDA events for timing
    cudaEvent_t simulation_start_event, simulation_stop_event;
    cudaEventCreate(&simulation_start_event);
    cudaEventCreate(&simulation_stop_event);
    
    // Record start event
    cudaEventRecord(simulation_start_event);
    
    // Main simulation loop
    for (int time_step = 0; time_step < number_of_time_steps; time_step++) {
        // Launch kernel
        update_temperature_field<<<blocks_per_grid, threads_per_block>>>(
            device_temperature_new, device_temperature_old, 
            grid_points_x, grid_points_y, grid_points_z, diffusion_alpha);
        
        // Check for kernel launch errors
        cudaError_t kernel_error = cudaGetLastError();
        if (kernel_error != cudaSuccess) {
            printf("CUDA kernel error: %s\n", cudaGetErrorString(kernel_error));
            break;
        }
        
        // Swap pointers
        double* temporary_pointer = device_temperature_old;
        device_temperature_old = device_temperature_new;
        device_temperature_new = temporary_pointer;
    }
    
    // Record stop event
    cudaEventRecord(simulation_stop_event);
    
    // Wait for all kernels to complete
    cudaEventSynchronize(simulation_stop_event);
    
    // Calculate elapsed time using CUDA events
    float elapsed_milliseconds = 0;
    cudaEventElapsedTime(&elapsed_milliseconds, simulation_start_event, simulation_stop_event);
    float elapsed_seconds = elapsed_milliseconds / 1000.0f; // Convert to seconds
    
    // Copy result back to host
    cudaMemcpy(host_temperature_result, device_temperature_old, memory_size, cudaMemcpyDeviceToHost);
    
    // Calculate final deviation
    float final_maximum_deviation = calculate_max_deviation_from_one(host_temperature_result, grid_points_x, grid_points_y, grid_points_z);
    
    printf("Total CUDA simulation time: %f seconds\n", elapsed_seconds);
    printf("Final max deviation from 1.0: %.6f\n", final_maximum_deviation);
    printf("Center value: %.6f\n", host_temperature_result[center_point_index]);
    // Resultado CPU 7.388917

    // Cleanup CUDA events
    cudaEventDestroy(simulation_start_event);
    cudaEventDestroy(simulation_stop_event);
    
    // Cleanup memory
    cudaFree(device_temperature_old);
    cudaFree(device_temperature_new);
    free(host_temperature_initial);
    free(host_temperature_result);
    
    return elapsed_seconds;
}
// Adicione esta função CPU ao seu navier.cu
void update_temperature_cpu(double* temp_new, double* temp_old, 
                           int nx, int ny, int nz, double alpha) {
    for (int z = 1; z < nz - 1; z++) {
        for (int y = 1; y < ny - 1; y++) {
            for (int x = 1; x < nx - 1; x++) {
                int idx = z * ny * nx + y * nx + x;
                int idx_xm = idx - 1;
                int idx_xp = idx + 1;
                int idx_ym = idx - nx;
                int idx_yp = idx + nx;
                int idx_zm = idx - nx * ny;
                int idx_zp = idx + nx * ny;
                
                temp_new[idx] = temp_old[idx] + alpha * (
                    temp_old[idx_xm] + temp_old[idx_xp] +
                    temp_old[idx_ym] + temp_old[idx_yp] +
                    temp_old[idx_zm] + temp_old[idx_zp] -
                    6.0 * temp_old[idx]);
            }
        }
    }
}

float run_cpu_simulation(int nx, int ny, int nz, int num_steps, double alpha) {
    size_t size = nx * ny * nz * sizeof(double);
    
    double *cpu_temp_curr = (double*)malloc(size);
    double *cpu_temp_next = (double*)malloc(size);
    
    if (!cpu_temp_curr || !cpu_temp_next) {
        printf("Failed to allocate CPU memory\n");
        return -1.0f;
    }
    
    // Initialize: all points to 1.0
    for (int i = 0; i < nx * ny * nz; i++) {
        cpu_temp_curr[i] = 1.0;
    }
    
    // Add perturbation at center
    int center_idx = (nz/2) * ny * nx + (ny/2) * nx + (nx/2);
    cpu_temp_curr[center_idx] += PERTURBATION_MAGNITUDE;
    
    // Start timing
    struct timeval start, end;
    gettimeofday(&start, NULL);
    
    // Simulation loop
    for (int t = 0; t < num_steps; t++) {
        update_temperature_cpu(cpu_temp_next, cpu_temp_curr, nx, ny, nz, alpha);
        
        // Swap pointers
        double* temp = cpu_temp_curr;
        cpu_temp_curr = cpu_temp_next;
        cpu_temp_next = temp;
    }
    
    // End timing
    gettimeofday(&end, NULL);
    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1e-6;
    
    // Calculate results
    float max_dev = calculate_max_deviation_from_one(cpu_temp_curr, nx, ny, nz);
    
    printf("CPU Results:\n");
    printf("Total CPU simulation time: %f seconds\n", elapsed);
    printf("Final max deviation from 1.0: %.6f\n", max_dev);
    printf("Center value: %.6f\n", cpu_temp_curr[center_idx]);
    
    free(cpu_temp_curr);
    free(cpu_temp_next);
    
    return elapsed;
}
bool compare_results(double* gpu_result, double* cpu_result, int nx, int ny, int nz, double tolerance = 1e-6) {
    int total_points = nx * ny * nz;
    int differences = 0;
    double max_error = 0.0;
    double sum_error = 0.0;
    
    for (int i = 0; i < total_points; i++) {
        double error = fabs(gpu_result[i] - cpu_result[i]);
        if (error > tolerance) {
            differences++;
        }
        if (error > max_error) {
            max_error = error;
        }
        sum_error += error;
    }
    
    double avg_error = sum_error / total_points;
    
    printf("\n=== VALIDATION RESULTS ===\n");
    printf("Total points: %d\n", total_points);
    printf("Points with differences > %.2e: %d (%.2f%%)\n", 
           tolerance, differences, 100.0 * differences / total_points);
    printf("Maximum error: %.2e\n", max_error);
    printf("Average error: %.2e\n", avg_error);
    printf("Tolerance: %.2e\n", tolerance);
    
    bool passed = (differences == 0) || (max_error < tolerance);
    printf("Validation: %s\n", passed ? "PASSED" : "FAILED");
    
    return passed;
}
int main() {
    // Simulation parameters
    int nx = GRID_SIZE_X;
    int ny = GRID_SIZE_Y;
    int nz = GRID_SIZE_Z;
    int nt = NUM_TIME_STEPS;
    double alpha = DIFFUSION_COEFFICIENT * TIME_STEP_INCREMENT / 
                   (CELL_SPACING_X * CELL_SPACING_X);
    
    printf("## Comparação CPU vs GPU - Difusão 3D ##\n");
    printf("Grid: %dx%dx%d\n", nx, ny, nz);
    printf("Time steps: %d\n", nt);
    printf("Alpha: %f\n", alpha);
    printf("========================================\n");
    
    // Check CUDA device
    int device_count;
    cudaGetDeviceCount(&device_count);
    if (device_count == 0) {
        printf("No CUDA devices found! Running CPU only.\n");
        run_cpu_simulation(nx, ny, nz, nt, alpha);
        return 1;
    }
    
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    printf("Using GPU: %s\n", prop.name);
    printf("========================================\n");
    
    // Run CPU simulation
    printf("\n>>> Running CPU simulation...\n");
    float cpu_time = run_cpu_simulation(nx, ny, nz, nt, alpha);
    
    printf("\n>>> Running GPU simulation...\n");
    float gpu_time = run_cuda_simulation(nx, ny, nz, nt, alpha);
    
    // Performance comparison
    if (cpu_time > 0 && gpu_time > 0) {
        printf("\n=== PERFORMANCE COMPARISON ===\n");
        printf("CPU time: %.6f seconds\n", cpu_time);
        printf("GPU time: %.6f seconds\n", gpu_time);
        printf("Speedup: %.2fx\n", cpu_time / gpu_time);
    }
    
    return 0;
}