#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <cuda_runtime.h>
#include <omp.h>
// --- Configuration Constants ---
#define GRID_SIZE_X 512        // Number of grid points in X dimension
#define GRID_SIZE_Y 512        // Number of grid points in Y dimension
#define GRID_SIZE_Z 512        // Number of grid points in Z dimension
#define NUM_TIME_STEPS 10      // Total number of simulation time steps
#define CELL_SPACING_X 1.0f     // Spatial step size (delta X)
#define CELL_SPACING_Y 1.0f     // Spatial step size (delta Y)
#define CELL_SPACING_Z 1.0f     // Spatial step size (delta Z)
#define TIME_STEP_INCREMENT 0.1f // Time step size (delta T)
#define DIFFUSION_COEFFICIENT 0.5f // Diffusion coefficient (Nu)
#define PERTURBATION_MAGNITUDE 100.0f // Magnitude of perturbation
#define RADIUS 1
#define BLOCK_SIZE_X 8 // Number of threads in X dimension per block
#define BLOCK_SIZE_Y 8 // Number of threads in Y dimension per block
#define BLOCK_SIZE_Z 8 // Number of threads in Z dimension per block

// Macro to convert 3D grid coordinates to a 1D array index
#define MAP_3D_TO_1D_INDEX(x_coord, y_coord, z_coord) ((z_coord) * GRID_SIZE_Y * GRID_SIZE_X + (y_coord) * GRID_SIZE_X + (x_coord))

__global__ void update_temperature_field(double* temperature_field_new, double* temperature_field_old, 
    int grid_points_x, int grid_points_y, int grid_points_z, double diffusion_alpha) {

    __shared__ double tile[BLOCK_SIZE_Z + 2*RADIUS][BLOCK_SIZE_Y + 2*RADIUS][BLOCK_SIZE_X + 2*RADIUS];
    
    // Coordenadas da thread (e do ponto no espaço)
    int thread_x = blockIdx.x * blockDim.x + threadIdx.x;
    int thread_y = blockIdx.y * blockDim.y + threadIdx.y;
    int thread_z = blockIdx.z * blockDim.z + threadIdx.z;
    
    // Posição da thread dentro do bloco
    int local_x = threadIdx.x;
    int local_y = threadIdx.y;
    int local_z = threadIdx.z;
    
    // Posição do tile (thread+raio)
    int tile_x = local_x + RADIUS;
    int tile_y = local_y + RADIUS;
    int tile_z = local_z + RADIUS;
    
    // Carregar o ponto da thread
    if (thread_x < grid_points_x && thread_y < grid_points_y && thread_z < grid_points_z) {
        int global_idx = thread_z * grid_points_y * grid_points_x + thread_y * grid_points_x + thread_x;
        tile[tile_z][tile_y][tile_x] = temperature_field_old[global_idx];
    } else {
        tile[tile_z][tile_y][tile_x] = 1.0;
    }
    
    // Carregar os halos
    // Eixo X
    // Apenas a primeira thread de cada dimensão carregam os halos, não a primeira e a ultima
    if (local_x < RADIUS) {
        // Da esquerda (posicao - raio)
        int halo_x = thread_x - RADIUS;
        if (halo_x >= 0 && thread_y < grid_points_y && thread_z < grid_points_z) {
            int halo_idx = thread_z * grid_points_y * grid_points_x + thread_y * grid_points_x + halo_x;
            tile[tile_z][tile_y][local_x] = temperature_field_old[halo_idx];
        } else {
            tile[tile_z][tile_y][local_x] = 1.0;
        }
        
        // Da dureita (posicao + tamanho)
        int right_x = thread_x + blockDim.x;
        if (right_x < grid_points_x && thread_y < grid_points_y && thread_z < grid_points_z) {
            int right_idx = thread_z * grid_points_y * grid_points_x + thread_y * grid_points_x + right_x;
            tile[tile_z][tile_y][local_x + blockDim.x + RADIUS] = temperature_field_old[right_idx];
        } else {
            tile[tile_z][tile_y][local_x + blockDim.x + RADIUS] = 1.0;
        }
    }
    // Eixo Y
    if (local_y < RADIUS) {
        // De frente (posicao - raio)
        int halo_y = thread_y - RADIUS;
        if (halo_y >= 0 && thread_x < grid_points_x && thread_z < grid_points_z) {
            int halo_idx = thread_z * grid_points_y * grid_points_x + halo_y * grid_points_x + thread_x;
            tile[tile_z][local_y][tile_x] = temperature_field_old[halo_idx];
        } else {
            tile[tile_z][local_y][tile_x] = 1.0;
        }
        
        // De trás (posicao + tamanho)
        int back_y = thread_y + blockDim.y;
        if (back_y < grid_points_y && thread_x < grid_points_x && thread_z < grid_points_z) {
            int back_idx = thread_z * grid_points_y * grid_points_x + back_y * grid_points_x + thread_x;
            tile[tile_z][local_y + blockDim.y + RADIUS][tile_x] = temperature_field_old[back_idx];
        } else {
            tile[tile_z][local_y + blockDim.y + RADIUS][tile_x] = 1.0;
        }
    }
    // Eixo Z
    if (local_z < RADIUS) {
        // De cima (posicao - raio)
        int halo_z = thread_z - RADIUS;
        if (halo_z >= 0 && thread_x < grid_points_x && thread_y < grid_points_y) {
            int halo_idx = halo_z * grid_points_y * grid_points_x + thread_y * grid_points_x + thread_x;
            tile[local_z][tile_y][tile_x] = temperature_field_old[halo_idx];
        } else {
            tile[local_z][tile_y][tile_x] = 1.0;
        }
        
        // De baixo (posicao - raio)
        int top_z = thread_z + blockDim.z;
        if (top_z < grid_points_z && thread_x < grid_points_x && thread_y < grid_points_y) {
            int top_idx = top_z * grid_points_y * grid_points_x + thread_y * grid_points_x + thread_x;
            tile[local_z + blockDim.z + RADIUS][tile_y][tile_x] = temperature_field_old[top_idx];
        } else {
            tile[local_z + blockDim.z + RADIUS][tile_y][tile_x] = 1.0;
        }
    }
    __syncthreads();
    // Calcular somente para pontos internos (nao fazer para os halos)
    if (thread_x > 0 && thread_x < grid_points_x-1 &&
        thread_y > 0 && thread_y < grid_points_y-1 &&
        thread_z > 0 && thread_z < grid_points_z-1) {
        
        int global_idx = thread_z * grid_points_y * grid_points_x + thread_y * grid_points_x + thread_x;
        
        double center = tile[tile_z][tile_y][tile_x];
        double neighbors = tile[tile_z][tile_y][tile_x-1] + tile[tile_z][tile_y][tile_x+1] +  // direita - esquerda
                          tile[tile_z][tile_y-1][tile_x] + tile[tile_z][tile_y+1][tile_x] +  // cima - baixo  
                          tile[tile_z-1][tile_y][tile_x] + tile[tile_z+1][tile_y][tile_x];   // frente - trás
        
        temperature_field_new[global_idx] = center + diffusion_alpha * (neighbors - 6.0 * center);
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
    #pragma omp parallel for num_threads(16)
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


// Adicione estas funções antes do main()

// Estrutura para armazenar resultados da simulação
typedef struct {
    double* temperature_field;
    float execution_time;
    float max_deviation;
    double center_value;
} SimulationResult;

// Função para executar simulação CPU e retornar resultados
SimulationResult run_cpu_simulation_with_results(int nx, int ny, int nz, int num_steps, double alpha) {
    SimulationResult result = {0};
    size_t size = nx * ny * nz * sizeof(double);
    
    result.temperature_field = (double*)malloc(size);
    double *cpu_temp_next = (double*)malloc(size);
    
    if (!result.temperature_field || !cpu_temp_next) {
        printf("Failed to allocate CPU memory\n");
        result.execution_time = -1.0f;
        return result;
    }
    
    // Initialize: all points to 1.0
    for (int i = 0; i < nx * ny * nz; i++) {
        result.temperature_field[i] = 1.0;
    }
    
    // Add perturbation at center
    int center_idx = (nz/2) * ny * nx + (ny/2) * nx + (nx/2);
    result.temperature_field[center_idx] += PERTURBATION_MAGNITUDE;
    
    // Start timing
    struct timeval start, end;
    gettimeofday(&start, NULL);
    
    // Simulation loop
    for (int t = 0; t < num_steps; t++) {
        update_temperature_cpu(cpu_temp_next, result.temperature_field, nx, ny, nz, alpha);
        
        // Swap pointers
        double* temp = result.temperature_field;
        result.temperature_field = cpu_temp_next;
        cpu_temp_next = temp;
    }
    
    // End timing
    gettimeofday(&end, NULL);
    result.execution_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1e-6;
    
    // Calculate results
    result.max_deviation = calculate_max_deviation_from_one(result.temperature_field, nx, ny, nz);
    result.center_value = result.temperature_field[center_idx];
    
    printf("CPU Results:\n");
    printf("Total CPU simulation time: %f seconds\n", result.execution_time);
    printf("Final max deviation from 1.0: %.6f\n", result.max_deviation);
    printf("Center value: %.6f\n", result.center_value);
    
    free(cpu_temp_next);
    return result;
}

// Função para executar simulação GPU e retornar resultados
SimulationResult run_cuda_simulation_with_results(int grid_points_x, int grid_points_y, int grid_points_z, int number_of_time_steps, double diffusion_alpha) {
    SimulationResult result = {0};
    size_t memory_size = grid_points_x * grid_points_y * grid_points_z * sizeof(double);
    
    // Allocate host memory
    double *host_temperature_initial = (double*)malloc(memory_size);
    result.temperature_field = (double*)malloc(memory_size);
    
    if (!host_temperature_initial || !result.temperature_field) {
        printf("Failed to allocate host memory\n");
        result.execution_time = -1.0f;
        return result;
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
    result.execution_time = elapsed_milliseconds / 1000.0f; // Convert to seconds
    
    // Copy result back to host
    cudaMemcpy(result.temperature_field, device_temperature_old, memory_size, cudaMemcpyDeviceToHost);
    
    // Calculate final deviation
    result.max_deviation = calculate_max_deviation_from_one(result.temperature_field, grid_points_x, grid_points_y, grid_points_z);
    result.center_value = result.temperature_field[center_point_index];
    
    printf("Total CUDA simulation time: %f seconds\n", result.execution_time);
    printf("Final max deviation from 1.0: %.6f\n", result.max_deviation);
    printf("Center value: %.6f\n", result.center_value);
    
    // Cleanup CUDA events
    cudaEventDestroy(simulation_start_event);
    cudaEventDestroy(simulation_stop_event);
    
    // Cleanup memory
    cudaFree(device_temperature_old);
    cudaFree(device_temperature_new);
    free(host_temperature_initial);
    
    return result;
}

// Função completa de validação com múltiplas métricas
typedef struct {
    bool validation_passed;
    int total_points;
    int points_with_differences;
    double max_absolute_error;
    double mean_absolute_error;
    double root_mean_square_error;
    double relative_error_percent;
    double correlation_coefficient;
    double tolerance_used;
} ValidationMetrics;

ValidationMetrics validate_simulation_results(SimulationResult cpu_result, SimulationResult gpu_result, 
                                             int nx, int ny, int nz, double tolerance = 1e-6) {
    ValidationMetrics metrics = {0};
    int total_points = nx * ny * nz;
    metrics.total_points = total_points;
    metrics.tolerance_used = tolerance;
    
    double sum_absolute_error = 0.0;
    double sum_squared_error = 0.0;
    double sum_cpu = 0.0;
    double sum_gpu = 0.0;
    double sum_cpu_squared = 0.0;
    double sum_gpu_squared = 0.0;
    double sum_cpu_gpu = 0.0;
    
    // Calculate all error metrics in one pass
    for (int i = 0; i < total_points; i++) {
        double cpu_val = cpu_result.temperature_field[i];
        double gpu_val = gpu_result.temperature_field[i];
        
        double absolute_error = fabs(gpu_val - cpu_val);
        
        // Count points exceeding tolerance
        if (absolute_error > tolerance) {
            metrics.points_with_differences++;
        }
        
        // Track maximum error
        if (absolute_error > metrics.max_absolute_error) {
            metrics.max_absolute_error = absolute_error;
        }
        
        // Accumulate for various metrics
        sum_absolute_error += absolute_error;
        sum_squared_error += absolute_error * absolute_error;
        
        // For correlation coefficient
        sum_cpu += cpu_val;
        sum_gpu += gpu_val;
        sum_cpu_squared += cpu_val * cpu_val;
        sum_gpu_squared += gpu_val * gpu_val;
        sum_cpu_gpu += cpu_val * gpu_val;
    }
    
    // Calculate derived metrics
    metrics.mean_absolute_error = sum_absolute_error / total_points;
    metrics.root_mean_square_error = sqrt(sum_squared_error / total_points);
    
    // Relative error as percentage
    double mean_cpu = sum_cpu / total_points;
    if (mean_cpu != 0.0) {
        metrics.relative_error_percent = (metrics.mean_absolute_error / fabs(mean_cpu)) * 100.0;
    }
    
    // Correlation coefficient (Pearson)
    double mean_cpu_calc = sum_cpu / total_points;
    double mean_gpu_calc = sum_gpu / total_points;
    
    double numerator = sum_cpu_gpu - total_points * mean_cpu_calc * mean_gpu_calc;
    double denominator = sqrt((sum_cpu_squared - total_points * mean_cpu_calc * mean_cpu_calc) * 
                             (sum_gpu_squared - total_points * mean_gpu_calc * mean_gpu_calc));
    
    if (denominator != 0.0) {
        metrics.correlation_coefficient = numerator / denominator;
    } else {
        metrics.correlation_coefficient = 1.0; // Perfect correlation if both are constant
    }
    
    // Determine if validation passed
    metrics.validation_passed = (metrics.max_absolute_error < tolerance) && 
                               (metrics.correlation_coefficient > 0.99999); // Very high correlation required
    
    return metrics;
}

void print_validation_report(ValidationMetrics metrics, SimulationResult cpu_result, SimulationResult gpu_result) {
    printf("\n=== VALIDATION REPORT ===\n");
    printf("========================================\n");
    
    printf("Dataset Statistics:\n");
    printf("  Total points analyzed: %d\n", metrics.total_points);
    printf("  Tolerance threshold: %.2e\n", metrics.tolerance_used);
    printf("\n");
    
    printf("Error Analysis:\n");
    printf("  Maximum absolute error: %.2e\n", metrics.max_absolute_error);
    printf("  Mean absolute error: %.2e\n", metrics.mean_absolute_error);
    printf("  Root mean square error: %.2e\n", metrics.root_mean_square_error);
    printf("  Relative error: %.6f%%\n", metrics.relative_error_percent);
    printf("  Points exceeding tolerance: %d (%.4f%%)\n", 
           metrics.points_with_differences, 
           100.0 * metrics.points_with_differences / metrics.total_points);
    printf("\n");
    
    printf("Statistical Correlation:\n");
    printf("  Pearson correlation coefficient: %.10f\n", metrics.correlation_coefficient);
    printf("\n");
    
    printf("Key Value Comparison:\n");
    printf("  CPU center value: %.10f\n", cpu_result.center_value);
    printf("  GPU center value: %.10f\n", gpu_result.center_value);
    printf("  Center value difference: %.2e\n", fabs(cpu_result.center_value - gpu_result.center_value));
    printf("  CPU max deviation: %.10f\n", cpu_result.max_deviation);
    printf("  GPU max deviation: %.10f\n", gpu_result.max_deviation);
    printf("  Max deviation difference: %.2e\n", fabs(cpu_result.max_deviation - gpu_result.max_deviation));
    printf("\n");
    
    printf("Performance Comparison:\n");
    printf("  CPU execution time: %.6f seconds\n", cpu_result.execution_time);
    printf("  GPU execution time: %.6f seconds\n", gpu_result.execution_time);
    if (cpu_result.execution_time > 0 && gpu_result.execution_time > 0) {
        printf("  Speedup (CPU/GPU): %.2fx\n", cpu_result.execution_time / gpu_result.execution_time);
    }
    printf("\n");
    
    
    if (!metrics.validation_passed) {
        printf("\nPossible causes of validation failure:\n");
        if (metrics.max_absolute_error >= metrics.tolerance_used) {
            printf("  - Numerical errors exceed tolerance\n");
        }
        if (metrics.correlation_coefficient < 0.99999) {
            printf("  - Low correlation between CPU and GPU results\n");
        }
        printf("  - Consider adjusting tolerance or checking implementation\n");
    }
    
    printf("========================================\n");
}

// Função para salvar resultados detalhados (opcional)
void save_detailed_comparison(SimulationResult cpu_result, SimulationResult gpu_result, 
                             ValidationMetrics metrics, int nx, int ny, int nz) {
    FILE *file = fopen("validation_report.txt", "w");
    if (!file) {
        printf("Warning: Could not create detailed report file\n");
        return;
    }
    
    fprintf(file, "CUDA Validation Report\n");
    fprintf(file, "======================\n\n");
    fprintf(file, "Grid dimensions: %dx%dx%d\n", nx, ny, nz);
    fprintf(file, "Total points: %d\n", metrics.total_points);
    
    fprintf(file, "Error Metrics:\n");
    fprintf(file, "Max absolute error: %.2e\n", metrics.max_absolute_error);
    fprintf(file, "Mean absolute error: %.2e\n", metrics.mean_absolute_error);
    fprintf(file, "RMSE: %.2e\n", metrics.root_mean_square_error);
    fprintf(file, "Correlation: %.10f\n", metrics.correlation_coefficient);
    
    fprintf(file, "Performance:\n");
    fprintf(file, "CPU time: %.6f s\n", cpu_result.execution_time);
    fprintf(file, "GPU time: %.6f s\n", gpu_result.execution_time);
    if (cpu_result.execution_time > 0 && gpu_result.execution_time > 0) {
        fprintf(file, "Speedup: %.2fx\n", cpu_result.execution_time / gpu_result.execution_time);
    }
    
    fclose(file);
    printf("Detailed report saved to: validation_report.txt\n");
}

// Modificar o main() para usar as novas funções
int main() {
    // Simulation parameters
    int nx = GRID_SIZE_X;
    int ny = GRID_SIZE_Y;
    int nz = GRID_SIZE_Z;
    int nt = NUM_TIME_STEPS;
    double alpha = DIFFUSION_COEFFICIENT * TIME_STEP_INCREMENT / 
                   (CELL_SPACING_X * CELL_SPACING_X);
    
    printf("## Validação CPU vs GPU - Difusão 3D ##\n");
    printf("Grid: %dx%dx%d\n", nx, ny, nz);
    printf("Time steps: %d\n", nt);
    printf("Alpha: %f\n", alpha);
    printf("========================================\n");
    
    // Check CUDA device
    int device_count;
    cudaGetDeviceCount(&device_count);
    if (device_count == 0) {
        printf("No CUDA devices found! Running CPU only.\n");
        SimulationResult cpu_only = run_cpu_simulation_with_results(nx, ny, nz, nt, alpha);
        free(cpu_only.temperature_field);
        return 1;
    }
    
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    printf("Using GPU: %s\n", prop.name);
    printf("========================================\n");
    
    // Run both simulations
    printf("\n>>> Running CPU simulation...\n");
    SimulationResult cpu_result = run_cpu_simulation_with_results(nx, ny, nz, nt, alpha);
    
    printf("\n>>> Running GPU simulation...\n");
    SimulationResult gpu_result = run_cuda_simulation_with_results(nx, ny, nz, nt, alpha);
    
    // Validate results
    if (cpu_result.execution_time > 0 && gpu_result.execution_time > 0) {
        ValidationMetrics validation = validate_simulation_results(cpu_result, gpu_result, nx, ny, nz, 1e-6);
        print_validation_report(validation, cpu_result, gpu_result);
        save_detailed_comparison(cpu_result, gpu_result, validation, nx, ny, nz);
    }
    
    // Cleanup
    free(cpu_result.temperature_field);
    free(gpu_result.temperature_field);
    
    return 0;
}