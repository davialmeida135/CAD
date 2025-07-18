#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <cuda_runtime.h>

// --- Configuration Constants ---
#define GRID_SIZE_X 512        // Number of grid points in X dimension
#define GRID_SIZE_Y 512        // Number of grid points in Y dimension
#define GRID_SIZE_Z 512        // Number of grid points in Z dimension
#define NUM_TIME_STEPS 5      // Total number of simulation time steps
#define CELL_SPACING_X 1.0f     // Spatial step size (delta X)
#define CELL_SPACING_Y 1.0f     // Spatial step size (delta Y)
#define CELL_SPACING_Z 1.0f     // Spatial step size (delta Z)
#define TIME_STEP_INCREMENT 0.1f // Time step size (delta T)
#define DIFFUSION_COEFFICIENT 0.1f // Diffusion coefficient (Nu)
#define PERTURBATION_MAGNITUDE 10.0f // Magnitude of perturbation

// Macro to convert 3D grid coordinates to a 1D array index
#define MAP_3D_TO_1D_INDEX(x_coord, y_coord, z_coord) ((z_coord) * GRID_SIZE_Y * GRID_SIZE_X + (y_coord) * GRID_SIZE_X + (x_coord))

__global__ void update_temperature_field(double* temperature_field_new, double* temperature_field_old, 
    int grid_points_x, int grid_points_y, int grid_points_z, double diffusion_alpha) {

    // Calcula as coordenadas 3D da thread atual dentro do grid
    // Cada thread processará um ponto específico da grade 3D
    // Block Idx = posição do bloco no grid
    // Block Dim = número de threads por bloco
    // Thread Idx = posição da thread dentro do bloco
    int thread_x = blockIdx.x * blockDim.x + threadIdx.x;  // Posição X da thread no grid global
    int thread_y = blockIdx.y * blockDim.y + threadIdx.y;  // Posição Y da thread no grid global
    int thread_z = blockIdx.z * blockDim.z + threadIdx.z;  // Posição Z da thread no grid global

    // Verifica se a thread está dentro dos limites válidos da grade
    // Exclui as bordas (pontos nas extremidades) pois eles não têm todos os vizinhos necessários
    // Para aplicar o operador de diferenças finitas, precisamos de 6 vizinhos (±X, ±Y, ±Z)
    if (thread_x > 0 && thread_x < grid_points_x-1 &&     // Não está na borda esquerda/direita
        thread_y > 0 && thread_y < grid_points_y-1 &&     // Não está na borda frente/trás  
        thread_z > 0 && thread_z < grid_points_z-1) {     // Não está na borda superior/inferior
        
        // Converte as coordenadas 3D (x,y,z) para índice 1D no array
        // Fórmula: z*altura*largura + y*largura + x
        // Isso mapeia a grade 3D para um array linear na memória
        int current_index = thread_z * grid_points_y * grid_points_x + thread_y * grid_points_x + thread_x;
        
        // Calcula os índices dos 6 vizinhos próximos (norte, sul, leste, oeste, cima, baixo)
        int neighbor_x_minus = current_index - 1;
        int neighbor_x_plus = current_index + 1;
        int neighbor_y_minus = current_index - grid_points_x;
        int neighbor_y_plus = current_index + grid_points_x;
        int neighbor_z_minus = current_index - grid_points_x * grid_points_y;
        int neighbor_z_plus = current_index + grid_points_x * grid_points_y;
        
        // Aplica a equação
        temperature_field_new[current_index] = temperature_field_old[current_index] + diffusion_alpha * (
            temperature_field_old[neighbor_x_plus] + temperature_field_old[neighbor_x_minus] +   
            temperature_field_old[neighbor_y_plus] + temperature_field_old[neighbor_y_minus] +    
            temperature_field_old[neighbor_z_plus] + temperature_field_old[neighbor_z_minus] -   
            6 * temperature_field_old[current_index]);                                           
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
    // Destino, origem, tamanho, direção da cópia
    cudaMemcpy(device_temperature_old, host_temperature_initial, memory_size, cudaMemcpyHostToDevice);
    
    // Define 3D block and grid
    int block_size_x = 8, block_size_y = 8, block_size_z = 8;
    dim3 threads_per_block(block_size_x, block_size_y, block_size_z);
    dim3 blocks_per_grid((grid_points_x + block_size_x - 1) / block_size_x, 
                        (grid_points_y + block_size_y - 1) / block_size_y, 
                        (grid_points_z + block_size_z - 1) / block_size_z);
    
    printf("Starting CUDA simulation...\n");
    printf("Grid size: %dx%dx%d\n", grid_points_x, grid_points_y, grid_points_z);
    printf("Block size: %dx%dx%d\n", block_size_x, block_size_y, block_size_z);
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
        
        // Swap pointers
        double* temporary_pointer = device_temperature_old;
        device_temperature_old = device_temperature_new;
        device_temperature_new = temporary_pointer;
        
        // Check for kernel launch errors
        cudaError_t kernel_error = cudaGetLastError();
        if (kernel_error != cudaSuccess) {
            printf("CUDA kernel error: %s\n", cudaGetErrorString(kernel_error));
            break;
        }
            
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

int main() {
    // Simulation parameters
    int grid_points_x = GRID_SIZE_X;
    int grid_points_y = GRID_SIZE_Y;
    int grid_points_z = GRID_SIZE_Z;
    int number_of_time_steps = NUM_TIME_STEPS;
    double diffusion_alpha = DIFFUSION_COEFFICIENT * TIME_STEP_INCREMENT / 
                           (CELL_SPACING_X * CELL_SPACING_X);
    
    printf("## Simulação de Difusão 3D com CUDA ##\n");
    printf("Grid: %dx%dx%d\n", grid_points_x, grid_points_y, grid_points_z);
    printf("Time steps: %d\n", number_of_time_steps);
    printf("Alpha: %f\n", diffusion_alpha);
    printf("----------------------------------------\n");
    
    // Check CUDA device
    int cuda_device_count;
    cudaGetDeviceCount(&cuda_device_count);
    if (cuda_device_count == 0) {
        printf("No CUDA devices found!\n");
        return 1;
    }
    
    cudaDeviceProp device_properties;
    cudaGetDeviceProperties(&device_properties, 0);
    printf("Using device: %s\n", device_properties.name);
    printf("----------------------------------------\n");
    
    // Run simulation
    float execution_time = run_cuda_simulation(grid_points_x, grid_points_y, grid_points_z, number_of_time_steps, diffusion_alpha);
    
    if (execution_time > 0) {
        printf("----------------------------------------\n");
        printf("Simulation completed successfully!\n");
    } else {
        printf("Simulation failed!\n");
        return 1;
    }
    
    return 0;
}