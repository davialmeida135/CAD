#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <cuda_runtime.h>

// --- Configuration Constants ---
#define GRID_SIZE_X 256        // Number of grid points in X dimension
#define GRID_SIZE_Y 256     // Number of grid points in Y dimension
#define GRID_SIZE_Z 256       // Number of grid points in Z dimension
#define NUM_TIME_STEPS 5      // Total number of simulation time steps
#define CELL_SPACING_X 1.0f     // Spatial step size (delta X)
#define CELL_SPACING_Y 1.0f     // Spatial step size (delta Y)
#define CELL_SPACING_Z 1.0f     // Spatial step size (delta Z)
#define TIME_STEP_INCREMENT 0.1f // Time step size (delta T)
#define DIFFUSION_COEFFICIENT 0.1f // Diffusion coefficient (Nu)
#define PERTURBATION_MAGNITUDE 10.0f // Magnitude of perturbation

// Macro to convert 3D grid coordinates to a 1D array index
#define MAP_3D_TO_1D_INDEX(x, y, z) ((z) * GRID_SIZE_Y * GRID_SIZE_X + (y) * GRID_SIZE_X + (x))

__global__ void atualiza(double* vnew, double* vold, 
    int nx, int ny, int nz, double alpha) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;
    
    if (x > 0 && x < nx-1 && y > 0 && y < ny-1 && z > 0 && z < nz-1) {
        int idx = z * ny * nx + y * nx + x;
        int xm = idx - 1;
        int xp = idx + 1;
        int ym = idx - nx;
        int yp = idx + nx;
        int zm = idx - nx * ny;
        int zp = idx + nx * ny;
        
        vnew[idx] = vold[idx] + alpha * (
            vold[xp] + vold[xm] +
            vold[yp] + vold[ym] +
            vold[zp] + vold[zm] - 6 * vold[idx]);
    }
}

float calculate_max_deviation_from_one(double* field, int nx, int ny, int nz) {
    float max_dev = 0.0f;
    for (int i = 0; i < nx * ny * nz; i++) {
        float dev = fabs(field[i] - 1.0);
        if (dev > max_dev) max_dev = dev;
    }
    return max_dev;
}

float loop_cuda(int nx, int ny, int nz, int num_steps, double alpha) {
    size_t size = nx * ny * nz * sizeof(double);
    
    // Allocate host memory
    double *h_vold = (double*)malloc(size);
    double *h_result = (double*)malloc(size);
    
    if (!h_vold || !h_result) {
        printf("Failed to allocate host memory\n");
        return -1.0f;
    }
    
    // Initialize host array: all points to 1.0
    for (int i = 0; i < nx * ny * nz; i++) {
        h_vold[i] = 1.0;
    }
    
    // Add perturbation at center
    int center_x = nx / 2;
    int center_y = ny / 2;
    int center_z = nz / 2;
    int center_idx = center_z * ny * nx + center_y * nx + center_x;
    h_vold[center_idx] += PERTURBATION_MAGNITUDE;
    
    // Allocate device memory
    double *d_vold, *d_vnew;
    cudaMalloc(&d_vold, size);
    cudaMalloc(&d_vnew, size);
    
    // Copy initial data to device
    cudaMemcpy(d_vold, h_vold, size, cudaMemcpyHostToDevice);
    
    // Define 3D block and grid
    int bx = 8, by = 8, bz = 8;
    dim3 threads(bx, by, bz);
    dim3 grid((nx + bx - 1) / bx, (ny + by - 1) / by, (nz + bz - 1) / bz);
    
    printf("Starting CUDA simulation...\n");
    printf("Grid size: %dx%dx%d\n", nx, ny, nz);
    printf("Block size: %dx%dx%d\n", bx, by, bz);
    printf("Grid dimensions: %dx%dx%d\n", grid.x, grid.y, grid.z);
    
    // Create CUDA events for timing
    cudaEvent_t start_event, stop_event;
    cudaEventCreate(&start_event);
    cudaEventCreate(&stop_event);
    
    // Record start event
    cudaEventRecord(start_event);
    
    // Main simulation loop
    for (int t = 0; t < num_steps; t++) {
        // Launch kernel
        atualiza<<<grid, threads>>>(d_vnew, d_vold, nx, ny, nz, alpha);
        
        // Check for kernel launch errors
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            printf("CUDA kernel error: %s\n", cudaGetErrorString(err));
            break;
        }
        
        // Swap pointers
        double* tmp = d_vold;
        d_vold = d_vnew;
        d_vnew = tmp;
    }
    
    // Record stop event
    cudaEventRecord(stop_event);
    
    // Wait for all kernels to complete
    cudaEventSynchronize(stop_event);
    
    // Calculate elapsed time using CUDA events
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start_event, stop_event);
    float elapsed = milliseconds / 1000.0f; // Convert to seconds
    
    // Copy result back to host
    cudaMemcpy(h_result, d_vold, size, cudaMemcpyDeviceToHost);
    
    // Calculate final deviation
    float final_max_dev = calculate_max_deviation_from_one(h_result, nx, ny, nz);
    
    printf("Total CUDA simulation time: %f seconds\n", elapsed);
    printf("Final max deviation from 1.0: %.6f\n", final_max_dev);
    printf("Center value: %.6f\n", h_result[center_idx]);
    
    // Cleanup CUDA events
    cudaEventDestroy(start_event);
    cudaEventDestroy(stop_event);
    
    // Cleanup memory
    cudaFree(d_vold);
    cudaFree(d_vnew);
    free(h_vold);
    free(h_result);
    
    return elapsed;
}

int main() {
    // Simulation parameters
    int nx = GRID_SIZE_X;
    int ny = GRID_SIZE_Y;
    int nz = GRID_SIZE_Z;
    int nt = NUM_TIME_STEPS;
    double alpha = DIFFUSION_COEFFICIENT * TIME_STEP_INCREMENT / 
                   (CELL_SPACING_X * CELL_SPACING_X);
    
    printf("## Simulação de Difusão 3D com CUDA ##\n");
    printf("Grid: %dx%dx%d\n", nx, ny, nz);
    printf("Time steps: %d\n", nt);
    printf("Alpha: %f\n", alpha);
    printf("----------------------------------------\n");
    
    // Check CUDA device
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount == 0) {
        printf("No CUDA devices found!\n");
        return 1;
    }
    
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    printf("Using device: %s\n", prop.name);
    printf("----------------------------------------\n");
    
    // Run simulation
    float exec_time = loop_cuda(nx, ny, nz, nt, alpha);
    
    if (exec_time > 0) {
        printf("----------------------------------------\n");
        printf("Simulation completed successfully!\n");
    } else {
        printf("Simulation failed!\n");
        return 1;
    }
    
    return 0;
}