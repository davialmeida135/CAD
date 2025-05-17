#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>

// --- Configuration Constants ---
#define GRID_SIZE_X 128         // Number of grid points in X dimension
#define GRID_SIZE_Y 128         // Number of grid points in Y dimension
#define GRID_SIZE_Z 128         // Number of grid points in Z dimension
#define NUM_TIME_STEPS 100      // Total number of simulation time steps
#define CELL_SPACING_X 1.0f     // Spatial step size (delta X)
#define CELL_SPACING_Y 1.0f     // Spatial step size (delta Y)
#define CELL_SPACING_Z 1.0f     // Spatial step size (delta Z)
#define TIME_STEP_INCREMENT 0.1f // Time step size (delta T)
#define DIFFUSION_COEFFICIENT 0.1f // Diffusion coefficient (Nu)
#define PERTURBATION_MAGNITUDE 10.0f // Magnitude of perturbation

// Macro to convert 3D grid coordinates to a 1D array index
#define MAP_3D_TO_1D_INDEX(ix, iy, iz) \
    ((ix) * GRID_SIZE_Y * GRID_SIZE_Z + (iy) * GRID_SIZE_Z + (iz))

float calculate_laplacian_periodic(float *scalar_field_data, int x_coord, int y_coord, int z_coord) {

    // Calculate neighbor indices with periodic wrap-around
    int x_plus_1 = (x_coord + 1) % GRID_SIZE_X;
    int x_minus_1 = (x_coord - 1 + GRID_SIZE_X) % GRID_SIZE_X;
    int y_plus_1 = (y_coord + 1) % GRID_SIZE_Y;
    int y_minus_1 = (y_coord - 1 + GRID_SIZE_Y) % GRID_SIZE_Y;
    int z_plus_1 = (z_coord + 1) % GRID_SIZE_Z;
    int z_minus_1 = (z_coord - 1 + GRID_SIZE_Z) % GRID_SIZE_Z;

    // Get values at current point and neighbors
    float val_center = scalar_field_data[MAP_3D_TO_1D_INDEX(x_coord, y_coord, z_coord)];
    float val_x_plus_1 = scalar_field_data[MAP_3D_TO_1D_INDEX(x_plus_1, y_coord, z_coord)];
    float val_x_minus_1 = scalar_field_data[MAP_3D_TO_1D_INDEX(x_minus_1, y_coord, z_coord)];
    float val_y_plus_1 = scalar_field_data[MAP_3D_TO_1D_INDEX(x_coord, y_plus_1, z_coord)];
    float val_y_minus_1 = scalar_field_data[MAP_3D_TO_1D_INDEX(x_coord, y_minus_1, z_coord)];
    float val_z_plus_1 = scalar_field_data[MAP_3D_TO_1D_INDEX(x_coord, y_coord, z_plus_1)];
    float val_z_minus_1 = scalar_field_data[MAP_3D_TO_1D_INDEX(x_coord, y_coord, z_minus_1)];

    // Discrete Laplacian formula
    float laplacian_val =
        (val_x_plus_1 - 2.0f * val_center + val_x_minus_1) / (CELL_SPACING_X * CELL_SPACING_X) +
        (val_y_plus_1 - 2.0f * val_center + val_y_minus_1) / (CELL_SPACING_Y * CELL_SPACING_Y) +
        (val_z_plus_1 - 2.0f * val_center + val_z_minus_1) / (CELL_SPACING_Z * CELL_SPACING_Z);

    return laplacian_val;
}

// Function to advance the simulation by one time step for a given scalar field
void perform_time_step(float *current_field_data, float *next_field_data) {
    for (int ix = 0; ix < GRID_SIZE_X; ix++) {
        for (int iy = 0; iy < GRID_SIZE_Y; iy++) {
            for (int iz = 0; iz < GRID_SIZE_Z; iz++) {
                int flat_index = MAP_3D_TO_1D_INDEX(ix, iy, iz);
                float laplacian_at_point = calculate_laplacian_periodic(current_field_data, ix, iy, iz);
                next_field_data[flat_index] = current_field_data[flat_index] +
                                           DIFFUSION_COEFFICIENT * TIME_STEP_INCREMENT * laplacian_at_point;
            }
        }
    }
}

// Function to calculate the maximum deviation of a field from a reference value (1.0f)
float calculate_max_deviation_from_one(float *scalar_field_data) {
    float max_abs_deviation = 0.0f;
    int total_grid_points = GRID_SIZE_X * GRID_SIZE_Y * GRID_SIZE_Z;
    for (int i = 0; i < total_grid_points; i++) {
        float current_deviation = fabsf(scalar_field_data[i] - 1.0f);
        if (current_deviation > max_abs_deviation) {
            max_abs_deviation = current_deviation;
        }
    }
    return max_abs_deviation;
}

int main() {
    struct timeval start, end, step_start, step_end;

    size_t single_field_mem_size = GRID_SIZE_X * GRID_SIZE_Y * GRID_SIZE_Z * sizeof(float);
    
    // Allocate memory for three scalar fields (u, v, w) and their "next step" buffers
    float *field_u_current = (float *)malloc(single_field_mem_size);
    float *field_u_next    = (float *)malloc(single_field_mem_size);
    
    if (!field_u_current || !field_u_next) {
        perror("Failed to allocate memory for fields");
        return EXIT_FAILURE;
    }
    
    FILE *output_csv_file = fopen("simulation_results.csv", "w");
    if (!output_csv_file) {
        perror("Failed to open output CSV file");
        // Free already allocated memory before exiting
        free(field_u_current); free(field_u_next);
        return EXIT_FAILURE;
    }
    // CSV Header: iter, u_center, v_center, w_center, max_dev_u
    fprintf(output_csv_file, "TimeStep,Field_U_Center,MaxDeviation_U,StepDuration\n");
    fprintf(output_csv_file, "-1,1.0,0.0,0.0\n");
    
    // Initialize scalar fields: all points to 1.0f
    for (int ix = 0; ix < GRID_SIZE_X; ix++) {
        for (int iy = 0; iy < GRID_SIZE_Y; iy++) {
            for (int iz = 0; iz < GRID_SIZE_Z; iz++) {
                int flat_index = MAP_3D_TO_1D_INDEX(ix, iy, iz);
                field_u_current[flat_index] = 1.0f;
            }
        }
    }
    
    // Introduce a perturbation at the center of field_u
    int center_x = GRID_SIZE_X / 2;
    int center_y = GRID_SIZE_Y / 2;
    int center_z = GRID_SIZE_Z / 2;
    field_u_current[MAP_3D_TO_1D_INDEX(center_x, center_y, center_z)] += PERTURBATION_MAGNITUDE;
    fprintf(output_csv_file, "0,%.6f,%.6f,0.0\n",
            field_u_current[MAP_3D_TO_1D_INDEX(GRID_SIZE_X / 2, GRID_SIZE_Y / 2, GRID_SIZE_Z / 2)],
            calculate_max_deviation_from_one(field_u_current));
    
    printf("Starting simulation...\n");
    
    // Simulação
    gettimeofday(&start, NULL);
    for (int time_step_count = 0; time_step_count < NUM_TIME_STEPS; time_step_count++) {
        gettimeofday(&step_start, NULL);
        perform_time_step(field_u_current, field_u_next);
        gettimeofday(&step_end, NULL);

        long step_seconds = step_end.tv_sec - step_start.tv_sec;
        long step_microseconds = step_end.tv_usec - step_start.tv_usec;
        if (step_microseconds < 0) {
            step_seconds--;
            step_microseconds += 1000000;
        }
        double step_elapsed = step_seconds + step_microseconds * 1e-6;

        // Inverte next e current
        float *temp_swap_pointer;
        temp_swap_pointer = field_u_current;
        field_u_current = field_u_next;
        field_u_next = temp_swap_pointer;

        // Get center values for logging
        float u_center_val = field_u_current[MAP_3D_TO_1D_INDEX(center_x, center_y, center_z)];
        
        // Calcula a máxima diferença de field_u em relação a 1.0
        float current_max_dev_u = u_center_val - 1.0f;

        fprintf(output_csv_file, "%d,%.6f,%.6f,%.6f\n",
                time_step_count+1, u_center_val, current_max_dev_u, step_elapsed);

    }

    float final_max_dev_u = calculate_max_deviation_from_one(field_u_current);
    gettimeofday(&end, NULL);
    long seconds = end.tv_sec - start.tv_sec;
    long microseconds = end.tv_usec - start.tv_usec;
    if (microseconds < 0) {
        seconds--;
        microseconds += 1000000;
    }
    double elapsed = seconds + microseconds * 1e-6;
    printf("Total simulation time: %f seconds\n", elapsed);
    printf("\nSimulation finished.\nFinal max deviation of field_u from 1.0: %.6f\n", final_max_dev_u);
    // Cleanup
    fclose(output_csv_file);
    free(field_u_current); free(field_u_next);

    printf("Results saved to simulation_results.csv\n");
    return EXIT_SUCCESS;
}