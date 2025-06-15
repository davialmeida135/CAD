#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Define the matrix dimensions to test (M=N=current_dim)
    int test_dimensions[] = {100,400,800,1000}; // Example sizes, ensure they are divisible by 'size' or adjust skip logic
    int num_test_dimensions = sizeof(test_dimensions) / sizeof(test_dimensions[0]);

    if (rank == 0) {
        printf("Testing on %d MPI processes.\n", size);
        printf("---------------------------------------------------------------------\n");
        printf("Matrix Dim | Comp Time (s) | Comm Time (s) | Total Time (s) (Rank 0)\n");
        printf("---------------------------------------------------------------------\n");
    }

    for (int k = 0; k < num_test_dimensions; ++k) {
        int current_dim = test_dimensions[k];
        int M_dim, N_dim; // M_dim and N_dim will be set to current_dim

        double *A_glob = NULL, *x_glob = NULL, *y_glob = NULL;
        double *local_A, *local_y;
        int rows_per_proc;
        int skip_this_iteration = 0;

        double iter_start_time, iter_end_time;
        double iter_computation_time = 0.0, iter_communication_time = 0.0, iter_total_time = 0.0;
        double t_op_start, t_op_end;

        if (rank == 0) {
            M_dim = current_dim;
            N_dim = current_dim;
            if (M_dim % size != 0) {
                fprintf(stderr, "Info: Dimension %d is not divisible by %d processes. Skipping this dimension.\n", M_dim, size);
                skip_this_iteration = 1;
            }
        }

        // Broadcast M_dim, N_dim, and the skip_this_iteration flag from rank 0
        MPI_Bcast(&M_dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&N_dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&skip_this_iteration, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (skip_this_iteration) {
            if (rank == 0) {
                printf("%-10d | SKIPPED       | SKIPPED       | SKIPPED\n", current_dim);
            }
            MPI_Barrier(MPI_COMM_WORLD); // Ensure all ranks acknowledge skip before continuing
            continue; // All processes skip to the next dimension
        }

        // If not skipping, proceed with allocations and computations
        rows_per_proc = M_dim / size;

        if (rank == 0) {
            A_glob = (double*)malloc((long)M_dim * N_dim * sizeof(double));
            x_glob = (double*)malloc((long)N_dim * sizeof(double));
            y_glob = (double*)malloc((long)M_dim * sizeof(double));

            if (!A_glob || !x_glob || !y_glob) {
                perror("Failed to allocate global matrices/vectors on rank 0");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            srand(time(NULL) + current_dim); // Vary seed for different dimensions
            for (long i = 0; i < (long)M_dim * N_dim; i++) A_glob[i] = (double)rand() / RAND_MAX;
            for (long i = 0; i < N_dim; i++) x_glob[i] = (double)rand() / RAND_MAX;
        } else {
            // Other ranks only need to allocate x_glob to receive the broadcast
            x_glob = (double*)malloc((long)N_dim * sizeof(double));
            if (!x_glob) {
                perror("Failed to allocate x_glob on non-zero rank");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        local_A = (double*)malloc((long)rows_per_proc * N_dim * sizeof(double));
        local_y = (double*)malloc((long)rows_per_proc * sizeof(double));
        if (!local_A || !local_y) {
            perror("Failed to allocate local matrices/vectors");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        MPI_Barrier(MPI_COMM_WORLD); // Synchronize before starting the timer for the iteration
        iter_start_time = MPI_Wtime();

        // Broadcast vector x from rank 0 to all processes
        t_op_start = MPI_Wtime();
        MPI_Bcast(x_glob, N_dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        t_op_end = MPI_Wtime();
        iter_communication_time += t_op_end - t_op_start;

        // Scatter matrix A from rank 0 to all processes
        t_op_start = MPI_Wtime();
        MPI_Scatter(A_glob, rows_per_proc * N_dim, MPI_DOUBLE, local_A, rows_per_proc * N_dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        t_op_end = MPI_Wtime();
        iter_communication_time += t_op_end - t_op_start;
        
        MPI_Barrier(MPI_COMM_WORLD); // Sync before computation timing
        t_op_start = MPI_Wtime();
        // Local computation: local_y = local_A * x_glob
        for (int i = 0; i < rows_per_proc; i++) {
            local_y[i] = 0.0;
            for (int j = 0; j < N_dim; j++) {
                local_y[i] += local_A[i * N_dim + j] * x_glob[j];
            }
        }
        MPI_Barrier(MPI_COMM_WORLD); // Sync after computation timing
        t_op_end = MPI_Wtime();
        iter_computation_time = t_op_end - t_op_start;

        // Gather results from all processes to y_glob on rank 0
        t_op_start = MPI_Wtime();
        MPI_Gather(local_y, rows_per_proc, MPI_DOUBLE, y_glob, rows_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        t_op_end = MPI_Wtime();
        iter_communication_time += t_op_end - t_op_start;

        iter_end_time = MPI_Wtime();
        iter_total_time = iter_end_time - iter_start_time;

        if (rank == 0) {
            printf("%-10d | %-13.6f | %-13.6f | %-17.6f\n",
                   current_dim, iter_computation_time, iter_communication_time, iter_total_time);
            
            // Optional: verification of results on rank 0
            // double y0_check = 0.0;
            // if (M_dim > 0 && N_dim > 0) { // Basic check for y_glob[0]
            //    for(int j=0; j<N_dim; j++) y0_check += A_glob[j] * x_glob[j]; // A_glob[0*N_dim + j]
            //    printf("Debug: y_glob[0] (calculated) = %f, y0_check = %f\n", y_glob[0], y0_check);
            // }
        }

        // Cleanup for the current iteration
        if (rank == 0) {
            free(A_glob);
            free(y_glob);
        }
        free(x_glob); // All ranks allocated and must free x_glob
        free(local_A);
        free(local_y);
        
        MPI_Barrier(MPI_COMM_WORLD); // Ensure all ranks finished cleanup before next iteration
    }

    if (rank == 0) {
        printf("---------------------------------------------------------------------\n");
    }

    MPI_Finalize();
    return 0;
}