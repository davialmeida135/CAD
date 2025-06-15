#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int test_dimensions[] ={100,400,800,1000};
    int num_test_dimensions = sizeof(test_dimensions) / sizeof(test_dimensions[0]);

    if (rank == 0) {
        printf("Testing on %d MPI processes (column-wise distribution).\n", size);
        printf("---------------------------------------------------------------------\n");
        printf("Matrix Dim | Comp Time (s) | Comm Time (s) | Total Time (s) (Rank 0)\n");
        printf("---------------------------------------------------------------------\n");
    }

    for (int k = 0; k < num_test_dimensions; ++k) {
        int current_dim = test_dimensions[k];
        int M_dim = current_dim, N_dim = current_dim;
        int cols_per_proc = N_dim / size;
        int skip = (N_dim % size != 0);

        double *A_glob = NULL, *x_glob = NULL, *y_glob = NULL;
        double *local_A = NULL, *local_x = NULL, *local_y = NULL;

        double t_start, t_end;
        double comm_time = 0.0, comp_time = 0.0, total_time = 0.0;

        MPI_Bcast(&M_dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&N_dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&skip, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (skip) {
        if (rank == 0)
            printf("%-10d | SKIPPED    | SKIPPED    | SKIPPED\n", current_dim);
            continue;
        }

        if (rank == 0) {
            A_glob = (double*)malloc((long)M_dim * N_dim * sizeof(double));
            x_glob = (double*)malloc((long)N_dim * sizeof(double));
            y_glob = (double*)calloc(M_dim, sizeof(double));

            srand(time(NULL) + current_dim);
            for (long i = 0; i < (long)M_dim * N_dim; i++) A_glob[i] = (double)rand() / RAND_MAX;
            for (long i = 0; i < N_dim; i++) x_glob[i] = (double)rand() / RAND_MAX;
        } else {
            x_glob = (double*)malloc((long)N_dim * sizeof(double));
        }

        local_A = (double*)malloc((long)M_dim * cols_per_proc * sizeof(double));
        local_x = (double*)malloc((long)cols_per_proc * sizeof(double));
        local_y = (double*)calloc(M_dim, sizeof(double));

        MPI_Barrier(MPI_COMM_WORLD);
        t_start = MPI_Wtime();

        double t_comm_start = MPI_Wtime();
        //MPI_Bcast(x_glob, N_dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        double t_comm_end = MPI_Wtime();
        comm_time += t_comm_end - t_comm_start;

        // Define derived type to scatter columns
        MPI_Datatype column_type, resized_column_type;
        MPI_Type_vector(M_dim, cols_per_proc, N_dim, MPI_DOUBLE, &column_type);
        MPI_Type_create_resized(column_type, 0, cols_per_proc * sizeof(double), &resized_column_type);
        MPI_Type_commit(&resized_column_type);

        t_comm_start = MPI_Wtime();
        MPI_Scatter(A_glob, 1, resized_column_type,
            local_A, M_dim * cols_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(x_glob, cols_per_proc, MPI_DOUBLE,
            local_x, cols_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        t_comm_end = MPI_Wtime();
        comm_time += t_comm_end - t_comm_start;

        MPI_Barrier(MPI_COMM_WORLD);
        double t_comp_start = MPI_Wtime();

        // Cada processo calcula y parcial
        for (int i = 0; i < M_dim; i++) {
            for (int j = 0; j < cols_per_proc; j++) {
                local_y[i] += local_A[i * cols_per_proc + j] * local_x[j];
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        double t_comp_end = MPI_Wtime();
        comp_time = t_comp_end - t_comp_start;

        t_comm_start = MPI_Wtime();
        MPI_Reduce(local_y, y_glob, M_dim, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        t_comm_end = MPI_Wtime();
        comm_time += t_comm_end - t_comm_start;

        t_end = MPI_Wtime();
        total_time = t_end - t_start;

        if (rank == 0) {
        printf("%-10d | %-13.6f | %-13.6f | %-17.6f\n",
            current_dim, comp_time, comm_time, total_time);
        }

        if (rank == 0) {
        free(A_glob);
        free(y_glob);
        }
        free(x_glob);
        free(local_A);
        free(local_x);
        free(local_y);

        MPI_Type_free(&resized_column_type);
        MPI_Type_free(&column_type);

        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (rank == 0) {
        printf("---------------------------------------------------------------------\n");
    }

    MPI_Finalize();
    return 0;
}
