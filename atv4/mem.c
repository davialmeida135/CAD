#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

// Um preenchendo
// Um somando
// Um somando em duas variáveis S1 = A[i], S2 = A[i+1]

// Compile with: gcc -fopenmp mem.c -o mem
#define N 1000000000 // Increased N for a more significant workload

int main(int argc, char *argv[]) {
    // Allocate array on the heap instead of the stack
    int *A = (int *)malloc(N * sizeof(int));
    if (A == NULL) {
        perror("Failed to allocate memory for array i");
        return 1; // Indicate failure
    }
    int *B = (int *)malloc(N * sizeof(int));
    if (B == NULL) {
        perror("Failed to allocate memory for array j");
        return 1; // Indicate failure
    }
    

    // Optional: Get number of threads from command line or set default
    int num_threads = 1; // Default number of threads
    printf("Using %d threads.\n", num_threads);

    // --- Preenchendo a lista ---
    struct timeval start, end;
    gettimeofday(&start, NULL);

    // Corrected pragma: parallel spelling
    #pragma omp parallel for
    for (int j = 0; j < N; j++) {
        if (j %2 == 0) {
            A[j] = j;
            B[j] = -j;
        } else {
            A[j] = -j; 
            B[j] = j;
        }
    }

    gettimeofday(&end, NULL);
    long seconds = end.tv_sec - start.tv_sec;
    long microseconds = end.tv_usec - start.tv_usec;
    // Handle potential negative microseconds
    if (microseconds < 0) {
        seconds--;
        microseconds += 1000000;
    }
    double elapsed = seconds + microseconds * 1e-6;
    printf("List fill took %f seconds\n", elapsed);

    // Sum
    gettimeofday(&start, NULL);

    long long int soma = 0;
    #pragma omp parallel for num_threads(num_threads)
    for (long int j = 0; j < N; j++) {
        soma = soma + A[j] + B[j];
    }

    gettimeofday(&end, NULL);
    seconds = end.tv_sec - start.tv_sec;
    microseconds = end.tv_usec - start.tv_usec;
    // Handle potential negative microseconds
    if (microseconds < 0) {
        seconds--;
        microseconds += 1000000;
    }
    elapsed = seconds + microseconds * 1e-6;

    printf("Sum in one variable took %f seconds with result %lld\n", elapsed, soma);

    // Free the allocated memory
    free(A);
    free(B);

    return 0; // Indicate success
}
