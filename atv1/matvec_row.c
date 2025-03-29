#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

// Function to get wall time (higher resolution than clock())
double get_wall_time() {
    struct timeval time;
    // Use gettimeofday, check for errors
    if (gettimeofday(&time, NULL)) {
        perror("gettimeofday failed");
        return 0;
    }
    // Return time in seconds as a double
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

// Função de multiplicação de matriz-vetor: result = A * x
void mat_vec_mul_int(int n, const int *A, const int *x, int *result_vector) {
    for (int i = 0; i < n; ++i) {
        int sum = 0;
        for (int j = 0; j < n; ++j) {
            // Row-major
            result_vector[i] += A[i * n + j] * x[j];
        }
    }
}


void test_matrix_size(int n) {

    long n_squared = (long)n * n;

    printf("===================================================================\n");
    printf("Testing n = %d (Matrix %ld elements, Vector %d elements)\n", n, n_squared, n);

    int *matriz = (int *)malloc(n_squared * sizeof(int));
    int *vetor = (int *)malloc(n * sizeof(int));
    int *result_vector = (int *)malloc(n * sizeof(int));

    printf("  Preenchendo a matriz\n");
    for (int i = 0; i < n; ++i) {
        vetor[i] = (i % 10) + 1; 
        for (int j = 0; j < n; ++j) {
            matriz[i * n + j] = (i - j + n) % 100; // Pra não ficar numeros grandes
        }
    }

    // Warm-up
    mat_vec_mul_int(n, matriz, vetor, result_vector);

    // Quantidade de repetições por tamanho
    int repetitions = 50; 
    if (n > 200) repetitions = 20;
    if (n > 500) repetitions = 10;
    if (n > 1000) repetitions = 5;
    if (n > 4000) repetitions = 2;  
    if (n > 10000) repetitions = 1;
    printf("  Executando %d multiplicações de matriz\n", repetitions);

    double start_time = get_wall_time();
    for (int r = 0; r < repetitions; ++r) {
        mat_vec_mul_int(n, matriz, vetor, result_vector);
    }
    double end_time = get_wall_time();

    double total_elapsed_time = end_time - start_time;
    double average_time_s = total_elapsed_time / repetitions; // Tempo médio por execução

    printf("  Repetições: %d\n", repetitions);
    printf("  Soma dos tempos de %d runs: %.6f segundos\n", repetitions, total_elapsed_time);
    printf("  Tempo médio por cálculo: %.9f segundos\n", average_time_s);

    free(matriz);
    free(vetor);
    free(result_vector);

    printf("===================================================================\n");
}


int main() {

    int sizes[] = {
        // Below L1
        32, 64, 80,
        // Around L1/L2
        88, 96, 112, 128, 192, 224,
        // Around L2/L3
        240, 256, 272, 384, 512, 768, 1024, 1280,
        // Around L3/RAM
        1400, 1448, 1536, 1792, 2048,
        // Larger (RAM bound - keep a few if desired)
        3000, 4000 // , 5000 // , 10000 // These get very slow
    };
    int num_sizes = sizeof(sizes) / sizeof(sizes[0]);

    for (int i = 0; i < num_sizes; i++) {
        test_matrix_size(sizes[i]);
    }

    return 0;
}