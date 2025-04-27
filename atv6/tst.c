#include <stdio.h>
#include <stdbool.h>
#include <sys/time.h>
#include <omp.h>
#include <math.h>
// #include <stdio.h> // Duplicate include removed
#include <stdlib.h> // Needed for calloc/free
#include <time.h>

int main() {
    int total_pontos = 10000;
    int dentro_circulo_total = 0; // Final result accumulator
    struct timeval start, end;

    gettimeofday(&start, NULL);
    int max_threads = omp_get_max_threads();

    printf("Maximum threads: %d\n", max_threads);

    // Lista que vai armazenar os contadores parciais de cada thread
    int *partial_counts = (int *)calloc(max_threads, sizeof(int));
    if (partial_counts == NULL) {
        perror("Failed to allocate memory for partial counts");
        return 1;
    }
    double x, y;
    int tid = 0; // Thread ID
    int seed; 
    int local_count = 0;
    #pragma omp parallel default(none) shared(total_pontos, partial_counts) private(x, y, tid, seed) firstprivate(local_count)
    {
        // Variáveis privadas para cada thread	
        
        tid = omp_get_thread_num();
        seed = time(NULL) + tid;

        #pragma omp for nowait // Threads podem sair do loop quando terminarem
        for (int i = 0; i < total_pontos; i++) {

            x = (double)rand_r(&seed) / RAND_MAX;
            y = (double)rand_r(&seed) / RAND_MAX;

            if (((x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5)) <= 0.25) {
                local_count++; 
            }
        }

        partial_counts[tid] = local_count;

    }

    gettimeofday(&end, NULL);

    
    for (int i = 0; i < max_threads; i++) {
        dentro_circulo_total += partial_counts[i];
    }

    free(partial_counts);

    long seconds = end.tv_sec - start.tv_sec;
    long microseconds = end.tv_usec - start.tv_usec;

    if (microseconds < 0) {
        seconds--;
        microseconds += 1000000;
    }

    double elapsed = seconds + microseconds * 1e-6;
    printf("Calculation took %f seconds\n", elapsed);

    // Use the final summed count
    double pi = 4.0 * dentro_circulo_total / total_pontos;
    printf("Aproximacao de Pi com %d pontos: %f\n", total_pontos, pi);

    double pi_real = M_PI;
    printf("Valor real de Pi: %f\n", pi_real);
    printf("Erro : %f\n", fabs(pi_real - pi));
    return 0;
}