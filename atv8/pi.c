#include <stdio.h>
#include <stdbool.h>
#include <sys/time.h>
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

int main() {
    long int total_pontos = 10000000000;
    
    struct timeval start, end;
    // x, y don't need to be defined here if declared inside parallel region

    gettimeofday(&start, NULL);
    int max_threads = omp_get_max_threads();
    printf("Maximum threads: %d\n", max_threads);

    int *partial_counts = (int *)calloc(max_threads, sizeof(int));
    if (partial_counts == NULL) {
        perror("Failed to allocate memory for partial counts");
        return 1;
    }
    double x, y;
    //srand(time(NULL));
    #pragma omp parallel default(none) shared(total_pontos, partial_counts) private(x, y)
    {
        int thread_num = omp_get_thread_num();
        unsigned int seed = time(NULL) + thread_num;
        #pragma omp for
        for (int i = 0; i < total_pontos; i++) {

            x = (double)rand_r(&seed) / RAND_MAX;
            y = (double)rand_r(&seed) / RAND_MAX; 
            //printf("Thread %d: x = %f, y = %f\n", thread_num, x, y);

            if (((x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5)) <= 0.25) {
                {
                    partial_counts[thread_num]++;
                }
            }
        }
    }

    int total_count = 0;
    for (int i = 0; i < max_threads; i++) {
        // Sum up the counts from each thread
        printf("Thread %d: %d\n", i, partial_counts[i]);
        total_count += partial_counts[i];
    }
    gettimeofday(&end, NULL);

    // ... rest of the calculation and printing ...
    long seconds = end.tv_sec - start.tv_sec;
    long microseconds = end.tv_usec - start.tv_usec;

    if (microseconds < 0) {
        seconds--;
        microseconds += 1000000;
    }

    double elapsed = seconds + microseconds * 1e-6;
    printf("Calculation took %f seconds\n", elapsed);

    double pi = 4.0 * total_count / total_pontos;
    printf("Aproximacao de Pi com %ld pontos: %f\n", total_pontos, pi);

    double pi_real = M_PI;
    printf("Erro : %f\n", fabs(pi_real - pi));
    return 0;
}