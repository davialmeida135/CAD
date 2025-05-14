#include <stdio.h>
#include <stdbool.h>
#include <sys/time.h>
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

int main() {
    long int total_pontos = 20000000;
    
    struct timeval start, end;
    // x, y don't need to be defined here if declared inside parallel region

    gettimeofday(&start, NULL);
    int max_threads = omp_get_max_threads();
    printf("Maximum threads: %d\n", max_threads);


    double x, y;
    srand(time(NULL));
    int total_dentro_circulo = 0;
    #pragma omp parallel default(none) shared(total_pontos, total_dentro_circulo) private(x, y)
    {
        unsigned int seed = time(NULL) + omp_get_thread_num();
        int thread_num = omp_get_thread_num();
        #pragma omp for reduction(+:total_dentro_circulo)
        for (int i = 0; i < total_pontos; i++) {
            x = (double)rand() / RAND_MAX;
            y = (double)rand() / RAND_MAX; 

            if (((x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5)) <= 0.25) {
                {
                    total_dentro_circulo++;
                }
            }
        }     
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

    double pi = 4.0 * total_dentro_circulo / total_pontos;
    printf("Aproximacao de Pi com %ld pontos: %f\n", total_pontos, pi);

    double pi_real = M_PI;
    printf("Erro : %f\n", fabs(pi_real - pi));
    return 0;
}