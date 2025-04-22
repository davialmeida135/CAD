#include <stdio.h>
#include <stdbool.h>
#include <sys/time.h>
#include <omp.h>
#include <math.h>
    gettimeofday(&start, NULL);
    #pragma omp parallel for num_threads(N) //reduction(+:soma)
    for (int i = 2; i <= N; i++) {
        
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
    printf("Number of primes between 2 and %d: %d\n", N, count);
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int main() {
    int total_pontos = 1000000; // número de pontos a serem gerados
    int dentro_circulo = 0;
    double x, y;
    
    // Semente para gerar números aleatórios
    srand(time(NULL));

    for (int i = 0; i < total_pontos; i++) {
        // Gera dois números aleatórios entre 0 e 1
        x = (double)rand() / RAND_MAX;
        y = (double)rand() / RAND_MAX;

        // Verifica se o ponto (x, y) está dentro do círculo
        if ((x * x + y * y) <= 1.0) {
            dentro_circulo++;
        }
    }

    // Calcula a aproximação de Pi
    double pi = 4.0 * dentro_circulo / total_pontos;
    printf("Aproximacao de Pi com %d pontos: %f\n", total_pontos, pi);

    return 0;
}
