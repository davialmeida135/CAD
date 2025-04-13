#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
// Um preenchendo
// Um somando
// Um somando em duas variáveis S1 = A[i], S2 = A[i+1]

// pragma omp paralel for()
// gcc -f openmp mem.c -o mem
#define N 10000000
int main(int argc, char *argv[]) {
    //int n = argv[1] ? atoi(argv[1]) : 1000000;
    // Preenchendo a lista
    
    struct timeval start, end;
    gettimeofday(&start, NULL);

    int i[N];
    #pragma omp paralel for num_threads()
    for (int j = 0; j < N; j++) {
        i[j] = j;
    }
    
    gettimeofday(&end, NULL);
    long seconds = end.tv_sec - start.tv_sec;
    long microseconds = end.tv_usec - start.tv_usec;
    double elapsed = seconds + microseconds * 1e-6;
    printf("List fill took %f seconds\n", elapsed);

    gettimeofday(&start, NULL);

    int soma = 0;
    #pragma omp paralel for numthreads()
    for (int j = 0; j < N; j++) {
        soma = soma + i[j];
    }

    gettimeofday(&end, NULL);
    seconds = end.tv_sec - start.tv_sec;
    microseconds = end.tv_usec - start.tv_usec;
    elapsed = seconds + microseconds * 1e-6;
    printf("Sum in one variable took %f seconds with result %d\n", elapsed, soma);
    return 0;
}