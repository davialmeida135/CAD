#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

// Um preenchendo
// Um somando
// Um somando em duas variáveis S1 = A[i], S2 = A[i+1]
#define N 500000000

int main(int argc, char *argv[]) {

    int *vetor = (int*)malloc(N * sizeof(int));
    if (vetor == NULL) {
        printf("Memory allocation failed\n");
        return 1;
    }
    
    struct timeval start, end;
    gettimeofday(&start, NULL);

    // Preenchendo a lista
    for (int j = 0; j < N; j++) {
        vetor[j] = j;
    }
    
    gettimeofday(&end, NULL);
    long seconds = end.tv_sec - start.tv_sec;
    long microseconds = end.tv_usec - start.tv_usec;
    double elapsed = seconds + microseconds * 1e-6;
    printf("List fill took %f seconds\n", elapsed);

    gettimeofday(&start, NULL);

    long long int soma = 0;
    for (int j = 0; j < N; j++) {
        soma = soma + vetor[j];
    }

    gettimeofday(&end, NULL);
    seconds = end.tv_sec - start.tv_sec;
    microseconds = end.tv_usec - start.tv_usec;
    elapsed = seconds + microseconds * 1e-6;
    printf("Sum in one variable took %f seconds with result %d\n", elapsed, soma);

    gettimeofday(&start, NULL);

    long long int S1 = 0;
    long long int S2 = 0;

    for (int j = 0; j < N - 1; j += 2) {
        S1 += vetor[j];
        S2 += vetor[j + 1];
    }

    gettimeofday(&end, NULL);
    seconds = end.tv_sec - start.tv_sec;
    microseconds = end.tv_usec - start.tv_usec;
    elapsed = seconds + microseconds * 1e-6;
    printf("Sum in two variables took %f seconds with result %d\n", elapsed, S1 + S2);

    free(vetor);

    return 0;
}