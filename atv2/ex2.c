#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
// Um preenchendo
// Um somando
// Um somando em duas variáveis S1 = A[i], S2 = A[i+1]
#define N 10000000
int main(int argc, char *argv[]) {
    //int n = argv[1] ? atoi(argv[1]) : 1000000;
    // Preenchendo a lista
    
    struct timeval start, end;
    gettimeofday(&start, NULL);

    int i[N];
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
    for (int j = 0; j < N; j++) {
        soma = soma + i[j];
    }

    gettimeofday(&end, NULL);
    seconds = end.tv_sec - start.tv_sec;
    microseconds = end.tv_usec - start.tv_usec;
    elapsed = seconds + microseconds * 1e-6;
    printf("Sum in one variable took %f seconds with result %d\n", elapsed, soma);

    gettimeofday(&start, NULL);

    int S1 = 0;
    int S2 = 0;
    for (int j = 0; j < N; j += 2) {
        S1 += i[j];
        S2 += i[j + 1];
    }

    gettimeofday(&end, NULL);
    seconds = end.tv_sec - start.tv_sec;
    microseconds = end.tv_usec - start.tv_usec;
    elapsed = seconds + microseconds * 1e-6;
    printf("Sum in two variables took %f seconds with result %d\n", elapsed, S1 + S2);


    return 0;
}