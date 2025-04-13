#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
// Um preenchendo
// Um somando
// Um somando em duas variáveis S1 = A[i], S2 = A[i+1]

// pragma omp paralel for()
// gcc -f openmp mem.c -o mem
int main(int argc, char *argv[]) {
    struct timeval start, end;
    gettimeofday(&start, NULL);

    int i[N];
    long int sum = 0;
    int i = 1000000;
    #pragma omp paralel for num_threads()
    for(int k=0;k<i;k++){
        sum += sin(i)+cos(i)*sqrt(i)+tan(sum)
    }

    gettimeofday(&end, NULL);
    long seconds = end.tv_sec - start.tv_sec;
    long microseconds = end.tv_usec - start.tv_usec;
    double elapsed = seconds + microseconds * 1e-6;
    printf("Sum took %f seconds\n", elapsed);

}
