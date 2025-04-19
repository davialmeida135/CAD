#include <stdio.h>
#include <stdbool.h>
#include <sys/time.h>
#include <omp.h>
#include <math.h>

bool is_prime(int num) {
    if (num < 2) return false;
    for (int i = 2; i * i <= num; i++) {
        if (num % i == 0)
            return false;
    }
    return true;
}

int main() {
    struct timeval start, end;
    int N, count = 0;
    printf("Enter the value of N: ");
    scanf("%d", &N);

    gettimeofday(&start, NULL);
    #pragma omp parallel for num_threads(N) reduction(+:soma)
    for (int i = 2; i <= N; i++) {
        if (is_prime(i)) {
            count++;
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
    printf("Number of primes between 2 and %d: %d\n", N, count);
    printf("Calculation took %f seconds\n", elapsed);
    return 0;
}
