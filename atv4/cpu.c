#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
#include <math.h> // Include for math functions

// Um preenchendo
// Um somando
// Um somando em duas variáveis S1 = A[i], S2 = A[i+1]

// pragma omp parallel for num_threads()
// gcc -fopenmp cpu.c -o cpu -lm  // Add -lm for linking math library
int main(int argc, char *argv[]) {
    struct timeval start, end;

    long int limit = 23000000; 
    double sum = 0.0; 
    int num_threads = 0;
    printf("Enter the value of N: ");
    scanf("%d", &num_threads);

    printf("Calculating sum with %d threads...\n", num_threads);

    gettimeofday(&start, NULL);

    #pragma omp parallel for num_threads(num_threads)
    for(long int k = 1.0; k <= limit; k++) { 
        double term1 = sin(k / 100.0) * cos(k / 50.0);
        double term2 = sqrt(fabs(tan((double)k))) + 1.0;
        double term3 = log(term2);
        double term4 = exp(-k / (double)limit);
        sum += pow(term1, 4) * term3 * term4 + cbrt((double)k); 
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
    printf("Sum = %f\n", sum); // Print the calculated sum
    printf("Calculation took %f seconds\n", elapsed);

    return 0; // Added return statement
}