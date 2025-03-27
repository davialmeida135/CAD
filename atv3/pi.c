#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>

int main(int argc, char *argv[])
{
    double pi = 0.0;
    int iterations = 100000000;
    double term = 0.0;
    float denom=1.0;
    
    // For measuring execution time
    struct timespec start_time, end_time;
    double execution_time;
    
    // CSV
    FILE *fp;
    fp = fopen("pi_results.csv", "w");
    if (fp == NULL) {
        printf("Error opening file!\n");
        return 1;
    }
    
    // CSV Colunas
    fprintf(fp, "Iteration,PI_Approximation,Error,Time(s)\n");
    
    // Get start time
    clock_gettime(CLOCK_MONOTONIC, &start_time);
    
    int log_interval = 100000; // Log every 10,000 iterations
    struct timespec current_time;
    double current_execution_time;

    for (int i = 1; i < iterations; i++) {
        if (i % 2 == 0) {
            term = -4.0/denom;
        } else {
            term = 4.0/denom;
        }
        pi += term;
        denom += 2;
        
        // Log at regular intervals or for the first few iterations
        if (i <= 100 || i % log_interval == 0 || i == iterations-1) {
            clock_gettime(CLOCK_MONOTONIC, &current_time);
            current_execution_time = (current_time.tv_sec - start_time.tv_sec) + 
                                    (current_time.tv_nsec - start_time.tv_nsec) / 1000000000.0;
            
            fprintf(fp, "%d,%.15f,%.15f,%.9f\n", i, pi, fabs(pi - M_PI), current_execution_time);
        }
    }
    
    // Get end time
    clock_gettime(CLOCK_MONOTONIC, &end_time);
    
    // Calculate execution time in seconds
    execution_time = (end_time.tv_sec - start_time.tv_sec) + 
                    (end_time.tv_nsec - start_time.tv_nsec) / 1000000000.0;

    printf("Approximation of pi using %d iterations: %.15f\n", iterations, pi);
    printf("Real value of pi (M_PI): %.15f\n", M_PI);
    printf("Absolute error: %.15f\n", fabs(pi - M_PI));
    printf("Relative error: %.15f%%\n", 100 * fabs(pi - M_PI) / M_PI);
    printf("Execution time: %.9f seconds\n", execution_time);
    
    fclose(fp);
    printf("Results saved to pi_results.csv\n");
    
    return 0;
}