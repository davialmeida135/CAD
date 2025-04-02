#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>  // Add this header for gettimeofday

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

void calcular_matriz(int n){
    // Input do tamanho da matriz
    struct timeval t0, t1;
    
    // Construção da matriz 
    int *matriz = (int *)malloc(n*n * sizeof(int));
    // Dynamic allocation for vectors
    int *vetor = (int *)malloc(n * sizeof(int));
    int *nova_matriz = (int *)malloc(n * sizeof(int));
    int c=1;
    for(int i=0;i<n;i++){
        vetor[i] = i%10 +1;
        for(int j=0;j<n;j++){
            matriz[i*n+j] = i-j;
            c++;
        }
    }
    // Warmup
    int soma = 0;
    for(int i=0;i<n;i++){
        soma = 0;
        for(int j=0;j<n;j++){
            soma += matriz[i*n+j] * vetor[j];
        }
        nova_matriz[i] = soma;
    }    
    soma = 0;
    gettimeofday(&t0, 0);
    for(int i=0;i<n;i++){
        soma = 0;
        for(int j=0;j<n;j++){
            soma += matriz[i*n+j] * vetor[i];
        }
        nova_matriz[i] = soma;
    }
    gettimeofday(&t1, 0);
    
    printf("===================================================================");
    printf("\nMultiplicação de matriz %d por %d por vetor de tamanho %d\n", n, n, n);
    long elapsed = ((t1.tv_sec * 1000000) + t1.tv_usec) - ((t0.tv_sec * 1000000) + t0.tv_usec);
    float time_in_seconds = elapsed / 1000000.0;  // Convert to seconds
    float n_float = (float)n;
    float time_per_iteration = time_in_seconds / (n_float * n_float/1000000.0);
    printf("Total time = %.9f seconds\n%.9f seconds per 1000000 iterations\n", time_in_seconds, time_per_iteration);
    printf("\n");
    
    // Free allocated memory
    free(matriz);
    free(vetor);
    free(nova_matriz);
}
int main(){
    int sizes[] = {80,89,90,91,100, 250,255,256,257,260,1440,1447,1448,1449,1460, 5000,10000,15000,20000,30000,35000};
    float n;
    for (int i = 0; i < sizeof(sizes)/sizeof(sizes[0]); i++) {
        n = sizes[i];
        calcular_matriz(n);
        printf("===================================================================\n");
    }
}