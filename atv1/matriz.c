#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void calcular_matriz(int n){
    // Input do tamanho da matriz
    //int n = 3500;
    
    // Construção da matriz 
    int *matriz = (int *)malloc(n*n * sizeof(int *));
    // Dynamic allocation for vectors
    int *vetor = (int *)malloc(n * sizeof(int));
    int *nova_matriz = (int *)malloc(n * sizeof(int));
    int c=1;
    for(int i=0;i<n;i++){
        vetor[i] = i%10;
        for(int j=0;j<n;j++){
            matriz[i*n+j] = i-j;
            c++;
            //printf("%d ", matriz[i][j]);
        }
        //printf("\n");
    }
    //printf("===================================================================");
    int soma = 0;
    struct timespec begin, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &begin);
    for(int i=0;i<n;i++){
        soma = 0;
        for(int j=0;j<n;j++){
            soma += matriz[i*n+j] * vetor[j];
            //printf("i = %d, j = %d\n", i, j);
            //printf("%d * %d = %d\n", matriz[i][j], vetor[j], matriz[i][j] * vetor[j]);
        }
        nova_matriz[i] = soma;
        soma=0;
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("===================================================================");
    printf("\nMultiplicação de matriz %d por %d por vetor de tamanho %d\n", n, n, n);
    printf ("Total time = %f seconds\n",
            (end.tv_nsec - begin.tv_nsec) / 1000000000.0 +
            (end.tv_sec  - begin.tv_sec));
    printf("\n");
}
int main(){
    int sizes[] = {100, 200, 350, 1000, 3000,4000, 5000,6000,7000};
    int n;
    for (int i = 0; i < sizeof(sizes)/sizeof(sizes[0]); i++) {
        n = sizes[i];
        calcular_matriz(n);
        printf("===================================================================\n");
    }
}