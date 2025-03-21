#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void main(){
    // Input do tamanho da matriz
    int n = 1000;
    
    // Construção da matriz 
    int matriz[n][n];
    int vetor[n];
    int c=1;
    for(int i=0;i<n;i++){
        vetor[i] = i%10;
        for(int j=0;j<n;j++){
            matriz[i][j] = i-j;
            c++;
            //printf("%d ", matriz[i][j]);
        }
        //printf("\n");
    }
    //printf("===================================================================");

    int nova_matriz[n];
    int soma = 0;
    struct timespec begin, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &begin);
    for(int i=0;i<n;i++){
        soma = 0;
        for(int j=0;j<n;j++){
            soma += matriz[i][j] * vetor[j];
            //printf("i = %d, j = %d\n", i, j);
            //printf("%d * %d = %d\n", matriz[i][j], vetor[j], matriz[i][j] * vetor[j]);
        }
        nova_matriz[i] = soma;
        soma=0;
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    //for(int i=0;i<n;i++){
    //    printf("%d ", vetor[i]);
    //}
    
    //printf("===================================================================");
    
    //for(int i=0;i<n;i++){
    //    printf("%d ", nova_matriz[i]);
    //}
    printf("===================================================================");
    printf("\nMultiplicação de matriz %d por %d por vetor de tamanho %d\n", n, n, n);
    printf ("Total time = %f seconds\n",
            (end.tv_nsec - begin.tv_nsec) / 1000000000.0 +
            (end.tv_sec  - begin.tv_sec));
    printf("\n");
}