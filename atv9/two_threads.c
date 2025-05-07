#include <stdio.h>
#include <stdbool.h>
#include <sys/time.h>
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

typedef struct Node{
    struct Node* next;
    double valor;
} Node;

#define N_THREADS 2
#define N_ITERATIONS 30

void print_list(Node* head) {
    printf("Lista "); // Imprime como Lista 1 e Lista 2
    Node* current = head;
    int count = 0;
    while (current != NULL) {
        printf("%.2f -> ", current->valor);
        current = current->next;
        count++;
        if (count > N_ITERATIONS + 5) { // Prevenção de loop infinito em caso de erro
            printf("... (lista muito longa ou ciclo detectado)\n");
            break;
        }
    }
    printf("NULL\n");
    printf("Tamanho da lista: %d\n", count-1);
}

void free_list(Node* head) {
    Node* current = head;
    Node* next_node;
    while (current != NULL) {
        next_node = current->next;
        free(current);
        current = next_node;
    }
}

int main(){
    Node* lista1 = (Node*)malloc(sizeof(Node));
    lista1->valor = 1;

    Node* lista2 = (Node*)malloc(sizeof(Node));
    lista2->valor = 1;
    unsigned int seed = 42;

    #pragma omp parallel num_threads(2)
    {
        #pragma omp single
        {
            for(int i=0; i<N_ITERATIONS; i++){
                #pragma omp task
                {
                    Node* next_node = (Node*)malloc(sizeof(Node));
                    double r = (double)rand_r(&seed) / RAND_MAX;
                    next_node->valor = r;
            
                    if(r>0.5){
                        # pragma omp critical (l1)
                        {
                        next_node->next = lista1;
                        lista1 = next_node;
                        }
                    }else{
                        # pragma omp critical (l2)
                        {
                        next_node->next = lista2;
                        lista2 = next_node; 
                        }  
                    }
                }
            }  
        }
    }
    print_list(lista1);
    print_list(lista2);
}