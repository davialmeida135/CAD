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

int main(){
    Node lista1;
    lista1.valor = 1;
    Node lista2;
    lista2.valor = 1;
    unsigned int seed = 42;
    #pragma omp parallel num_threads(2)
    {
        Node next_node;
        next_node.valor = r;
        r = rand_r(&seed)
        if(r>0.5){
            lista1.next = &next_node;
        }
        else{
            lista2.next = &next_node;
        }
    }
}