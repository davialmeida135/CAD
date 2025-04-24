#include <stdio.h>
#include <stdbool.h>
#include <sys/time.h>
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

struct Node{
    char[10] nome;
    struct Node* next;
}

int main(){
    Node head = Node()
    head.nome = "Capitao"
    printf(head.nome)
}