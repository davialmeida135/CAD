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
    char nome[10];
} Node;

int main(){
    Node capitao;
    strcpy(capitao.nome, "Capitaoo\n");

    Node recruta;
    strcpy(recruta.nome, "Recrutaa\n");
    capitao.next = &recruta;

    Node kowalski;
    strcpy(kowalski.nome, "Kowalski\n");
    recruta.next = &kowalski;
    Node rico;
    strcpy(rico.nome, "Rico\n");
    kowalski.next = &rico;

    Node alex;
    strcpy(alex.nome, "Alex\n");
    rico.next = &alex;

    Node gloria;
    strcpy(gloria.nome, "Gloria\n");
    alex.next = &gloria;

    Node melman;
    strcpy(melman.nome, "Melman\n");
    gloria.next = &melman;

    Node marty;
    strcpy(marty.nome, "Marty\n");
    melman.next = &marty;

    Node julian;
    strcpy(julian.nome, "Julian\n");
    marty.next = &julian;

    Node mort;
    strcpy(mort.nome, "Mort\n");
    julian.next = &mort;
    mort.next = NULL;

    printf("%s", capitao.next->nome);
    // ordem de impressão [capitao, recruta, kowalski, rico, alex, gloria, melman, marty, julian, mort]
}