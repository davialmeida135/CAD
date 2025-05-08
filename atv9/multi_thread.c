#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

typedef struct Node {
    struct Node* next;
    double valor;
} Node;

#define N_ITERATIONS 30 // Total number of nodes to add across all lists

// Function to print a linked list
void print_list(Node* head, int list_num, int total_iterations) {
    printf("Lista %d: ", list_num);
    Node* current = head;
    int count = 0;
    while (current != NULL) {
        printf("%.2f -> ", current->valor);
        current = current->next;
        count++;
        // Safety break for very long lists or cycles
        if (count > total_iterations + 5) {
            printf("... (lista muito longa ou ciclo detectado)\n");
            break;
        }
    }
    printf("NULL\n");
    printf("Tamanho da lista %d: %d\n", list_num, count);
}

// Function to free a linked list
void free_list(Node* head) {
    Node* current = head;
    Node* next_node;
    while (current != NULL) {
        next_node = current->next;
        free(current);
        current = next_node;
    }
}

int main() {
    int num_lists;
    printf("Enter the number of lists: ");
    if (scanf("%d", &num_lists) != 1 || num_lists <= 0) {
        fprintf(stderr, "Invalid number of lists. Please enter a positive integer.\n");
        return 1;
    }

    // Allocate an array for list heads
    Node** list_heads = (Node**)malloc(num_lists * sizeof(Node*));
    if (list_heads == NULL) {
        perror("Failed to allocate memory for list heads");
        return 1;
    }
    for (int i = 0; i < num_lists; i++) {
        list_heads[i] = NULL; // Initialize all lists as empty
    }

    // Allocate and initialize an array of OpenMP locks, one for each list
    omp_lock_t* locks = (omp_lock_t*)malloc(num_lists * sizeof(omp_lock_t));
    if (locks == NULL) {
        perror("Failed to allocate memory for locks");
        free(list_heads);
        return 1;
    }
    for (int i = 0; i < num_lists; i++) {
        omp_init_lock(&locks[i]);
    }

    unsigned int seed = (unsigned int)time(NULL);

    #pragma omp parallel
    {
        #pragma omp single
        {
            for (int i = 0; i < N_ITERATIONS; i++) {
                #pragma omp task
                {
                    Node* new_node = (Node*)malloc(sizeof(Node));
                    if (new_node != NULL) {
                        // Numero aleatorio entre 0 e 1
                        new_node->valor = (double)rand_r(&seed) / RAND_MAX;
                        new_node->next = NULL;

                        // Numero aleatorio entre 0 e num_lists - 1
                        int list_index = rand_r(&seed) % num_lists;

                        // Lock na lista
                        omp_set_lock(&locks[list_index]);

                        new_node->next = list_heads[list_index];
                        list_heads[list_index] = new_node;

                        // Unlock na lista
                        omp_unset_lock(&locks[list_index]);
                    }
                }
            }
        }
    }

    printf("\nLists after processing:\n");
    for (int i = 0; i < num_lists; i++) {
        print_list(list_heads[i], i, N_ITERATIONS);
    }

    // Cleanup
    for (int i = 0; i < num_lists; i++) {
        free_list(list_heads[i]); // Free each list
        omp_destroy_lock(&locks[i]); // Destroy locks
    }
    free(list_heads);
    free(locks);

    return 0;
}