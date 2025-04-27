#include <stdio.h>
#include <stdlib.h>
#include <time.h>       // For time(), srand()
#include <sys/time.h>   // For gettimeofday()
#include <omp.h>        // For OpenMP

// --- Configuration ---
// Make NUM_ELEMENTS large enough to likely exceed your CPU's L3 cache
// L3 cache might be 8MB, 16MB, 32MB, etc.
// Size of Node is 16 bytes (on 64-bit: 8 for pointer, 8 for long long).
// 100 million elements * 16 bytes/element = 1.6 GB. Adjust if needed.
#define NUM_ELEMENTS 10000000 // Number of nodes in our linked list pool
#define CHASES_PER_THREAD 50000000 // How many pointers each thread will chase

// --- Data Structure ---
typedef struct Node {
    struct Node* next;
    long long data; // Some payload data
} Node;

// --- Fisher-Yates Shuffle ---
// Used to randomize the pointer chain
void shuffle_indices(int* array, size_t n) {
    if (n > 1) {
        for (size_t i = n - 1; i > 0; i--) {
            size_t j = rand() % (i + 1); // Generate random index j such that 0 <= j <= i
            // Swap array[i] and array[j]
            int temp = array[i];
            array[i] = array[j];
            array[j] = temp;
        }
    }
}

int main() {
    int num_threads = 1;
    printf("Enter the number of threads: ");
    scanf("%d", &num_threads);
    printf("Using %d threads.\n", num_threads);
    printf("Total elements: %ld\n", (long)NUM_ELEMENTS);
    printf("Chases per thread: %ld\n", (long)CHASES_PER_THREAD);

    // --- Allocation ---
    printf("Allocating memory...\n");
    Node* nodes = (Node*)malloc(NUM_ELEMENTS * sizeof(Node));
    if (nodes == NULL) {
        perror("Failed to allocate memory for nodes");
        return 1;
    }

    // --- Initialization and Random Chaining ---
    printf("Initializing and creating random chain...\n");
    srand(time(NULL)); // Seed random number generator

    // Create an array of indices [0, 1, 2, ..., N-1]
    int* indices = (int*)malloc(NUM_ELEMENTS * sizeof(int));
    if (indices == NULL) {
        perror("Failed to allocate memory for indices");
        free(nodes);
        return 1;
    }
    for (long i = 0; i < NUM_ELEMENTS; ++i) {
        indices[i] = i;
        nodes[i].data = i; // Initialize data payload
    }

    // Shuffle the indices
    shuffle_indices(indices, NUM_ELEMENTS);

    // Link nodes based on shuffled order to create a pseudo-random chain
    // nodes[indices[i]] points to nodes[indices[i+1]]
    for (long i = 0; i < NUM_ELEMENTS - 1; ++i) {
        nodes[indices[i]].next = &nodes[indices[i + 1]];
    }
    // Make the last element point back to the first to form a cycle
    nodes[indices[NUM_ELEMENTS - 1]].next = &nodes[indices[0]];

    free(indices); // No longer need the index array

    // --- Pointer Chasing Measurement ---
    printf("Starting pointer chasing...\n");
    struct timeval start, end;
    volatile long long total_sum = 0; // Use volatile to reduce chance of optimization

    gettimeofday(&start, NULL);
    
    // Each thread starts at a different point in the cycle
    long long local_sum = 0; // Avoid false sharing on a shared sum
    
    // The core latency-bound loop
    #pragma omp parallel num_threads(num_threads)
    for (long long i = 0; i < CHASES_PER_THREAD; ++i) {
        Node* current = &nodes[i%NUM_ELEMENTS]; // Start at a different node for each thread
        current = current->next; // <<< THE LATENCY-BOUND STEP
        local_sum += (long long)current; // Add pointer address to prevent optimization
    }

    // Ensure the compiler doesn't optimize away the loop
    // by using the final 'current' pointer and adding to a thread-local sum.


    gettimeofday(&end, NULL);

    // --- Results ---
    long seconds = end.tv_sec - start.tv_sec;
    long microseconds = end.tv_usec - start.tv_usec;
    if (microseconds < 0) {
        seconds--;
        microseconds += 1000000;
    }
    double elapsed = seconds + microseconds * 1e-6;

    printf("Pointer chasing took %f seconds\n", elapsed);
    // printf("Dummy checksum (volatile): %lld\n", total_sum); // Optional print

    // --- Cleanup ---
    printf("Freeing memory...\n");
    free(nodes);

    return 0;
}