#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
    
    printf("Hello, World! %d\n", omp_get_num_devices());
    return 0;
}