#include <iostream>
#include <vector>
#include <ctime>
#include <sys/time.h>
#include <iomanip>

void calcular_matriz(int n) {
    // Input do tamanho da matriz
    struct timeval t0, t1;
    
    // Use C++ vectors instead of malloc
    std::vector<int> matriz(n * n);
    std::vector<int> vetor(n);
    std::vector<int> nova_matriz(n);
    
    int c = 1;
    for (int i = 0; i < n; i++) {
        vetor[i] = i % 10 + 1;
        for (int j = 0; j < n; j++) {
            matriz[i*n + j] = i % 10 + 1;
            c++;
        }
    }
    
    // Warmup
    int soma = 0;
    for (int i = 0; i < n; i++) {
        soma = 0;
        for (int j = 0; j < n; j++) {
            soma += matriz[j*n + i] * vetor[i];
        }
        nova_matriz[i] = soma;
    }
    
    soma = 0;
    gettimeofday(&t0, 0);
    for (int i = 0; i < n; i++) {
        soma = 0;
        for (int j = 0; j < n; j++) {
            soma += matriz[j*n + i] * vetor[i];
        }
        nova_matriz[i] = soma;
    }
    gettimeofday(&t1, 0);
    
    std::cout << "===================================================================" << std::endl;
    std::cout << "Multiplicação de matriz " << n << " por " << n << " por vetor de tamanho " << n << std::endl;
    
    long elapsed = ((t1.tv_sec * 1000000) + t1.tv_usec) - ((t0.tv_sec * 1000000) + t0.tv_usec);
    float time_in_seconds = elapsed / 1000000.0;  // Convert to seconds
    float n_float = static_cast<float>(n);
    float time_per_iteration = time_in_seconds / (n_float * n_float / 1000000.0);
    
    std::cout << std::fixed << std::setprecision(9);
    std::cout << "Total time = " << time_in_seconds << " seconds" << std::endl;
    std::cout << time_per_iteration << " seconds per 1000000 iterations" << std::endl;
    std::cout << std::endl;
    
}

int main() {
    std::vector<int> sizes = {80, 89, 90, 91, 100, 250, 255, 256, 257, 260, 
                             1440, 1447, 1448, 1449, 1460, 5000, 10000, 15000, 
                             20000, 30000, 35000};
    
    for (size_t i = 0; i < sizes.size(); i++) {
        int n = sizes[i];
        calcular_matriz(n);
        std::cout << "===================================================================" << std::endl;
    }
    
    return 0;
}