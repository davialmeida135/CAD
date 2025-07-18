# Relatório atividade 23
- Aluno: Davi Diógenes Ferreira de Almeida
- Disciplina: IMD1116 Computação de Alto Desempenho/2025.1
Para esta atividade, implementamos um programa em Cuda que simula a propagação de uma perturbação em um espaço tridimensional por meio da equação de navier-strokes. O programa foi executado com uma GPU com 8GB de VRAM e em um espaço tridimensional de 512x512x512 floats e 10 time steps. 
Para cada bloco de threads, foi utilizada a função __shared__ para se aproveitar da memória compartilhada entre as threads dentro do bloco. Cada thread é responsável por inserir o valor da matriz em que ela vai operar na memória compartilhada, de modo que, após todas as threads terem o feito (após syncthreads), teremos todos os dados que serão trabalhados pelas threads do bloco dentro da memória compartilhada do bloco.

O tempo de execução do programa, em CPU e com 16 threads, foi de 2.629803 segundos. Já o programa em Cuda, com o mesmo espaço tridimensional, alcançou um tempo de 0.104353 segundos. Para validar os cálculos realizados, foram calculadas algumas métricas comparando os espaços finais das duas execuções:
Maximum absolute error: 1.78e-15
Mean absolute error: 3.13e-18
Root mean square error: 1.87e-17

Taxas de erro tão pequenas validam o programa e os tempos de execução mostram que, para este problema, a computação em GPU é um bom caminho.

Abaixo está a função de atualização do espaço:

```c
__global__ void update_temperature_field(double* temperature_field_new, double* temperature_field_old, 
    int grid_points_x, int grid_points_y, int grid_points_z, double diffusion_alpha) {

    __shared__ double tile[BLOCK_SIZE_Z + 2*RADIUS][BLOCK_SIZE_Y + 2*RADIUS][BLOCK_SIZE_X + 2*RADIUS];
    
    // Coordenadas da thread (e do ponto no espaço)
    int thread_x = blockIdx.x * blockDim.x + threadIdx.x;
    int thread_y = blockIdx.y * blockDim.y + threadIdx.y;
    int thread_z = blockIdx.z * blockDim.z + threadIdx.z;
    
    // Posição da thread dentro do bloco
    int local_x = threadIdx.x;
    int local_y = threadIdx.y;
    int local_z = threadIdx.z;
    
    // Posição do tile (thread+raio)
    int tile_x = local_x + RADIUS;
    int tile_y = local_y + RADIUS;
    int tile_z = local_z + RADIUS;
    
    // Carregar o ponto da thread
    if (thread_x < grid_points_x && thread_y < grid_points_y && thread_z < grid_points_z) {
        int global_idx = thread_z * grid_points_y * grid_points_x + thread_y * grid_points_x + thread_x;
        tile[tile_z][tile_y][tile_x] = temperature_field_old[global_idx];
    } else {
        tile[tile_z][tile_y][tile_x] = 1.0;
    }
    
    // Carregar os halos
    // Eixo X
    // Apenas a primeira thread de cada dimensão carregam os halos, não a primeira e a ultima
    if (local_x < RADIUS) {
        // Da esquerda (posicao - raio)
        int halo_x = thread_x - RADIUS;
        if (halo_x >= 0 && thread_y < grid_points_y && thread_z < grid_points_z) {
            int halo_idx = thread_z * grid_points_y * grid_points_x + thread_y * grid_points_x + halo_x;
            tile[tile_z][tile_y][local_x] = temperature_field_old[halo_idx];
        } else {
            tile[tile_z][tile_y][local_x] = 1.0;
        }
        
        // Da dureita (posicao + tamanho)
        int right_x = thread_x + blockDim.x;
        if (right_x < grid_points_x && thread_y < grid_points_y && thread_z < grid_points_z) {
            int right_idx = thread_z * grid_points_y * grid_points_x + thread_y * grid_points_x + right_x;
            tile[tile_z][tile_y][local_x + blockDim.x + RADIUS] = temperature_field_old[right_idx];
        } else {
            tile[tile_z][tile_y][local_x + blockDim.x + RADIUS] = 1.0;
        }
    }
    // Eixo Y
    if (local_y < RADIUS) {
        // De frente (posicao - raio)
        int halo_y = thread_y - RADIUS;
        if (halo_y >= 0 && thread_x < grid_points_x && thread_z < grid_points_z) {
            int halo_idx = thread_z * grid_points_y * grid_points_x + halo_y * grid_points_x + thread_x;
            tile[tile_z][local_y][tile_x] = temperature_field_old[halo_idx];
        } else {
            tile[tile_z][local_y][tile_x] = 1.0;
        }
        
        // De trás (posicao + tamanho)
        int back_y = thread_y + blockDim.y;
        if (back_y < grid_points_y && thread_x < grid_points_x && thread_z < grid_points_z) {
            int back_idx = thread_z * grid_points_y * grid_points_x + back_y * grid_points_x + thread_x;
            tile[tile_z][local_y + blockDim.y + RADIUS][tile_x] = temperature_field_old[back_idx];
        } else {
            tile[tile_z][local_y + blockDim.y + RADIUS][tile_x] = 1.0;
        }
    }
    // Eixo Z
    if (local_z < RADIUS) {
        // De cima (posicao - raio)
        int halo_z = thread_z - RADIUS;
        if (halo_z >= 0 && thread_x < grid_points_x && thread_y < grid_points_y) {
            int halo_idx = halo_z * grid_points_y * grid_points_x + thread_y * grid_points_x + thread_x;
            tile[local_z][tile_y][tile_x] = temperature_field_old[halo_idx];
        } else {
            tile[local_z][tile_y][tile_x] = 1.0;
        }
        
        // De baixo (posicao - raio)
        int top_z = thread_z + blockDim.z;
        if (top_z < grid_points_z && thread_x < grid_points_x && thread_y < grid_points_y) {
            int top_idx = top_z * grid_points_y * grid_points_x + thread_y * grid_points_x + thread_x;
            tile[local_z + blockDim.z + RADIUS][tile_y][tile_x] = temperature_field_old[top_idx];
        } else {
            tile[local_z + blockDim.z + RADIUS][tile_y][tile_x] = 1.0;
        }
    }
    __syncthreads();
    // Calcular somente para pontos internos (nao fazer para os halos)
    if (thread_x > 0 && thread_x < grid_points_x-1 &&
        thread_y > 0 && thread_y < grid_points_y-1 &&
        thread_z > 0 && thread_z < grid_points_z-1) {
        
        int global_idx = thread_z * grid_points_y * grid_points_x + thread_y * grid_points_x + thread_x;
        
        double center = tile[tile_z][tile_y][tile_x];
        double neighbors = tile[tile_z][tile_y][tile_x-1] + tile[tile_z][tile_y][tile_x+1] +  // direita - esquerda
                          tile[tile_z][tile_y-1][tile_x] + tile[tile_z][tile_y+1][tile_x] +  // cima - baixo  
                          tile[tile_z-1][tile_y][tile_x] + tile[tile_z+1][tile_y][tile_x];   // frente - trás
        
        temperature_field_new[global_idx] = center + diffusion_alpha * (neighbors - 6.0 * center);
    }
}
```

