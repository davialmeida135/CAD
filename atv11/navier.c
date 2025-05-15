#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 32
#define NY 32
#define NZ 32
#define NT 100
#define DX 1.0f
#define DY 1.0f
#define DZ 1.0f
#define DT 0.1f
#define NU 0.1f

// Macros para acessar 3D como 1D
#define IDX(i, j, k) ((i)*NY*NZ + (j)*NZ + (k))

// Função para aplicar o laplaciano com condições periódicas
float laplacian(float *f, int i, int j, int k) {
    int ip = (i + 1) % NX;
    int im = (i - 1 + NX) % NX;
    int jp = (j + 1) % NY;
    int jm = (j - 1 + NY) % NY;
    int kp = (k + 1) % NZ;
    int km = (k - 1 + NZ) % NZ;

    float lap = (f[IDX(ip, j, k)] - 2.0f * f[IDX(i, j, k)] + f[IDX(im, j, k)]) / (DX * DX)
              + (f[IDX(i, jp, k)] - 2.0f * f[IDX(i, j, k)] + f[IDX(i, jm, k)]) / (DY * DY)
              + (f[IDX(i, j, kp)] - 2.0f * f[IDX(i, j, k)] + f[IDX(i, j, km)]) / (DZ * DZ);
    return lap;
}

void step(float *u, float *u_new) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                int idx = IDX(i, j, k);
                u_new[idx] = u[idx] + NU * DT * laplacian(u, i, j, k);
            }
        }
    }
}

float max_deviation(float *u) {
    float max_dev = 0.0f;
    for (int i = 0; i < NX * NY * NZ; i++) {
        float dev = fabsf(u[i] - 1.0f);
        if (dev > max_dev) {
            max_dev = dev;
        }
    }
    return max_dev;
}

int main() {
    size_t size = NX * NY * NZ * sizeof(float);
    float *u  = (float *)malloc(size);
    float *u2 = (float *)malloc(size);

    float *v  = (float *)malloc(size);
    float *v2 = (float *)malloc(size);

    float *w  = (float *)malloc(size);
    float *w2 = (float *)malloc(size);

    // Inicializar campos com valor constante + perturbação no centro
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                int idx = IDX(i, j, k);
                u[idx] = v[idx] = w[idx] = 1.0f;
            }
        }
    }
    // Perturbação no centro
    u[IDX(NX/2, NY/2, NZ/2)] += 15f;

    // Loop de tempo
    for (int t = 0; t < NT; t++) {
        step(u, u2);
        step(v, v2);
        step(w, w2);

        // Swap buffers
        float *tmp;

        tmp = u;  u = u2;  u2 = tmp;
        tmp = v;  v = v2;  v2 = tmp;
        tmp = w;  w = w2;  w2 = tmp;

        if (t % 10 == 0) {
            float dev = max_deviation(u);
            printf("Passo %d: desvio máximo da velocidade u = %.6f\n", t, dev);
        }
    }

    printf("\nDesvio final: %.6f\n", max_deviation(u));

    free(u);  free(u2);
    free(v);  free(v2);
    free(w);  free(w2);
    return 0;
}
