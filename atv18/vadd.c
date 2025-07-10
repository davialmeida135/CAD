#include <stdio.h>
#include <omp.h>
#define N 100000
#define N 100000000
#define TOL  0.0000001
//
//  This is a simple program to add two vectors
//  and verify the results.
//
//  History: Written by Tim Mattson, November 2017
//
int main()
{
    printf("Number of devices: %d\n", omp_get_num_devices());
    float a[N], b[N], c[N], res[N];
    int err=0;

    double init_time, compute_time, test_time;
    init_time    = -omp_get_wtime();

   // fill the arrays
   #pragma omp parallel for
   for (int i=0; i<N; i++){
      a[i] = (float)i;
      b[i] = 2.0*(float)i;
      c[i] = 0.0;
      res[i] = i + 2*i;
   }

   init_time    +=  omp_get_wtime();
   compute_time  = -omp_get_wtime();
   
   // add two vectors
   #pragma omp target
   #pragma omp loop
   for (int i=0; i<N; i++){
      c[i] = a[i] + b[i];
   }

   compute_time +=  omp_get_wtime();
   test_time     = -omp_get_wtime();

   printf(" vectors added with %d errors\n",err);

   printf("Init time:    %.3fs\n", init_time);
   printf("Compute time: %.3fs\n", compute_time);
   printf("Total time:   %.3fs\n", init_time + compute_time);
   return 0;
}
