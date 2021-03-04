#include "matmul.cuh"
#include <cstdio>

int main(int argc, char *argv[]) {
  unsigned int n = atol(argv[1]);
  unsigned int block_dim = atol(argv[2]);

  float *A, *B, *C;

  cudaMallocManaged((void **)&A, n * n * sizeof(float));
  cudaMallocManaged((void **)&B, n * n * sizeof(float));
  cudaMallocManaged((void **)&C, n * n * sizeof(float));

  for (size_t i = 0; i < n * n; ++i) {
    A[i] = 2;
    B[i] = 0.5;
  }

  cudaEvent_t start;
  cudaEvent_t stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start);
  matmul(A, B, C, n, block_dim);
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (C[i * n + j] != (float) n){
        printf("Error: (%d, %d) is %f\n", i, j, C[i * n + j]);
      }
    }
  }

  float ms;
  cudaEventElapsedTime(&ms, start, stop);
  printf("%f\n%f\n%f\n", C[0], C[n * n - 1], ms);

  cudaFree(A);
  cudaFree(B);
  cudaFree(C);
}
