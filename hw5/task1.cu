#include "reduce.cuh"
#include <cstdio>

int main(int argc, char *argv[]) {
  int n = atol(argv[1]);
  int threads_per_block = atol(argv[2]);

  auto arr = new int[n];

  for (int i = 0; i < n; i++) {
    arr[i] = 1;
  }

  cudaEvent_t start;
  cudaEvent_t stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start);
  int result = reduce(arr, n, threads_per_block);
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);

  if (result != n)
    printf("Error: result is %d instead of %d\n", result, n);

  float ms;
  cudaEventElapsedTime(&ms, start, stop);
  printf("%d\n%f\n", result, ms);
}
