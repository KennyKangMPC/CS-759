#include "reduce.cuh"

__global__ void reduce_kernel(const int *g_idata, int *g_odata,
                              unsigned int n) {
  extern __shared__ int sdata[];

  int i = blockIdx.x * blockDim.x + threadIdx.x;
  sdata[threadIdx.x] = i < n ? g_idata[i] : 0;

  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (threadIdx.x < s)
      sdata[threadIdx.x] += sdata[threadIdx.x + s];
    __syncthreads();
  }

  if (threadIdx.x == 0)
    g_odata[blockIdx.x] = sdata[0];
}

__host__ int reduce(const int *arr, unsigned int N, unsigned int t_len) {
  int *g_idata, *g_odata;

  cudaMalloc(&g_idata, N * sizeof(int));
  cudaMalloc(&g_odata, (N + t_len - 1) / t_len * sizeof(int));
  cudaMemcpy(g_idata, arr, N * sizeof(int), cudaMemcpyHostToDevice);

  for (int n = N; n != 1; n = (n + t_len - 1) / t_len) {
    size_t n_blk = (n + t_len - 1) / t_len;
    reduce_kernel<<<n_blk, t_len, t_len * sizeof(int)>>>(g_idata, g_odata, n);
    cudaMemcpy(g_idata, g_odata, n_blk * sizeof(int), cudaMemcpyDeviceToDevice);
  }
  cudaDeviceSynchronize();

  int result;
  cudaMemcpy(&result, g_odata, sizeof(int), cudaMemcpyDeviceToHost);

  cudaFree(g_idata);
  cudaFree(g_odata);

  return result;
}