#include "scan.cuh"
#include "stdio.h"


__global__ void inScan() {
		
}

// here the code could also handle more blocks. In lecture note, only 1
__global__ void hs(float *g_odata, float *g_idata, float *int n) {
	extern volatile __shared__  float temp[]; // allocated on invocation
}

// "inclusive scan". Use lecture notes
__host__ void scan(const float* input, float* output, unsigned int n, unsigned int threads_per_block) {
	// allocate cuda memory
	float *g_id, *g_od
	cudaMalloc(&g_id, n * sizeof(float));
  	cudaMalloc(&g_od, n * sizeof(float));
  	
  	// allocate sum cuda memory for each iteration and overall sum
  	// this is needed due to various number of blocks
  	float *g_is, *g_os
  	// need to calculate number of block
	int num_block = (n - 1 + threads_per_block) / threads_per_block;
	cudaMalloc(&g_is, num_block * sizeof(float));
	cudaMalloc(&g_os, num_block * sizeof(float));
	
	// map to device input  for g_id
	cudaMemcpy(g_id, input, n * sizeof(float), cudaMemcpyHostToDevice);
	
	// number of memory for shared size
	int shared_size = 2 * threads_per_block * sizeof(float);
	
	
}
