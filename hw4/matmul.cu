#include "matmul.cuh"

// Computes the matrix product of A and B, storing the result in C.
// Each thread should compute _one_ element of output.
// Does not use shared memory.
//
// A, B, and C are row major representations of nxn matrices.
//
// Assumptions:
// - 1D kernel configuration
__global__ void matmul_kernel(const float* A, const float* B, float* C, size_t n) {
	// check if the element is oversized
	size_t pos = blockIdx.x * blockDim.x + threadIdx.x;
	size_t size = n * n;
	if (pos < size) {
		for (size_t k = 0; k < n; k++) {
			C[pos] += A[pos + k] * B[k * n + pos % n];
		}
	} else {
		return;
	}
}


// Makes one call to matmul_kernel with threads_per_block threads per block.
// The kernel call should be followed by a call to cudaDeviceSynchronize for timing purposes.
void matmul(const float* A, const float* B, float* C, size_t n, unsigned int threads_per_block) {
	// fills  memory
	memset(C, 0, n*n*sizeof(float));
	
	// calculate num_block
	size_t num_block = (threads_per_block - 1 + n * n) / threads_per_block;
	
	// call kernel to do multiplications 
	matmul_kernel<<<num_block, threads_per_block>>>(A, B, C, n);
	
	// for timing purposes
	cudaDeviceSynchronize();
}
