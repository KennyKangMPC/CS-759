#include "vscale.cuh"

__global__ void vscale(const float *a, float *b, unsigned int n) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (i < n) {
		b[i] *= a[i];
	}
}
