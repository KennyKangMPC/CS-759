#include "stencil.cuh"

// Reference from lecture 11 notess

__global__ void stencil_kernel(const float* image, const float* mask, float* output, unsigned int n, unsigned int R) {
	extern __shared__ float shMemArray[];
	
	int t_id = threadIdx.x;
	int t_len = blockDim.x;
	int i = blockIdx.x * t_len + t_id;
	
	if (i >= n) {
		return;
	}
	
	// Copy mask to shared memory
	float *shared_mask = shMemArray;
	if (t_id < 2*R+1) {
		shared_mask[t_id] = mask[t_id]; 
	}
	
	// Initialize shared output
	float *shared_output = shared_mask + 2*R + 1;
	shared_output[t_id] = 0;
	
	// Copy image to shared memory
	float *shared_image = shared_output + t_len + R;
	shared_image[t_id] = image[i];
	
	if (t_id < R) {
		shared_image[t_id - R] = i - R > 0 ? image[i - R] : 0;
	} else if (t_len - t_id < R) {
		shared_image[t_id + R] = i + R < n ? image[i + R] : 0;
	}
	
	__syncthreads();
	
	for (int j = -R; j <= R; j++) {
		shared_output[t_id] += shared_image[t_id + j] * shared_mask[j + R];
	}
	
	// copy back
	output[i] = shared_output[t_id];
}

__host__ void stencil(const float* image, const float* mask, float* output, unsigned int n, unsigned int R, unsigned int threads_per_block) {
	int numBlock = (n - 1 + threads_per_block) / threads_per_block;
	int shared_size = (1 + 4*R + 2*threads_per_block) * sizeof(float);
	stencil_kernel<<<numBlock, threads_per_block, shared_size>>>(image, mask, output, n, R);
	cudaDeviceSynchronize();
}
