#include "stencil.cuh"
#include <cuda.h>
#include <stdio.h>
#include <random>

int main(int argc, char *argv[]) {
	
	size_t n = atol(argv[1]);
	size_t R = atol(argv[2]);
	size_t threads_per_block = atol(argv[3]);
	
	// set up random number from -1 to 1 generator
	std::random_device entropy_source;
	std::mt19937_64 generator(entropy_source()); 
	const int min = -1.0, max = 1.0; // The range for the random number generator is -1.0 to 1.0
	// there are tons of oter distributino that could be found from https://en.cppreference.com/w/cpp/header/random
	std::uniform_real_distribution<float> dist(min, max);
	
	float *image, *output, *mask;
	size_t mask_size = 2 * R + 1;
	
	// allocate array 
	cudaMallocManaged((void **)&image, n * sizeof(float));
  	cudaMallocManaged((void **)&output, n * sizeof(float));
  	cudaMallocManaged((void **)&mask, mask_size * sizeof(float));
  	
  	for (size_t i = 0; i < n; i++) {
  		image[i] = dist(generator);
  	}
	
	for (size_t i = 0; i < mask_size; ++i) {
    	mask[i] = dist(generator);
  	}
	
	//set up timer
 	cudaEvent_t start;
  	cudaEvent_t stop;
  	cudaEventCreate(&start);
  	cudaEventCreate(&stop);
	
	// start timing and test
	cudaEventRecord(start);
	stencil(image, mask, output, n, R, threads_per_block);
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	float ms;
  	cudaEventElapsedTime(&ms, start, stop);
  	printf("%f\n%f\n", output[n - 1], ms);
  	
  	cudaFree(image);
  	cudaFree(output);
  	cudaFree(mask);
}
