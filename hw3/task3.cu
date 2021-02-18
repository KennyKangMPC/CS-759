#include "vscale.cuh"
#include <cuda.h>
#include <stdio.h>
#include <random>

#define NUM_THREADS 512 // another option is 16 based on the problem statement

// reference code is: https://github.com/DanNegrut/ME759/blob/main/2021Spring/Assignments/general/timing.md
int main(int argc, char *argv[]) {
	int n = atoi(argv[1]);
	
	// set up random value generator for both matrix
	std::random_device entropy_source;
	std::mt19937_64 generator(entropy_source()); 
	const float minA = -10.0, maxA = 10.0;
	const float minB = 0.0, maxB = 1.0;
	// there are tons of oter distributino that could be found from https://en.cppreference.com/w/cpp/header/random
	std::uniform_real_distribution<float> distA(minA, maxA);
	std::uniform_real_distribution<float> distB(minB, maxB);
	
	// allocate arrays
	float *a, *b;
	cudaMallocManaged((void **)&a, sizeof(float) * n);
  	cudaMallocManaged((void **)&b, sizeof(float) * n);
  	
  	// initialize with appropriate random value
  	for (int i = 0; i < n; i++) {
  		a[i] = distA(generator);
  		b[i] = distB(generator);
  	}
  	
  	// set up timing variables for cuda
  	cudaEvent_t start;
	cudaEvent_t stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	
	// allocate device
	int device = -1;
	cudaGetDevice(&device);
	cudaMemPrefetchAsync(a, sizeof(float) * n, device, NULL);
  	cudaMemPrefetchAsync(b, sizeof(float) * n, device, NULL);
	
	// calculate number of blocks needed for the number of threads
	int NUM_BLOCKS = (n + NUM_THREADS - 1) / NUM_THREADS;
	
	// timing the kernel call
	cudaEventRecord(start);
	vscale<<<NUM_BLOCKS, NUM_THREADS>>>(a, b, n);
	cudaDeviceSynchronize();
  	cudaEventRecord(stop);
  	cudaEventSynchronize(stop);
		
	// Get the elapsed time in milliseconds
	float ms;
	cudaEventElapsedTime(&ms, start, stop);
	
	// print out time
	std::printf("%f\n%f\n%f\n", ms, b[0], b[n - 1]);
	
	// clean memory
	cudaFree(a);
	cudaFree(b);
}

