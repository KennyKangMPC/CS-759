#include "scan.cuh"
#include <cuda.h>
#include <stdio.h>
#include <random>

int main(int argc, char *argv[]) {
	int n = atol(argv[1]);
 	int threads_per_block = atol(argv[2]);
 	
 	// set up random number from -1 to 1 generator
	std::random_device entropy_source;
	std::mt19937_64 generator(entropy_source()); 
	const float min = -1.0, max = 1.0; // The range for the random number
	std::uniform_real_distribution<float> dist(min, max);
	
	float *input, *output;
	cudaMallocManaged((void **)&input, sizeof(float) * n);
  	cudaMallocManaged((void **)&output, sizeof(float) * n);
	
	for (int i = 0; i < n; i++) {
		input[i] = dist(generator);
	}
	
	/// time for the operations.
	// set up timer
  	cudaEvent_t start;
  	cudaEvent_t stop;
  	cudaEventCreate(&start);
  	cudaEventCreate(&stop);
  	
  	// start timing
  	cudaEventRecord(start);
  	scan(input, output, n, threads_per_block);
  	cudaEventRecord(stop);
  	cudaEventSynchronize(stop);
  	
  	// Get the elapsed time in milliseconds
	float ms;
	cudaEventElapsedTime(&ms, start, stop);
	
	// test
//	for (int i = 0; i < n; i++) {
//		printf("%f, ", input[i]);
//	}
//	printf("\n");
//	for (int i = 0; i < n; i++) {
//		printf("%f, ", output[i]);
//	}
//	printf("\n");
	printf("%f\n%f\n", output[n - 1], ms);
	
	//free memory
	cudaFree(input);
	cudaFree(output);
}
