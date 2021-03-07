#include "mmul.h"
#include <cuda.h>
#include <stdio.h>
#include <random>
#include <cublas_v2.h>

int main(int argc, char *argv[]) {
	// obtain commandline input
	int n = atol(argv[1]);
  	int n_tests = atol(argv[2]);
  	
  	// set up random number from -1 to 1 generator
	std::random_device entropy_source;
	std::mt19937_64 generator(entropy_source()); 
	const float min = -1.0, max = 1.0; // The range for the random number
	std::uniform_real_distribution<float> dist(min, max);
	
  	// allocate array
	float *a, *b, *c;
	cudaMallocManaged((void **)&a, sizeof(float) * n * n);
  	cudaMallocManaged((void **)&b, sizeof(float) * n * n);
  	cudaMallocManaged((void **)&c, sizeof(float) * n * n);
  	
  	// insert random initial value into it
	for (int i = 0; i < n * n; i++) {
		a[i] = dist(generator);
		b[i] = dist(generator);
		c[i] = dist(generator);
	}
	
	// setup the use of cublas
	// Create a handle for CUBLAS
	cublasHandle_t handle;
	cublasCreate(&handle);
	
	/// time for the operations.
	// set up timer
  	cudaEvent_t start;
  	cudaEvent_t stop;
  	cudaEventCreate(&start);
  	cudaEventCreate(&stop);
  	
  	cudaEventRecord(start);
  	for (int i = 0; i < n_tests; i++) {
  		mmul(handle, a, b, c, n);
  	}
	cudaEventRecord(stop);
  	cudaEventSynchronize(stop);
  	// Get the elapsed time in milliseconds
	float ms;
	cudaEventElapsedTime(&ms, start, stop);
	float averageTime = ms/n_tests;
	printf("%f\n", averageTime);
	
	//clean up everything
	cublasDestroy(handle);
	cudaFree(a);
	cudaFree(b);
	cudaFree(c);
}
