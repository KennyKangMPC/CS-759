#include "matmul.cuh"
#include <cuda.h>
#include <stdio.h>
#include <random>

int main(int argc, char *argv[]) {
	// obtain user input
	size_t n = atol(argv[1]);
	size_t threads_per_block = atol(argv[2]);
	
	// set up random number from -1 to 1 generator
	std::random_device entropy_source;
	std::mt19937_64 generator(entropy_source()); 
	const int min = -1.0, max = 1.0; // The range for the random number generator is -1.0 to 1.0
	// there are tons of oter distributino that could be found from https://en.cppreference.com/w/cpp/header/random
	std::uniform_real_distribution<float> distA(min, max);
	std::uniform_real_distribution<float> distB(min, max);

	// allocate array
	float *ma, *mb, *mc;
	// device array

	cudaMallocManaged((void **)&ma, sizeof(float) * n * n);
  	cudaMallocManaged((void **)&mb, sizeof(float) * n * n);
  	cudaMallocManaged((void **)&mc, sizeof(float) * n * n);

	// insert random initial value into it
	for (size_t i = 0; i < n * n; i++) {
		ma[i] = distA(generator);
		mb[i] = distB(generator);
	}
	
	// set up timer
  	cudaEvent_t start;
  	cudaEvent_t stop;
  	cudaEventCreate(&start);
  	cudaEventCreate(&stop);

  	// record time
  	cudaEventRecord(start);
  	matmul(ma, mb, mc, n, threads_per_block);
  	cudaEventRecord(stop);
  	cudaEventSynchronize(stop);
  	
  	// Get the elapsed time in milliseconds
	float ms;
	cudaEventElapsedTime(&ms, start, stop);
	
	// print out the last element of c and the time
  	printf("%f\n%f\n", mc[n * n - 1], ms);
  	
  	// clearn memory
  	cudaFree(ma);
  	cudaFree(mb);
  	cudaFree(mc); 	
}
