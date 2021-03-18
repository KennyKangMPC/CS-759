#include <stdio.h>
#include <cuda.h>
#include <random>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <iostream>
#include "count.cuh"

int main(int argc, char *argv[]) {
	int n = atol(argv[1]);
	
	// set up random number from -1 to 1 generator
	std::random_device entropy_source;
	std::mt19937_64 generator(entropy_source()); 
	const int min = 0, max = 500; // The range for the random number
	std::uniform_int_distribution<int> dist(min, max);
	
	thrust::host_vector<int> h_in(n);
	for (int i = 0; i < n; i++) { // fill host vector based on random gen
		h_in[i] = dist(generator);
	}
	thrust::device_vector<int> d_in = h_in;
	thrust::device_vector<int> values(n);
  	thrust::device_vector<int> counts(n);
  	
  	// set up timer
  	cudaEvent_t start;
  	cudaEvent_t stop;
  	cudaEventCreate(&start);
  	cudaEventCreate(&stop);
  	
  	cudaEventRecord(start);
  	count(d_in, values, counts);
  	cudaEventRecord(stop);
  	cudaEventSynchronize(stop);
  	
  	float ms;
  	cudaEventElapsedTime(&ms, start, stop);
  	
  	std::cout << values.back() << std::endl;
  	std::cout << counts.back() << std::endl;
  	std::cout << ms << std::endl;
}
