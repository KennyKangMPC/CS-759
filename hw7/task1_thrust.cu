#include <stdio.h>
#include <cuda.h>
#include <random>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <iostream>

int main(int argc, char *argv[]) {
	int n = atol(argv[1]);
	
	// set up random number from -1 to 1 generator
	std::random_device entropy_source;
	std::mt19937_64 generator(entropy_source()); 
	const float min = -1.0, max = 1.0; // The range for the random number
	std::uniform_real_distribution<float> dist(min, max);
	
	thrust::host_vector<float> hostVec(n); // creat host vector
	for (int i = 0; i < n; i++) { // fill host vector based on random gen
		hostVec[i] = dist(generator);
	}
	thrust::device_vector<float> deviceVec = hostVec; // copy to device
	//start timer
	// set up timer
  	cudaEvent_t start;
  	cudaEvent_t stop;
  	cudaEventCreate(&start);
  	cudaEventCreate(&stop);
  	
  	//do the timing
  	cudaEventRecord(start);
  	float result = thrust::reduce(deviceVec.begin() , deviceVec.end(), 0.0, thrust::plus<float>());
  	cudaEventRecord(stop);
  	cudaEventSynchronize(stop);
  	
  	float ms;
  	cudaEventElapsedTime(&ms, start, stop);
  	std::cout << result << std::endl;
  	std::cout << ms << std::endl;
}
