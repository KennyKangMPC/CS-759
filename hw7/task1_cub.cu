#define CUB_STDERR // print CUDA runtime errors to console
#include <stdio.h>
#include <cub/device/device_reduce.cuh>
#include <cub/util_allocator.cuh>
#include <random>
//#include "cub/util_debug.cuh" // this is not needed here
using namespace cub;
CachingDeviceAllocator  g_allocator(true);  // Caching allocator for device memory

// note that I copied and pasted the code from the given github page and then modified the code based on the instruction directly on the code.
int main(int argc, char *argv[]) {
    int n = atol(argv[1]);
    
    // set up random number from -1 to 1 generator
	std::random_device entropy_source;
	std::mt19937_64 generator(entropy_source()); 
	const float min = -1.0, max = 1.0; // The range for the random number
	std::uniform_real_distribution<float> dist(min, max);
    
    // Set up host arrays
    float *h_in = new float[n];
    for (int i = 0; i < n; i++)
    	h_in[i] = dist(generator);

    // Set up device arrays
    float* d_in = NULL;
    CubDebugExit(g_allocator.DeviceAllocate((void **)&d_in, sizeof(float) * n));
    // Initialize device input
    CubDebugExit(cudaMemcpy(d_in, h_in, sizeof(float) * n, cudaMemcpyHostToDevice));
    
    // Setup device output array
    float* d_sum = NULL;
    CubDebugExit(g_allocator.DeviceAllocate((void **)&d_sum, sizeof(float) * 1));
    
    // Request and allocate temporary storage
    void* d_temp_storage = NULL;
    size_t temp_storage_bytes = 0;
    CubDebugExit(DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, d_in, d_sum, n));
    CubDebugExit(g_allocator.DeviceAllocate(&d_temp_storage, temp_storage_bytes));

    // Do the actual reduce operation
	cudaEvent_t start;
  	cudaEvent_t stop;
  	cudaEventCreate(&start);
  	cudaEventCreate(&stop);
  	cudaEventRecord(start);
	CubDebugExit(DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, d_in, d_sum, n));
	cudaEventRecord(stop);
  	cudaEventSynchronize(stop);
	
	float result;
	CubDebugExit(cudaMemcpy(&result, d_sum, sizeof(float), cudaMemcpyDeviceToHost));
	
	float ms;
  	cudaEventElapsedTime(&ms, start, stop);
  	
  	printf("%f\n%f\n", result, ms);
  	
    // Cleanup
    if (d_in) CubDebugExit(g_allocator.DeviceFree(d_in));
    if (d_sum) CubDebugExit(g_allocator.DeviceFree(d_sum));
    if (d_temp_storage) CubDebugExit(g_allocator.DeviceFree(d_temp_storage));
    
    return 0;
}
