#include <cuda.h>
#include <stdio.h>
#include <random>

__global__ void randSumKernel(int *arr, int a) {
	// threadIdx.x is x and blockIdx.x is y
	arr[blockIdx.x * blockDim.x + threadIdx.x] = a * threadIdx.x + blockIdx.x;
}
	
// reference is https://github.com/DanNegrut/ME759/blob/main/2021Spring/GPU/setArray.cu
int main(){
	const int numBlocks = 2;
	const int numThreads = 8;
	const int numElement = 16;
	
	// set up random number generator
	std::random_device entropy_source;
	std::mt19937_64 generator(entropy_source()); 
	const int min = 0, max = 10; // The range for the random number generator is 0 to 10
	// there are tons of oter distributino that could be found from https://en.cppreference.com/w/cpp/header/random
	std::uniform_int_distribution<> dist(min, max);
	
	// use random number generator to generate integer a
	int a = dist(generator);
	
	// initialize device array dA
	int *dA;
	// allocate memory on the device; zero out all entries in this device array
  	cudaMalloc((void **)&dA, sizeof(int) * numElement);
  	cudaMemset(dA, 0, numElement * sizeof(int));
  	
	// initialize host array hA
	int hA[numElement];
	
	// invoke GPU kernel with 2 blocks that has eight threads
	randSumKernel<<<numBlocks, numThreads>>>(dA, a);
	cudaDeviceSynchronize();
	
	// bring the result back from the GPU into the hostArray
	cudaMemcpy(hA, dA, sizeof(int) * numElement, cudaMemcpyDeviceToHost);
	
	// free array
	cudaFree(dA);
	
	// print out element
  	for (int i = 0; i < numElement - 1; i++) {
    	std::printf("%d ", hA[i]);
  	}
  	std::printf("%d\n", hA[numElement - 1]);
	
	return 0;
}
