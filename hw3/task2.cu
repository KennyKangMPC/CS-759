#include <cuda.h>
#include <stdio.h>
#include <random>

__global__ void randSumKernel(size_t *arr, int a) {
	
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
	
	
	// initialize host array hA
	int hA[numElement];
	
	
	
	return 0;
}
