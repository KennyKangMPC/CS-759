#include <cuda.h>
#include <stdio.h>


__global__ void  factKernel(){
	int a = threadIdx.x + 1;
	int b = 1;
	for (int i = 1; i <= a; i++) {
		b *= i;
	}
	printf("%d != %d\n", a, b);		
}

int main() {
	const int numThreads = 8;
	const int numBlocks = 1;
	// invoke GPU kernel, with one block that has eight threads
  	factKernel<<<numBlocks, numThreads>>>();
	cudaDeviceSynchronize();
	return 0;
}
