#include "matmul.h"

void mmul(const float* A, const float* B, float* C, const std::size_t n){
	// Initialize C contents with zeros
	for (size_t i = 0; i < n * n; i++) {
		C[i] = 0;
	}
	//Another easier way is to use: memset(C, 0x00, n*n);
	
	// collapse(n) collapses the following n nested loops into a single parallel loop shared by the threads
#pragma omp parallel for collapse(2)
	for (size_t i = 0; i < n; ++i) {
		 for (size_t j = 0; j < n; ++j) {
		 	for (size_t k = 0; k < n; ++k) {
		 		C[i * n + j] += A[i * n + k] * B[j * n + k];
		 	}
		 }
	}
}
