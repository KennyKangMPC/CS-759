#include "mmul.h""

void mmul(cublasHandle_t handle, const float* A, const float* B, float* C, int n) {
	
	// The most approriate function I found from https://docs.nvidia.com/cuda/cublas/index.html#cublas-lt-t-gt-gemm
	// is level three 2.7.1 cublas<t>gemm() function with \alpha = \beta = 1, <t> = S
	
	// notice here, since op(A) = A, transa == CUBLAS_OP_N same for op(B) = B
	const float ALPHA = 1;
	const float BETA = 1;
	const float *alpha = &ALPHA;
	const float *beta = &BETA;
	cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, alpha, A, n, B, n, beta, C, n);
	cudaDeviceSynchronize();
}
