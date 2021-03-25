#include <omp.h>
#include "matmul.h"
#include <iostream>
#include <cstddef>
#include <random>
#include <stdio.h>

using namespace std;

int main(int argc, char *argv[]) {
	int n = atol(argv[1]);
	int t = atol(argv[2]); // number of threads to be used
	
	random_device entropy_source;
	mt19937_64 generator(entropy_source()); 
	const float min = -1.0, max = 1.0; // The range for the random number generator is -1.0 to 1.0
	uniform_real_distribution<float> dist(min, max);
	
	float *A = (float*)malloc(n * n *sizeof(float));
	float *B = (float*)malloc(n * n * sizeof(float));
	float *C = (float*)malloc(n * n * sizeof(float));
	
	// initialize with random number
	for (int i = 0; i < n * n; i++) {
		A[i] = dist(generator);
		B[i] = dist(generator);
	}
	
	// set number of threads to be used
	omp_set_num_threads(t);
	// let me use omp for timing instead
	double startTime = omp_get_wtime();
	mmul(A, B, C, n);
	double endTime = omp_get_wtime();
	double duration = (endTime - startTime) * 1000;
	printf("%d\n%d\n%d\n", C[0], C[n*n-1], duration);
}
