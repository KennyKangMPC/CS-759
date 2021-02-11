#include "matmul.h"
#include <iostream>

void mmul1(const double* A, const double* B, double* C, const unsigned int n) {
	unsigned int i, j, k;
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++) {
			for (k = 0; k < n; k++) {
				C[i * n + j] += A[i * n + k] * B[k * n + j];
			}
		}
	}
}


void mmul2(const double* A, const double* B, double* C, const unsigned int n) {
	unsigned int i, j, k;
	for (i = 0; i < n; i++){
		for (k = 0; k < n; k++) {
			for (j = 0; j < n; j++) {
				C[i * n + j] += A[i * n + k] * B[k * n + j];
			}
		}
	}
}


void mmul3(const double* A, const double* B, double* C, const unsigned int n) {
	unsigned int i, j, k;
	for (j = 0; j < n; j++) {
		for (k = 0; k < n; k++) {
			for (i = 0; i < n; i++){
				C[i * n + j] += A[i * n + k] * B[k * n + j];
			}
		}
	}
}


void mmul4(const std::vector<double>& A, const std::vector<double>& B, double* C, const unsigned int n) {
	unsigned int i, j, k;
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++) {
			for (k = 0; k < n; k++) {
				
				C[i * n + j] += A[i * n + k] * B[k * n + j];
			}
		}
	}
}


