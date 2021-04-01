#include "montecarlo.h"

int montecarlo(const size_t n, const float *x, const float *y, const float radius) {
	int isIn = 0; 
	// below is used in case it happens to be Single thread
#ifdef SIMD
#pragma omp parallel for simd reduction(+ : isIn)
#else
#pragma omp parallel for reduction(+ : isIn)
#endif
	for (size_t i = 0; i < n; i++) {
		if (x[i] * x[i] + y[i] * y[i] < radius * radius)
			isIn += 1;
	}
	return isIn;
}
