#include "montecarlo.h"

int montecarlo(const size_t n, const float *x, const float *y, const float radius) {
	int isIn = 0; 
	// below is used in case it happens to be 
#ifdef SIMD
#pragma omp parallel for simd reduction(+ : count)
#else
#pragma omp parallel for reduction(+ : isIn)
#endif
	for (size_t i = 0; i < n; i++) {
		if (x[i] * x[i] + y[i] * y[i] < r * r)
			isIn += 1;
	}
	return count;
}
