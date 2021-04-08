#include "reduce.h"

float reduce(const float* arr, const size_t l, const size_t r) {
	float output = 0;
#pragma omp parallel for simd reduction(+ : output)
	for (size_t i = l; i < r; i++)
		output += arr[i];
	return output;
}
