#include "scan.h"

void scan(const float *arr, float *output, std::size_t n) {
	output[0] = arr[0];
	for (size_t i = 1; i < n; i++){
		output[i] = arr[i] + output[i-1];
	}
}


