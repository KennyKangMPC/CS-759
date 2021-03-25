#include "convolution.h"

float getfval(const float *f, std::size_t i, std::size_t j, std::size_t n) {
	// padding
	//return (0 <= i && i < n && 0 <= j && j < n) : f[i*n+j] ? ((0 <= i && i < n) || (0 <= j && j < n)) : 1 ? 0;
	if (0 <= i && i < n && 0 <= j && j < n)
		return f[i*n+j];
	if ((0 <= i && i < n) || (0 <= j && j < n))
		return 1;
	return 0;
	//return (0 <= i && i < n && 0 <= j && j < n) : f[i*n+j] ? ((0 <= i && i < n) || (0 <= j && j < n)) : 1 ? 0;
}
