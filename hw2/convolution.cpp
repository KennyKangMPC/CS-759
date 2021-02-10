#include "convolution.h"

// Computes the result of applying a mask to an image as in the convolution process described in HW02.pdf.
// image is an nxn grid stored in row-major order.
// mask is an mxm grid stored in row-major order.
// Stores the result in output, which is an nxn grid stored in row-major order.
void convolve(const float *image, float *output, std::size_t n, const float *mask, std::size_t m) {
	// recall that f is image and w is mask
	// here we use 1D array to store the matrix
	for (std::size_t x = 0; x < n; x++) {
		for (std::size_t y = 0; y < n; y++) {
			output[x * n + y] = 0; // initialize value to be 0
			for (std::size_t i = 0; i < m - 1; i++) {
				for (std::size_t j = 0; j < m - 1; j++) {
					output[x * n + y] += mask[i * m + j] * getfval(image, x + i - (m-1)/2, y + j - (m-1)/2, n);
				}
			}
		}
	}
}

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
