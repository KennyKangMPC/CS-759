#include "msort.h"
#include <algorithm>

void msortRecursive(int *arr, const std::size_t n, const std::size_t threshold, int nT) {
	if (n <= 1) {return;} // base case size of array is 1 or less
	
	// If the size of array goes below	
	// the threshold, a serial sort algorithm will be used to avoid overhead
	// of task scheduling
	if (n < threshold) {
		// we need to search separate half part of the array
		for (auto pt = arr; pt != arr + n; pt++) { // stop at the mid point
			// find the first element larger than pt and then stop there
			// and then shift to the rest part of the array
			std::rotate(std::upper_bound(arr, pt, *pt), pt, pt+1);
		}
		return;
	}
	
	if (nT == 1) { // only one thread
		msortRecursive(arr, n/2, threshold, nT);
        msortRecursive(arr+n/2, n - n/2, threshold, nT);
	} else {
 		//half threads sort first half, the rest use for another part
	#pragma omp task
		msortRecursive(arr, n/2, threshold, nT/2);
	#pragma omp task
		msortRecursive(arr+n/2, n- n/2, threshold, nT - nT/2);
	#pragma omp taskwait
	}
	// merge
	std::inplace_merge(arr, arr + n/2, arr + n);
}

// "threshold" is the lower limit of array size where your function would 
// start making parallel recursive calls. If the size of array goes below
// the threshold, a serial sort algorithm will be used to avoid overhead
// of task scheduling
void msort(int* arr, const std::size_t n, const std::size_t threshold) {
#pragma omp parallel
#pragma omp single
	msortRecursive(arr, n, threshold, omp_get_num_threads());
}
