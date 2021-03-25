#include "msort.h"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <random>
#include <cassert>
#include <stdio.h>

int main(int argc, char *argv[]) {
	int n = atol(argv[1]);
  	int t = atol(argv[2]); // number of threads
  	int ts = atol(argv[3]); // threshold
	
	// set up random number from -1000 to 1000 generator
	std::random_device entropy_source;
	std::mt19937_64 generator(entropy_source()); 
	const int min = -1000, max = 1000; // The range for the random number
	std::uniform_int_distribution<int> dist(min, max);
	
	int *arr = new int[n];
    for (int i = 0; i < n; i++) {
        arr[i] = dist(generator);
    }
    
    omp_set_num_threads(t);
  	omp_set_nested(1);
  	
  	double startTime= omp_get_wtime();
  	msort(arr, n, ts);
  	double endTime= omp_get_wtime();
  	assert(std::is_sorted(arr, arr+n)); // To check if the sort is correct
  	double duration = (endTime - startTime) * 1000;
	printf("%d\n%d\n%d\n", arr[0], arr[n-1], duration);
}

