#include <chrono>
#include <cstring>
#include <random>
#include <algorithm>
#include <random>
#include <stdio.h>
#include "cluster.h"

using namespace std;
using namespace chrono;

int main(int argc, char *argv[]) {
    int n = atol(argv[1]);
    int t = atol(argv[2]);	
    
    random_device entropy_source;
	mt19937_64 generator(entropy_source()); 
	const float min = 0.0, max = n; // The range for the random number 		generator is 0.0 to n
	uniform_real_distribution<float> dist(min, max);
    
	float *arr = new float[n];
	float *centers = new float[t];
	float *dists = new float[t];
	
	for (int i = 0; i < n; i++)
		arr[i] = dist(generator);	
	std::sort(arr, arr + n);
	
	for (int i = 1; i <= t; i++) {
		centers[i-1] = (2.0 * i - 1) * n / (2.0 * t);
		dists[i-1] = 0.0;
	}
	
	double T = 0.0;
	int iteration = 10;
	// start timing:
	for (int itr = 0; itr < iteration; itr++) {
		omp_set_num_threads(t);
		auto start = high_resolution_clock::now();
  		cluster(n, t, arr, centers, dists);
  		auto end = high_resolution_clock::now();
  		T += duration_cast<duration<double, std::milli>>(end - start).count();
	}
	
  	
	// looking for the largest elements in the array  	
  	float maxVal = -999.0;
  	int maxPos = 0;
  	for (int i = 0; i < t; i++) {
  		if (dists[i] > maxVal) {
  			maxVal = dists[i];
  			maxPos = i;
  		}
  	}
  	printf("%f\n%d\n%f\n", maxVal, maxPos, T/iteration);
}
