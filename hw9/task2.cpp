#include <algorithm>
#include <chrono>
#include <iostream>
#include <omp.h>
#include <random>
#include <stdio.h>
#include "montecarlo.h"

using namespace std;

int main(int argc, char *argv[]) {
    int n = atol(argv[1]);
    int t = atol(argv[2]);
    
    random_device entropy_source;
	mt19937_64 generator(entropy_source()); 
	const float min = -r, max = r; // The range for the random number 		generator is 0.0 to n
	uniform_real_distribution<float> dist(min, max);
    
    float r = 1.0;
    float *x = new float[n];
    float *y = new float[n];
    
    for (int i = 0; i < n; i++) {
    	x[i] = dist(generator);
    	y[i] = dist(generator);
    }
    
    // start timing
    auto pi_estimate;
    double time = 0.0;
    int iteration = 10;
    for (int itr = 0; itr < iteration; itr++) {
    	omp_set_num_threads(t);
    	double startTime = omp_get_wtime();
		pi_estimate = 4.0 * montecarlo(n, x, y, r) / n;
		double endTime = omp_get_wtime();
		time += (endTime - startTime) * 1000;
    }
    
    printf("%f\n%f\n", pi_estimate, time/iteration);
}
