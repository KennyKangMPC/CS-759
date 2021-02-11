#include "scan.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <random>
// I searched how to do timing from timing.md
using namespace std;
using namespace chrono;

int main(int argc, char *argv[]){
	// set up random value generator for both matrix
	random_device entropy_source;
	mt19937_64 generator(entropy_source()); 
	const float min = -1.0, max = 1.0; // The range for the random number generator is -1.0 to 1.0
	// there are tons of oter distributino that could be found from https://en.cppreference.com/w/cpp/header/random
	uniform_real_distribution<float> dist(min, max);
	
	size_t n = atol(argv[1]);
	float *arr = (float*)malloc(n * sizeof(float));
	float *output = (float*)malloc(n * sizeof(float));
	
	for (size_t i = 0; i < n; i++) {
		arr[i] = dist(generator);
	}
	
	auto start = high_resolution_clock::now();
	scan(arr, output, n);
	auto stop = high_resolution_clock::now();
	
	//auto duration = duration_cast<milliseconds>(stop - start);
	auto duration_sec = duration_cast<duration<double, std::milli>>(stop - start);
	cout << "Time taken by scan function: " << duration_sec.count() << " milliseconds" << endl;
	cout << "First element in output is " << output[0] << endl;
	cout << "Last element in output is " << output[n-1] << endl;
	
//	// some test functions
//	for (size_t i=0; i < n; i++) {
//		cout << arr[i] << ' ';
//	}
//	
//	cout << endl;
//	for (size_t i=0; i < n; i++) {
//		cout << output[i] << ' ';
//	}
//	cout << endl;
	
	free(arr);
	free(output);
}
