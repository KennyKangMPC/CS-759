#include "convolution.h"
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <iostream>
#include <random>
#include <omp.h>

using namespace std;
using namespace chrono;


struct squashedMatrix {//no need this for homework but I wanna use this as a practice
	// 2D matrix stored row-wise in 1D array
	size_t height; // number of rows in matrix
	size_t width; // number of columns in matrix
	float *pMatVal;
};

int main(int argc, char *argv[]){
	size_t n = atoi(argv[1]); //image size is n*n;
	size_t m = 3; //mask size is 3*3;
	int t = atol(argv[2]); // number of threads to use
	
	// set up random value generator for both matrix
	random_device entropy_source;
	mt19937_64 generator(entropy_source()); 
	const float minIm = -10.0, maxIm = 10.0;
	const float minMa = -1.0, maxMa = 1.0;
	uniform_real_distribution<float> distIm(minIm, maxIm);
	uniform_real_distribution<float> distMa(minMa, maxMa); 
	
	struct squashedMatrix image;
	struct squashedMatrix mask;
	
	image.height = n;
	image.width = n;
	image.pMatVal = (float*)malloc(image.height * image.width * sizeof(float));
	
	mask.height = m;
	mask.width = m;
	mask.pMatVal = (float*)malloc(mask.height * mask.width * sizeof(float));
	
	//initialize image with random float value from -10.0 to 10.0
	for (size_t i = 0; i < image.height; i++) {
		for (size_t j = 0; j < image.width; j++) {
			image.pMatVal[i * image.width + j] = distIm(generator);
		}
	}
	
	//initialize mask with random float value from -1.0 to 1.0
	for (size_t i = 0; i < mask.height; i++) {
		for (size_t j = 0; j < mask.width; j++) {
			mask.pMatVal[i * mask.width + j] = distMa(generator);
		}
	}
	
	float *output = (float*)malloc(n * n * sizeof(float));
	
	omp_set_num_threads(t);
	//timing for the convolution
  	auto start = high_resolution_clock::now();
  	convolve(image.pMatVal, output, n, mask.pMatVal, m);
  	auto end = high_resolution_clock::now();
	
	auto duration_sec = duration_cast<duration<double, std::milli>>(end - start);
	cout << "First element in output is " << output[0] << endl;
	cout << "Last element in output is " << output[n*n-1] << endl;
	cout << "Time taken by scan function: " << duration_sec.count() << " milliseconds" << endl;
	//Done with the matrix then delete;
	delete[] image.pMatVal;
	delete[] mask.pMatVal;
	delete[] output;
}
