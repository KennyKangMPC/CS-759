#include "convolution.h"
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <iostream>
#include <random>

using namespace std;
using namespace chrono;


struct squashedMatrix {//no need this for homework but I wanna use this as a practice
	// 2D matrix stored row-wise in 1D array
	size_t height; // number of rows in matrix
	size_t width; // number of columns in matrix
	float *pMatVal;
};

int main(int argc, char *argv[]){
	size_t n = atoi(argv[1]);
	size_t m = atoi(argv[2]);
	
	// set up random value generator for both matrix
	random_device entropy_source;
	mt19937_64 generator(entropy_source()); 
	const float minIm = -10.0, maxIm = 10.0;
	const float minMa = -1.0, maxMa = 1.0;
	// there are tons of oter distributino that could be found from https://en.cppreference.com/w/cpp/header/random
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
	
	//timing for the convolution
  	auto start = high_resolution_clock::now();
  	convolve(image.pMatVal, output, n, mask.pMatVal, m);
  	auto end = high_resolution_clock::now();
	
	//print out image
//	for (size_t i = 0; i < n; i++) {
//		for (size_t j = 0; j < n; j++) {
//			cout << image.pMatVal[i * n + j] << "\t"; // print out in 2D format
//		}
//		cout << endl;
//	}
	
	//print out mask
//	for (size_t i = 0; i < m; i++) {
//		for (size_t j = 0; j < m; j++) {
//			cout << mask.pMatVal[i * m + j] << "\t"; // print out in 2D format
//		}
//		cout << endl;
//	}
	
	//print out output
//	for (size_t i = 0; i < n; i++) {
//		for (size_t j = 0; j < n; j++) {
//			cout << output[i * n + j] << "\t"; // print out in 2D format
//		}
//		cout << endl;
//	}
	
	auto duration_sec = duration_cast<duration<double, std::milli>>(end - start);
	cout << "Time taken by scan function: " << duration_sec.count() << " milliseconds" << endl;
	cout << "First element in output is " << output[0] << endl;
	cout << "Last element in output is " << output[n-1] << endl;
	
	//Done with the matrix then delete;
	delete[] image.pMatVal;
	delete[] mask.pMatVal;
	delete[] output;
}
