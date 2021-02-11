#include "matmul.h"
#include <iostream>
#include <chrono>
#include <cstddef>
#include <vector>
#include <random>

using namespace std;
using namespace chrono;

int main(int argc, char *argv[]) {
	// set up random value generator for both matrix
	random_device entropy_source;
	mt19937_64 generator(entropy_source()); 
	const double min = -1.0, max = 1.0; // The range for the random number generator is -1.0 to 1.0
	// there are tons of oter distributino that could be found from https://en.cppreference.com/w/cpp/header/random
	uniform_real_distribution<double> dist(min, max);
	
	size_t n = 1024;
	
	double *A = (double*)malloc(n * n *sizeof(double));
	double *B = (double*)malloc(n * n * sizeof(double));
	auto result = new double[n * n]();
	// initialize with random number
	for (size_t i = 0; i < n * n; i++) {
		A[i] = dist(generator);
		B[i] = dist(generator);
	}
	
	vector<double> Av, Bv; // initialize the vector
	for (size_t i = 0; i < n * n; i++) {
		Av.push_back(A[i]);
		Bv.push_back(B[i]);
	}
	
	// I know there is a way that I don't have to rewrite below 4 times.
	cout << n << endl;
	//test 1
	auto start = high_resolution_clock::now();
  	mmul1(A, B, result, n);
  	auto end = high_resolution_clock::now();

 	auto duration_sec = duration_cast<duration<double, std::milli>>(end - start);
  	cout << duration_sec.count() << endl;
 	cout << result[n * n - 1] << endl;
	
	//test 2
	delete[] result;
	result = new double[n * n]();
	start = high_resolution_clock::now();
  	mmul2(A, B, result, n);
  	end = high_resolution_clock::now();

 	duration_sec = duration_cast<duration<double, std::milli>>(end - start);
  	cout << duration_sec.count() << endl;
 	cout << result[n * n - 1] << endl;
	
	//test 3
	delete[] result;
	result = new double[n * n]();
	start = high_resolution_clock::now();
  	mmul3(A, B, result, n);
  	end = high_resolution_clock::now();

 	duration_sec = duration_cast<duration<double, std::milli>>(end - start);
  	cout << duration_sec.count() << endl;
 	cout << result[n * n - 1] << endl;
	
	//test 4
	delete[] result;
	result = new double[n * n]();
	start = high_resolution_clock::now();
  	mmul4(Av, Bv, result, n);
  	end = high_resolution_clock::now();

 	duration_sec = duration_cast<duration<double, std::milli>>(end - start);
  	cout << duration_sec.count() << endl;
 	cout << result[n * n - 1] << endl;
 	
 	//delete
 	delete[] A;
	delete[] B;
	delete[] result;
}
