#include "scan.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <>
// I searched how to do timing from https://www.geeksforgeeks.org/measure-execution-time-function-cpp/

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[]){
	size_t n = atol(argv[1]);
	float *arr = (float*)malloc(n * sizeof(float));
	float *output = (float*)malloc(n * sizeof(float));
	
	for (int t = 0; i < n; i++) {
		arr[i] = 2 * ((float) rand()) / RAND_MAX - 1;
	}
	
	auto start = high_resolution_clock::now();
	scan(arr, output, n);
	auto stop = high_resolution_clock::now();
	
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << "Time taken by scan function: " << duration.count() << " milliseconds" << endl;
	cout << "First element in output is " << output[0] << endl;
	cout << "Last element in output is " << output[n-1] << endl;
	free(arr);
	free(output);
}
