#include <algorithm>
#include <chrono>
#include <iostream>
#include <omp.h>
#include <random>
#include <stdio.h>
#include <mpi.h>

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int n = atol(argv[1]);
    int t = atol(argv[2]);
    
    // random generator
    random_device entropy_source;
	mt19937_64 generator(entropy_source()); 
	const float min = -1.0, max = 1.0; // The range for the random number 		generator is 0.0 to n
	uniform_real_distribution<float> dist(min, max);
    
    float *arr = new float[n];
    for (int i = 0; i < n; i++) {
    	arr[i] = dist(generator)
    }
    float global_res = 0;
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Barrier(MPI_COMM_WORLD);
    
    
}

