#include <algorithm>
#include <chrono>
#include <iostream>
#include <omp.h>
#include <random>
#include <stdio.h>
#include <mpi.h>

#include "reduce.h"

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int n = atol(argv[1]);
    int t = atol(argv[2]);
    omp_set_num_threads(t);
    // random generator
    random_device entropy_source;
	mt19937_64 generator(entropy_source()); 
	const float min = -1.0, max = 1.0; // The range for the random number 		generator is 0.0 to n
	uniform_real_distribution<float> dist(min, max);
    
    float *arr = new float[n];
    for (int i = 0; i < n; i++) {
    	arr[i] = dist(generator);
    }
    
    float global_res = 0;
    int rank;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    reduce(arr, 0, n);
    MPI_Barrier(MPI_COMM_WORLD);
    // start timing
    auto start = high_resolution_clock::now();
    float result = reduce(arr, 0, n);
    MPI_Reduce(&result, &global_res, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    auto end = high_resolution_clock::now();
    double ms = duration_cast<duration<double, std::milli>>(end - start).count();
    
    if (rank == 0) {
    	printf("%f\n%f\n", global_res, ms);
    }
    
    MPI_Finalize();
}
