#include <chrono>
#include <cstring>
#include <random>
#include <algorithm>
#include <random>
#include <stdio.h>
#include <mpi.h>

using namespace std;
using std::sort;
using std::chrono::duration;
using std::chrono::high_resolution_clock;

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	
	int n = atol(argv[1]);
	float sendBuff = new float[n];
	float recevBuff = new float[n];
	
	random_device entropy_source;
	mt19937_64 generator(entropy_source()); 
	const float min = -5.0, max = 5.0; // The range for the random number 		generator is 0.0 to n
	uniform_real_distribution<float> dist(min, max);
	
	for (int i = 0; i < n; i++) {
		sendBuff[i] = dist(generator);
		recevBuff[i] = dist(generator);
	}
	
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double t0 = 0;
	double t1 = 0;
	
	if (rank == 0) {
		auto start = high_resolution_clock::now();
		MPI_Send(sendBuff, n, MPI_FLOAT, 1, 0, MPI_COMM_WORLD);
		MPI_Recv(recevBuff, n, MPI_FLOAT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		auto end = high_resolution_clock::now();
		t0 = duration_cast<duration<double, std::milli>>(end - start).count();
		MPI_Recv(&t1, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("%f\n", t0+t1);
	} else if (rank == 1) {
		auto start = high_resolution_clock::now();
		MPI_Recv(recevBuff, n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Send(sendBuff, n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
		auto end = high_resolution_clock::now();
		t1 = duration_cast<duration<double, std::milli>>(end - start).count();
		MPI_Send(&t1, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	
	MPI_Finalize();
}


