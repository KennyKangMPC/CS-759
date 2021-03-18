#include <iostream>
#include <omp.h>
#include <stdio.h>

int factorial(int n) {

	int result = 1;
	//#pragma omp parallel for reduction(*:result)
	for (int i = 2; i <= n; i++) {
		result *= i;
	}
	return result;
}

int main() {

	const int nThreads = 4;
	omp_set_num_threads(nThreads);
	
	printf("Number of threads: %d\n", nThreads);
#pragma omp parallel
	{	
		int myId = omp_get_thread_num();
		
		printf("I am thread No.: %d\n", myId);
	}
	
#pragma omp parallel for
	for (int i = 1; i <= 8; i++) {
		printf("%d!=%d\n", i, factorial(i));
	}
	return 0;
}
