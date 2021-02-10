#include "convolution.h"
#include <iostream>

using namespace std;

int main(int argc, char *argv[]){
	// test cases from piazza. Credit to Katherine Heidi Fehr
//	const size_t = 4;
//	const size_t = 3;
//	float image[n * n]{1.2, 3.4, 5.6, 7.8, 9.1, 2.3, 4.5, 6.7, 8.9,1.2, 3.4, 5.6, 7.8, 9.1, 2.3, 4.5};;
//	float mask[m * m]{0, 0, 1, 0, 1, 0, 1, 0, 0};
//	float output[n * n];
	
	const size_t = 4;
	const size_t = 3;
	float image[n * n]{1, 3, 4, 8, 6, 5, 2, 4, 3, 4, 6, 8, 1, 4, 5, 2};
	float mask[m * m]{0, 0, 1, 0, 1, 0, 1, 0, 0};
	float output[n * n];
	
	convolve(image, output, n, mask, m);
	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << output[i * n + j] << "\t"; // print out in 2D format
		}
		cout << endl;
	}
}

