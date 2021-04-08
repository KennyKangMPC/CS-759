#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <stdio.h>
#include "optimize.h"

using namespace std;
using namespace chrono;

double time1(vec *v, data_t *dest) {	
	auto start = high_resolution_clock::now();
  	optimize1(v, dest);
  	auto end = high_resolution_clock::now();
  	double ds = duration_cast<duration<double, std::milli>>(end - start).count();
  return ds;
}

double time2(vec *v, data_t *dest) {
	auto start = high_resolution_clock::now();
	optimize2(v, dest);
	auto end = high_resolution_clock::now();
	double ds = duration_cast<duration<double, std::milli>>(end - start).count();
	return ds;
}

double time3(vec *v, data_t *dest) {
	auto start = high_resolution_clock::now();
	optimize3(v, dest);
	auto end = high_resolution_clock::now();
	double ds = duration_cast<duration<double, std::milli>>(end - start).count();
	return ds;
}

double time4(vec *v, data_t *dest) {
	auto start = high_resolution_clock::now();
	optimize4(v, dest);
	auto end = high_resolution_clock::now();
	double ds = duration_cast<duration<double, std::milli>>(end - start).count();
	return ds;
}

double time5(vec *v, data_t *dest) {
	auto start = high_resolution_clock::now();
	optimize5(v, dest);
	auto end = high_resolution_clock::now();
	double ds = duration_cast<duration<double, std::milli>>(end - start).count();
	return ds;
}

int main(int argc, char *argv[]) {
	int n = atol(argv[1]);
	data_t *data = new data_t[n];
	vec v = vec(n);
	
	for (int i = 0; i < n; i++) {
		// to make things simple, I choose to fill all 1
		data[i] = 1;
	}
	v.data = data;
	data_t dest;
	
	int iterations = 10;
	
	// next is for optimize 1
	double time = 0;
	for (int i = 0; i < iterations; i++) {
		time += time1(&v, &dest);
	}
    time /= iterations;
    cout << dest << endl;
    cout << time << endl;
    
    // next is for optimize 2
	time = 0;
	for (int i = 0; i < iterations; i++) {
		time += time2(&v, &dest);
	}
    time /= iterations;
    cout << dest << endl;
    cout << time << endl;
    
    // next is for optimize 3
	time = 0;
	for (int i = 0; i < iterations; i++) {
		time += time3(&v, &dest);
	}
    time /= iterations;
    cout << dest << endl;
    cout << time << endl;
    
    // next is for optimize 4
	time = 0;
	for (int i = 0; i < iterations; i++) {
		time += time4(&v, &dest);
	}
    time /= iterations;
    cout << dest << endl;
    cout << time << endl;
    
    // next is for optimize 5
	time = 0;
	for (int i = 0; i < iterations; i++) {
		time += time5(&v, &dest);
	}
    time /= iterations;
    cout << dest << endl;
    cout << time << endl;
}

