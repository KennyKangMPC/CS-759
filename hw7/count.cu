#include <thrust/device_vector.h>
#include "count.cuh"

void count(const thrust::device_vector<int> &d_in,
           thrust::device_vector<int> &values,
           thrust::device_vector<int> &counts) {
	
	thrust::device_vector<int> dVal(int(d_in.size()));
    thrust::fill(dVal.begin(), dVal.end(), 1);
	values = d_in;
	thrust::sort(values.begin(), values.end());
	auto back = thrust::reduce_by_key(values.begin(), values.end(), dVal.begin(), values.begin(), counts.begin());
	values.resize(back.first - values.begin());
	counts.resize(back.second - counts.begin());
}
