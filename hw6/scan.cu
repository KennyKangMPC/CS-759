#include "scan.cuh"
#include "stdio.h"


__host__ void scan(const float* input, float* output, unsigned int n, unsigned int threads_per_block);
