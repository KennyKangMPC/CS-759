#include "scan.cuh"
#include "stdio.h"

__global__ void hs_scanl(float *g_odata. float *g_idata, int n){// another possible version
	extern volatile __shared__  float temp[]; // allocated on invocation
	int thid = threadIdx.x;
	int thlen = blockDim.x;
	int blid = blockIdx.x;
	int id = thlen * blid + thid;
    int pout = 0, pin = 1;
    
    if (id >= n) {
    	return;
    }
    temp[thid] = g_idata[id]; // this is for inclusive scan. In lecture it is for exclusive
    __syncthreads();
    
    for( int offset = 1; offset < thlen; offset *= 2 ) {
        pout = 1 - pout; // swap double buffer indices
        pin  = 1 - pout;

        if (thid >= offset)
            temp[pout * thlen + thid] = temp[pin * thlen + thid] + temp[pin * thlen + thid - offset];
        else
            temp[pout * thlen + thid] = temp[pin * thlen + thid];

        __syncthreads(); // I need this here before I start next iteration
    }
    
    if (pout *thlen + thid < thlen)
        g_odata[id] = temp[pout * n + thid];
}

// here the code could also handle more blocks. In lecture note, only 1
__global__ void hs_scan(float *g_od, float *g_id, float *g_bs, int n, int isNull) {
	extern volatile __shared__  float temp[]; // allocated on invocation
	int thid = threadIdx.x;
	int thlen = blockDim.x;
	int blid = blockIdx.x;
	int id = thlen * blid + thid;
	int pout = 0, pin = 1;
	
	if (id < n) {
		temp[thid] = g_id[id];
		for( int offset = 1; offset < n; offset *= 2 ) { // modified code from lecture notes
        	pout = 1 - pout; // swap double buffer indices
        	pin  = 1 - pout;
        	
        	temp[pout * thlen + thid] = temp[pin * thlen + thid];
        	if (thid >= offset) {
        		temp[pout * thlen + thid] += temp[pin * thlen + thid - offset];
        	}
			
        	__syncthreads(); // I need this here before I start next iteration 
   		}
   		int t_id = pout * thlen + thid;
   		if (pout *thlen + thid < thlen) {
   			g_od[id] = temp[pout * n + thid];
   		}
   		
   		if (isNull == 0 && thid == thlen - 1) { // add block
   			g_bs[blid] = temp[pout * n + thid];
   		}
	}
}

__global__ void inAdd(float *g_od, float *g_os, int n) {
	int thid = threadIdx.x;
	int thlen = blockDim.x;
	int blid = blockIdx.x;
	int id = thlen * blid + thid;
	if (id < n && blid != 0) { //inclusive add
		g_od[id] += g_os[blid - 1];
	}
}

// "inclusive scan". Use lecture notes
__host__ void scan(const float* input, float* output, unsigned int n, unsigned int threads_per_block) {
	// allocate cuda memory
	float *g_id, *g_od
	cudaMalloc(&g_id, n * sizeof(float));
  	cudaMalloc(&g_od, n * sizeof(float));
  	
  	// allocate sum cuda memory for each iteration and overall sum
  	// this is needed due to various number of blocks
  	float *g_is, *g_os, *g_em
  	// need to calculate number of block
	int num_block = (n - 1 + threads_per_block) / threads_per_block;
	cudaMalloc(&g_is, num_block * sizeof(float));
	cudaMalloc(&g_os, num_block * sizeof(float));
	cudaMalloc(&g_em, num_block * sizeof(float));
	
	// map to device input  for g_id
	cudaMemcpy(g_id, input, n * sizeof(float), cudaMemcpyHostToDevice);
	
	// number of memory for shared size
	int size_shareM = 2 * threads_per_block * sizeof(float);
	hs_scan<<<num_block, threads_per_block, size_shareM>>>(g_od, g_id, g_is, n, 0);
	hs_scan<<<1, threads_per_block, size_shareM>>>(g_os, g_is, g_em, num_block, 1);
	inAdd<<<num_block, threads_per_block>>>(g_od, g_os, n);
	
	// copy back from device to host to host
	cudaMemcpy(out, g_od, n * sizeof(float), cudaMemcpyDeviceToHosts);
	cudaDeviceSynchronize();
	
	// free array
	cudaFree(g_id);
  	cudaFree(g_od);
  	cudaFree(g_is);
  	cudaFree(g_os);
  	cudaFree(g_em);
}
