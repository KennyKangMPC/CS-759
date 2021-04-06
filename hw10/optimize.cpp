#include "optimize.h"

int vec_length(vec *v) {
	return v->len;
}

data_t *get_vec_start(vec *v) {
	return v->data;
}

void optimize1(vec *v, data_t *dest) {
	int length = vec_length(v);
	data_t *d = get_vec_start(v);
	data_t temp = IDENT;
	for(int i = 0; i < length; i++) {
		temp = temp OP d[i];
    
	}
    *dest= temp;
}
void optimize2(vec *v, data_t *dest) {
	int length = vec_length(v);
	int limit = length - 1;
	data_t *d = get_vec_start(v);
	data_t x = IDENT;
	int i;
	
	// Combine 2 elements at a time
	for (i = 0; i < limit; i+= 2) {
		x = (x OP d[i]) OP d[i+1];
	}
	
	// Finish any remaining elements;
	for (; i < length; i++) {
		x = x OP d[i];
	}
	*dest = x;
}

void optimize3(vec *v, data_t *dest) {
	int length = vec_length(v);
	int limit = length - 1;
	data_t x = IDENT;
	int i;
  	/* Combine 2 elements at a time */
  	for (i = 0; i < limit; i += 2) {
    	x = x OP(d[i] OP d[i + 1]);
  	}
  	
    /* Finish any remaining elements */
  	for (; i < length; i++) {
    	x = x OP d[i];
  	}
  	*dest = x;
}

void optimize4(vec *v, data_t *dest) {
	long length = vec_length(v);
	long limit = length - 1;
	data_t *d = get_vec_start(v);
	data_t x0 = IDENT;
	data_t x1 = IDENT;
	long i;
	// Combine 2 elements at a time
	for (i = 0; i < limit; i += 2) {
    	x0 = x0 OP d[i];
    	x1 = x1 OP d[i + 1];
  	}
  	
  	// Finish any remaining elements
  	for (; i < length; i++) {
  	  x0 = x0 OP d[i];
  	}
  	*dest = x0 OP x1;
}

void optimize5(vec *v, data_t *dest) {
	long length = vec_length(v);
	long limit = length - 1;
	data_t *d = get_vec_start(v);
	data_t x0 = IDENT;
	data_t x1 = IDENT;
	data_t x2 = IDENT;
	
	long i;
	// Combine 3 elements at a time
	for (i = 0; i < limit; i += 3) {
		x0 = x0 OP d[i];
		x1 = x1 OP d[i+1];
		x2 = x2 OP d[i+2];
	}		
	
	// Finish any remaining elements
    for (; i < length; i++) {
        x0 = x0 OP d[i];
    }
    
    *dest = x0 OP x1 OP x2;
}
