#ifndef NCOMP_UTIL_H
#define NCOMP_UTIL_H

#include <ncomp/types.h>

ncomp_array* ncomp_array_alloc(void*, int, int, size_t*);
ncomp_array* ncomp_array_alloc_with_data(int, int, size_t*);
void         ncomp_array_free(ncomp_array*, int);
int sizeof_ncomp_array_data(int);

void _to_double(void*, size_t, size_t, int, double*);
void _to_float(void*, size_t, size_t, int, float*);

void coerce_missing(int, int, ncomp_missing*, double*, float*);
double *coerce_input_double(void*, int, size_t, int, void*, double*);
void _ncomp_coerce(void*, int, void*, void*, int, void*, size_t);
void _ncomp_coerce_to_double(void*, int, void*, double*, double*, size_t);
void _ncomp_coerce_to_float(void*, int, void*, float*, float*, size_t);

void coerce_subset_input_double_step(void*, double*, size_t, size_t, int, size_t, int, void*, double*);
void coerce_subset_input_float_step(void*, float*, size_t, size_t, int, size_t, int, void*, float*);
void set_subset_output_missing_step(void*, size_t, size_t, int, size_t, double);
void coerce_output_float_or_double_step(void*, double*, int, size_t, size_t, size_t);
#endif
