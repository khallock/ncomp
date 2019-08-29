#ifndef NCOMP_UTIL_H
#define NCOMP_UTIL_H

#include <ncomp/types.h>

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
extern "C" {
#endif

ncomp_array* ncomp_array_alloc(void*, int, int, size_t*);
void ncomp_array_free(ncomp_array*, int);
size_t sizeof_ncomp_array_data(int);
void set_subset_output_missing(void *x, size_t index_x, int type_x,
                               size_t size_x, const double &missing_x);

#ifdef __cplusplus /* If this is a C++ compiler, end C linkage */
}
#endif

#endif
