#ifndef NCOMP_UTIL_H
#define NCOMP_UTIL_H

#include <ncomp/types.h>

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
extern "C" {
#endif

ncomp_array* ncomp_array_alloc(void*, int, int, size_t*);
void ncomp_array_copy(ncomp_array * from, ncomp_array * to);
void ncomp_array_free(ncomp_array*, int);
size_t sizeof_ncomp_array_data(int);

#ifdef __cplusplus /* If this is a C++ compiler, end C linkage */
}
#endif

#endif
