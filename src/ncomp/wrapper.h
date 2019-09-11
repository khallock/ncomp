#ifndef NCOMP_WRAPPER_H
#define NCOMP_WRAPPER_H

#include <ncomp/types.h>

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
extern "C" {
#endif

int linint2(const ncomp_array *, const ncomp_array *, const ncomp_array *,
            const ncomp_array *, const ncomp_array *, ncomp_array *, int, int);

int eofunc(
    const ncomp_array * x_in, const int neval_in,
    const attributes & options_in,
    ncomp_array * x_out, attributes * attr_out)

#ifdef __cplusplus /* If this is a C++ compiler, end C linkage */
}
#endif

#endif
