#ifndef NCOMP_WRAPPER_H
#define NCOMP_WRAPPER_H

#include <ncomp/types.h>

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
extern "C" {
#endif

int linint2(const ncomp_array *, const ncomp_array *, const ncomp_array *,
            const ncomp_array *, const ncomp_array *, ncomp_array *, int, int);

int moc_globe_atl( const ncomp_array *, const ncomp_array *, const ncomp_array *,
                   const ncomp_array *, const ncomp_array *, const ncomp_array *,
                   ncomp_array * )

#ifdef __cplusplus /* If this is a C++ compiler, end C linkage */
}
#endif

#endif
