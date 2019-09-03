#ifndef NCOMP_INTERNAL_UTIL_H
#define NCOMP_INTERNAL_UTIL_H

#include <ncomp/types.h>

template <typename T>
void convert_to(void *in_arr, size_t arr_size, size_t arr_offset, int arr_type,
                T *out_arr);

void coerce_missing(int, int, const ncomp_missing *, double *, float *);

template <typename T>
T *coerce_input_T(void *, int, size_t, int, void *, T *);

template <typename T>
void _ncomp_coerce_internal(void *from_ptr, int from_type, void *from_missing,
                            T *to_ptr, T *to_missing, bool has_missing,
                            size_t num);

void _ncomp_coerce(void *from_ptr, int from_type, void *from_missing,
                   void *to_ptr, int to_type, void *to_missing, size_t num);

int hasAttribute(const attributes& attributeList, const char* attributeName, int& attributePosInList);
void getAttributeOrDefault(const attributes& attributeList, const char* attributeName, const single_attribute* defaultValue, single_attribute* output);
void getAttribute(const attributes& attributeList, const char* attributeName, single_attribute* output);

#endif
