#ifndef NCOMP_INTERNAL_UTIL_H
#define NCOMP_INTERNAL_UTIL_H

#include <ncomp/types.h>
#include <algorithm>
#include <vector>

template <typename T>
void convert_to(void *in_arr, size_t arr_size, size_t arr_offset, int arr_type,
                T *out_arr);

template <typename T>
T* convert_to_with_copy_avoiding(void *in_arr, size_t arr_size,
       size_t arr_offset, int arr_type, NcompTypes intendedType);

void coerce_missing(int, int, const ncomp_missing *, double *, float *);

template <typename T>
T *coerce_input_T(void *, int, size_t, int, void *, T *);

template <typename T>
void _ncomp_coerce_internal(void *from_ptr, int from_type, void *from_missing,
                            T *to_ptr, T *to_missing, bool has_missing,
                            size_t num);

void _ncomp_coerce(void *from_ptr, int from_type, void *from_missing,
                   void *to_ptr, int to_type, void *to_missing, size_t num);

template <typename T>
T* allocateAndInit(size_t size, T initValue);

size_t prod(const size_t* shape, int ndim);
size_t prod(const size_t* shape, int start_inclusive_idx, int end_exclusive_idx);

ncomp_attributes * collectAttributeList(std::vector<ncomp_single_attribute *> attrVector);
void collectAttributeList(std::vector<ncomp_single_attribute*> attrVector, ncomp_attributes * collectedAttributedList);

#endif
