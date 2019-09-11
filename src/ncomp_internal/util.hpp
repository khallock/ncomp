#ifndef NCOMP_INTERNAL_UTIL_H
#define NCOMP_INTERNAL_UTIL_H

#include <ncomp/types.h>
#include <algorithm>
#include <vector>

template <typename T>
void convert_to(void *in_arr, size_t arr_size, size_t arr_offset, int arr_type,
                T *out_arr);

template <typename T>
T* convert_to_with_copy_avoiding(const void *in_arr, size_t arr_size,
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

int hasAttribute(
  const ncomp_attributes * attributeList,
  const char* attributeName,
  int& attributePosInList);

int getAttributeOrDefault(
  const ncomp_attributes * attributeList,
  const char* attributeName,
  const ncomp_single_attribute* defaultValue,
  ncomp_single_attribute* output);

int getAttribute(
  const ncomp_attributes * attributeList,
  const char* attributeName,
  ncomp_single_attribute* output);

void* getAttributeOrDefault(
  const ncomp_attributes * attributeList,
  const char* attributeName,
  const void * defaultValue);

size_t prod(const size_t* shape, int ndim);

ncomp_single_attribute* create_ncomp_single_attribute(char * name, void * data, NcompTypes type, int ndim, size_t * dims);

ncomp_attributes * collectAttributeList(std::vector<ncomp_single_attribute *> attrVector);

#endif
