#ifndef NCOMP_INTERNAL_UTIL_H
#define NCOMP_INTERNAL_UTIL_H

#include "ncomp/util.h"
#include <ncomp/types.h>
#include <type_traits>
#include <iostream>
#include <limits.h> // for INT_MAX
#include <algorithm>
#include <vector>

void _ncomp_coerce(void *from_ptr, int from_type, void *from_missing,
                   void *to_ptr, int to_type, void *to_missing, size_t num);

template <typename T>
T* convert_to_with_copy_avoiding(void *in_arr, size_t arr_size,
       size_t arr_offset, int arr_type, NcompTypes intendedType);

void coerce_missing(int, int, const ncomp_missing *, double *, float *);

void set_subset_output_missing(void *x, size_t index_x, int type_x,
                               size_t size_x, const double &missing_x);

template <typename T>
void convert_to(void *in_arr, size_t in_arr_size, size_t in_arr_offset,
                int in_arr_type, T *out_arr) {
  switch (in_arr_type) {
  case NCOMP_DOUBLE: {
    for (size_t i = 0; i < in_arr_size; i++) {
      out_arr[i] =
          static_cast<T>(static_cast<double *>(in_arr)[i + in_arr_offset]);
    }
    break;
  }
  case NCOMP_FLOAT: {
    for (size_t i = 0; i < in_arr_size; i++) {
      out_arr[i] =
          static_cast<T>(static_cast<float *>(in_arr)[i + in_arr_offset]);
    }
    break;
  }
  case NCOMP_BOOL: {
    for (size_t i = 0; i < in_arr_size; i++) {
      out_arr[i] =
          static_cast<T>(static_cast<char *>(in_arr)[i + in_arr_offset]);
    }
    break;
  }
  case NCOMP_BYTE: {
    for (size_t i = 0; i < in_arr_size; i++) {
      out_arr[i] =
          static_cast<T>(static_cast<signed char *>(in_arr)[i + in_arr_offset]);
    }
    break;
  }
  case NCOMP_UBYTE: {
    for (size_t i = 0; i < in_arr_size; i++) {
      out_arr[i] = static_cast<T>(
          static_cast<unsigned char *>(in_arr)[i + in_arr_offset]);
    }
    break;
  }
  case NCOMP_SHORT: {
    for (size_t i = 0; i < in_arr_size; i++) {
      out_arr[i] =
          static_cast<T>(static_cast<short *>(in_arr)[i + in_arr_offset]);
    }
    break;
  }
  case NCOMP_USHORT: {
    for (size_t i = 0; i < in_arr_size; i++) {
      out_arr[i] = static_cast<T>(
          static_cast<unsigned short *>(in_arr)[i + in_arr_offset]);
    }
    break;
  }
  case NCOMP_INT: {
    for (size_t i = 0; i < in_arr_size; i++) {
      out_arr[i] =
          static_cast<T>(static_cast<int *>(in_arr)[i + in_arr_offset]);
    }
    break;
  }
  case NCOMP_UINT: {
    for (size_t i = 0; i < in_arr_size; i++) {
      out_arr[i] = static_cast<T>(
          static_cast<unsigned int *>(in_arr)[i + in_arr_offset]);
    }
    break;
  }
  case NCOMP_LONG: {
    for (size_t i = 0; i < in_arr_size; i++) {
      out_arr[i] =
          static_cast<T>(static_cast<long *>(in_arr)[i + in_arr_offset]);
    }
    break;
  }
  case NCOMP_ULONG: {
    for (size_t i = 0; i < in_arr_size; i++) {
      out_arr[i] = static_cast<T>(
          static_cast<unsigned long *>(in_arr)[i + in_arr_offset]);
    }
    break;
  }
  case NCOMP_LONGLONG: {
    for (size_t i = 0; i < in_arr_size; i++) {
      out_arr[i] =
          static_cast<T>(static_cast<long long *>(in_arr)[i + in_arr_offset]);
    }
    break;
  }
  case NCOMP_ULONGLONG: {
    for (size_t i = 0; i < in_arr_size; i++) {
      out_arr[i] = static_cast<T>(
          static_cast<unsigned long long *>(in_arr)[i + in_arr_offset]);
    }
    break;
  }
  case NCOMP_LONGDOUBLE: {
    for (size_t i = 0; i < in_arr_size; i++) {
      out_arr[i] =
          static_cast<T>(static_cast<long double *>(in_arr)[i + in_arr_offset]);
    }
    break;
  }
  default:
    break;
  }
  return;
}

/*
 * Coerce data to double, or just return a pointer to it if
 * it is already double.
 *
 * FILE:"ncl/ni/src/lib/nfp/wrapper.h"
 * IMPLEMENTATION:"./ni/src/lib/nfp/wrapper.c:9503"
 *
 * extern double coerce_input_double(void,NclBasicDataTypes,ng_size_t,int,
 * 								  NclScalar,NclScalar);
 * extern float coerce_input_float(void,NclBasicDataTypes,ng_size_t,int,
 * 								NclScalar,NclScalar);
 * extern int coerce_input_int(void,NclBasicDataTypes,ng_size_t,int,
 * 							NclScalar,NclScalar);
 * extern unsigned int coerce_input_uint(void,NclBasicDataTypes,ng_size_t,int,
 * 									  NclScalar,NclScalar);
 * extern unsigned long coerce_input_ulong(void,NclBasicDataTypes,ng_size_t,int,
 * 										NclScalar,NclScalar);
 */

template <typename T,
          typename = typename std::enable_if<
	    std::is_same<T, double>::value || std::is_same<T, float>::value ||
	    std::is_same<T, int>::value ||
	    std::is_same<T, unsigned int>::value ||
	    std::is_same<T, long>::value,
	    T>::type>
  T *coerce_input_T(void *x, int type_x, size_t size_x, int has_missing_x,
		    void *missing_x, T *missing_dx) {
  T *dx = nullptr;

  if ((std::is_same<T, double>::value && type_x == NCOMP_DOUBLE) ||
      (std::is_same<T, float>::value && type_x == NCOMP_FLOAT) ||
      (std::is_same<T, int>::value && type_x == NCOMP_INT) ||
      (std::is_same<T, unsigned int>::value && type_x == NCOMP_UINT) ||
      (std::is_same<T, unsigned long>::value && type_x == NCOMP_ULONG)) {
    dx = static_cast<T *>(x);
  } else {
    dx = new (std::nothrow) T[size_x];
    if (dx == nullptr) return nullptr;

    NcompTypes ncompType;
    if (std::is_same<T, double>::value)
      ncompType = NCOMP_DOUBLE;
    else if (std::is_same<T, float>::value)
      ncompType = NCOMP_FLOAT;
    else if (std::is_same<T, int>::value)
      ncompType = NCOMP_INT;
    else if (std::is_same<T, unsigned int>::value)
      ncompType = NCOMP_UINT;
    else if (std::is_same<T, unsigned long>::value)
      ncompType = NCOMP_ULONG;
    else
      return nullptr;

    if (has_missing_x) {
      _ncomp_coerce(x, type_x, missing_x, (void *)dx, ncompType, (void *)missing_dx, size_x);
    } else {
      _ncomp_coerce(x, type_x, nullptr, (void *)dx, ncompType, nullptr, size_x);
    }
  }

  return dx;
}

template <typename T>
void _ncomp_coerce_internal(void *from_ptr, int from_type, void *from_missing,
                            T *to_ptr, T *to_missing, bool has_missing,
                            size_t num) {
  switch (from_type) {
  case NCOMP_DOUBLE: {
    for (size_t i = 0; i < num; i++) {
      if (has_missing && ((double *)from_ptr)[i] == *((double *)from_missing)) {
        to_ptr[i] = *to_missing;
      } else {
        to_ptr[i] = static_cast<T>(static_cast<double *>(from_ptr)[i]);
      }
    }
    break;
  }
  case NCOMP_FLOAT: {
    for (size_t i = 0; i < num; i++) {
      if (has_missing && ((float *)from_ptr)[i] == *((float *)from_missing)) {
        to_ptr[i] = *to_missing;
      } else {
        to_ptr[i] = static_cast<T>(static_cast<float *>(from_ptr)[i]);
      }
    }
    break;
  }
  case NCOMP_BOOL: {
    for (size_t i = 0; i < num; i++) {
      if (has_missing && ((char *)from_ptr)[i] == *((char *)from_missing)) {
        to_ptr[i] = *to_missing;
      } else {
        to_ptr[i] = static_cast<T>(static_cast<char *>(from_ptr)[i]);
      }
    }
    break;
  }
  case NCOMP_BYTE: {
    for (size_t i = 0; i < num; i++) {
      if (has_missing &&
          ((signed char *)from_ptr)[i] == *((signed char *)from_missing)) {
        to_ptr[i] = *to_missing;
      } else {
        to_ptr[i] = static_cast<T>(static_cast<signed char *>(from_ptr)[i]);
      }
    }
    break;
  }
  case NCOMP_UBYTE: {
    for (size_t i = 0; i < num; i++) {
      if (has_missing &&
          ((unsigned char *)from_ptr)[i] == *((unsigned char *)from_missing)) {
        to_ptr[i] = *to_missing;
      } else {
        to_ptr[i] = static_cast<T>(static_cast<unsigned char *>(from_ptr)[i]);
      }
    }
    break;
  }
  case NCOMP_SHORT: {
    for (size_t i = 0; i < num; i++) {
      if (has_missing && ((short *)from_ptr)[i] == *((short *)from_missing)) {
        to_ptr[i] = *to_missing;
      } else {
        to_ptr[i] = static_cast<T>(static_cast<short *>(from_ptr)[i]);
      }
    }
    break;
  }
  case NCOMP_USHORT: {
    for (size_t i = 0; i < num; i++) {
      if (has_missing && ((unsigned short *)from_ptr)[i] ==
                             *((unsigned short *)from_missing)) {
        to_ptr[i] = *to_missing;
      } else {
        to_ptr[i] = static_cast<T>(static_cast<unsigned short *>(from_ptr)[i]);
      }
    }
    break;
  }
  case NCOMP_INT: {
    for (size_t i = 0; i < num; i++) {
      if (has_missing && ((int *)from_ptr)[i] == *((int *)from_missing)) {
        to_ptr[i] = *to_missing;
      } else {
        to_ptr[i] = static_cast<T>(static_cast<int *>(from_ptr)[i]);
      }
    }
    break;
  }
  case NCOMP_UINT: {
    for (size_t i = 0; i < num; i++) {
      if (has_missing &&
          ((unsigned int *)from_ptr)[i] == *((unsigned int *)from_missing)) {
        to_ptr[i] = *to_missing;
      } else {
        to_ptr[i] = static_cast<T>(static_cast<unsigned int *>(from_ptr)[i]);
      }
    }
    break;
  }
  case NCOMP_LONG: {
    for (size_t i = 0; i < num; i++) {
      if (has_missing && ((long *)from_ptr)[i] == *((long *)from_missing)) {
        to_ptr[i] = *to_missing;
      } else {
        to_ptr[i] = static_cast<T>(static_cast<long *>(from_ptr)[i]);
      }
    }
    break;
  }
  case NCOMP_ULONG: {
    for (size_t i = 0; i < num; i++) {
      if (has_missing &&
          ((unsigned long *)from_ptr)[i] == *((unsigned long *)from_missing)) {
        to_ptr[i] = *to_missing;
      } else {
        to_ptr[i] = static_cast<T>(static_cast<unsigned long *>(from_ptr)[i]);
      }
    }
    break;
  }
  case NCOMP_LONGLONG: {
    for (size_t i = 0; i < num; i++) {
      if (has_missing &&
          ((long long *)from_ptr)[i] == *((long long *)from_missing)) {
        to_ptr[i] = *to_missing;
      } else {
        to_ptr[i] = static_cast<T>(static_cast<long long *>(from_ptr)[i]);
      }
    }
    break;
  }
  case NCOMP_ULONGLONG: {
    for (size_t i = 0; i < num; i++) {
      if (has_missing && ((unsigned long long *)from_ptr)[i] ==
                             *((unsigned long long *)from_missing)) {
        to_ptr[i] = *to_missing;
      } else {
        to_ptr[i] =
            static_cast<T>(static_cast<unsigned long long *>(from_ptr)[i]);
      }
    }
    break;
  }
  case NCOMP_LONGDOUBLE: {
    for (size_t i = 0; i < num; i++) {
      if (has_missing &&
          ((long double *)from_ptr)[i] == *((long double *)from_missing)) {
        to_ptr[i] = *to_missing;
      } else {
        to_ptr[i] = static_cast<T>(static_cast<long double *>(from_ptr)[i]);
      }
    }
    break;
  }
  default:
    break;
  }
  return;
}

template <typename T>
T* allocateAndInit(size_t size, T initValue);

size_t prod(const size_t* shape, int ndim);

ncomp_attributes * collectAttributeList(std::vector<ncomp_single_attribute *> attrVector);
void collectAttributeList(std::vector<ncomp_single_attribute*> attrVector, ncomp_attributes * collectedAttributedList);

#endif
