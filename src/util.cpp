#include "ncomp/constants.h"
#include "ncomp/types.h"
#include "ncomp/util.h"
#include "ncomp_internal/util.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <type_traits>
#include <vector>
#include <iostream>

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

// If the intended type and the actual type are the same, the array is just casted
// So nothing is copied. Hence, you need to make sure that in_arr is not scoped
// out. However, if the arr_type and the intended are not the same, then a brand
// new memory locaiton is assigned and the in_arr is copied and converted to it.
// Hence, if the variable is no longer needed you would need to delete it manually.
template <typename T>
T* convert_to_with_copy_avoiding(void *in_arr, unsigned long arr_size,
        unsigned long arr_offset, int arr_type, NcompTypes intendedType) {
  if ( arr_type == intendedType) {
    // no need to copy/convert just cast it
    return static_cast<T *>(in_arr);
  } else {
    T * tmpStorage = new T[arr_size];
    convert_to<T>(in_arr, arr_size, 0, arr_type, tmpStorage);
    return tmpStorage;
  }
}

// NOTE: array_ptr is not copied, but array_shape is copied. So make sure that
//       array_ptr and the memory associated to it is not scoped out. Otherwise,
//       ncomp_array->addr is pointing to a memory location that is freed.
extern "C" ncomp_array *ncomp_array_alloc(void *array_ptr, int array_type, int array_ndim,
                               size_t *array_shape) {
  ncomp_array * new_array = (ncomp_array *) malloc(sizeof(ncomp_array));
  new_array->type = array_type;
  new_array->addr = array_ptr;
  new_array->ndim = array_ndim;
  new_array->has_missing = 0;
  new_array->msg.msg_bool = 0;

  // We are still going to copy the array_shape because some of our codes assumes
  // that this copy is happening.
  new_array->shape = (size_t *) malloc(sizeof(size_t)*array_ndim);
  for (int i = 0; i < array_ndim; i++)
    new_array->shape[i] = array_shape[i];
  return new_array;
}

extern "C" ncomp_array *ncomp_array_alloc_scalar(void *array_ptr, int array_type) {
  int array_ndim = 1;
  size_t array_shape[1] = {1};

  return ncomp_array_alloc(array_ptr, array_type, array_ndim, array_shape);
}

extern "C" void ncomp_array_free(ncomp_array* old_array, int keep_data_ptr) {
    /* free ptr by default */
    if (!keep_data_ptr)
        free(old_array->addr);

    free(old_array->shape);
    free(old_array);

    return;
}

extern "C" size_t sizeof_ncomp_array_data(int array_type) {
  size_t array_type_size = 0;
  switch (array_type) {
  case NCOMP_DOUBLE:
    array_type_size = sizeof(double);
    break;
  case NCOMP_FLOAT:
    array_type_size = sizeof(float);
    break;
  case NCOMP_BOOL:
  case NCOMP_BYTE:
  case NCOMP_UBYTE:
    array_type_size = sizeof(char);
    break;
  case NCOMP_SHORT:
  case NCOMP_USHORT:
    array_type_size = sizeof(short);
    break;
  case NCOMP_INT:
  case NCOMP_UINT:
    array_type_size = sizeof(int);
    break;
  case NCOMP_LONG:
  case NCOMP_ULONG:
    array_type_size = sizeof(long);
    break;
  case NCOMP_LONGLONG:
  case NCOMP_ULONGLONG:
    array_type_size = sizeof(long long);
    break;
  case NCOMP_LONGDOUBLE:
    array_type_size = sizeof(long double);
    break;
  }
  return array_type_size;
}

/*
 * Coerce a missing value to double.  Also, set a default missing
 * value and set a float missing value for the return (in case the
 * return type is a float).
 */
void coerce_missing(int type_x, int has_missing_x,
                    const ncomp_missing *missing_x, double *missing_dx,
                    float *missing_fx) {
  /*
   * Check for missing value and coerce if neccesary.
   */
  if (has_missing_x) {
    switch (type_x) {
    case NCOMP_DOUBLE:
      (*missing_dx) = missing_x->msg_double;
      break;
    case NCOMP_FLOAT:
      (*missing_dx) = static_cast<double>(missing_x->msg_float);
      break;
    case NCOMP_BOOL:
      (*missing_dx) = static_cast<double>(missing_x->msg_bool);
      break;
    case NCOMP_BYTE:
      (*missing_dx) = static_cast<double>(missing_x->msg_byte);
      break;
    case NCOMP_UBYTE:
      (*missing_dx) = static_cast<double>(missing_x->msg_ubyte);
      break;
    case NCOMP_SHORT:
      (*missing_dx) = static_cast<double>(missing_x->msg_short);
      break;
    case NCOMP_USHORT:
      (*missing_dx) = static_cast<double>(missing_x->msg_ushort);
      break;
    case NCOMP_INT:
      (*missing_dx) = static_cast<double>(missing_x->msg_int);
      break;
    case NCOMP_UINT:
      (*missing_dx) = static_cast<double>(missing_x->msg_uint);
      break;
    case NCOMP_LONG:
      (*missing_dx) = static_cast<double>(missing_x->msg_long);
      break;
    case NCOMP_ULONG:
      (*missing_dx) = static_cast<double>(missing_x->msg_ulong);
      break;
    case NCOMP_LONGLONG:
      (*missing_dx) = static_cast<double>(missing_x->msg_longlong);
      break;
    case NCOMP_ULONGLONG:
      (*missing_dx) = static_cast<double>(missing_x->msg_ulonglong);
      break;
    case NCOMP_LONGDOUBLE:
      (*missing_dx) = static_cast<double>(missing_x->msg_longdouble);
      break;
    default:
      break;
    }

    if (type_x != NCOMP_DOUBLE && missing_fx != nullptr) {
      switch (type_x) {
      case NCOMP_FLOAT:
        (*missing_fx) = missing_x->msg_float;
        break;
      case NCOMP_BOOL:
        (*missing_fx) = static_cast<float>(missing_x->msg_bool);
        break;
      case NCOMP_BYTE:
        (*missing_fx) = static_cast<float>(missing_x->msg_byte);
        break;
      case NCOMP_UBYTE:
        (*missing_fx) = static_cast<float>(missing_x->msg_ubyte);
        break;
      case NCOMP_SHORT:
        (*missing_fx) = static_cast<float>(missing_x->msg_short);
        break;
      case NCOMP_USHORT:
        (*missing_fx) = static_cast<float>(missing_x->msg_ushort);
        break;
      case NCOMP_INT:
        (*missing_fx) = static_cast<float>(missing_x->msg_int);
        break;
      case NCOMP_UINT:
        (*missing_fx) = static_cast<float>(missing_x->msg_uint);
        break;
      case NCOMP_LONG:
        (*missing_fx) = static_cast<float>(missing_x->msg_long);
        break;
      case NCOMP_ULONG:
        (*missing_fx) = static_cast<float>(missing_x->msg_ulong);
        break;
      case NCOMP_LONGLONG:
        (*missing_fx) = static_cast<float>(missing_x->msg_longlong);
        break;
      case NCOMP_ULONGLONG:
        (*missing_fx) = static_cast<float>(missing_x->msg_ulonglong);
        break;
      case NCOMP_LONGDOUBLE:
        (*missing_fx) = static_cast<float>(missing_x->msg_longdouble);
        break;
      }
    }
  } else { // assign mising value
    if (missing_dx != nullptr) {
      if (type_x != NCOMP_DOUBLE) {
        *missing_dx = NC_FILL_DOUBLE;
        if (missing_fx != nullptr)
          *missing_fx = NC_FILL_FLOAT;
      } else {
        *missing_dx = NC_FILL_DOUBLE;
      }
    }
  }
}

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
    std::vector<T> dxVec(size_x);
    dx = dxVec.data();

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
      _ncomp_coerce(x, type_x, missing_x, (void *)dx, ncompType,
                    (void *)missing_dx, size_x);
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

void _ncomp_coerce(void *from_ptr, int from_type, void *from_missing,
                   void *to_ptr, int to_type, void *to_missing, size_t num) {
  // find if a missing value exists
  bool has_missing = false;
  if (from_missing != nullptr && to_missing != nullptr)
    has_missing = true;

  switch (to_type) {
  case NCOMP_DOUBLE:
    _ncomp_coerce_internal<double>(from_ptr, from_type, from_missing,
                                   (double *)to_ptr, (double *)to_missing,
                                   has_missing, num);
    break;
  case NCOMP_FLOAT:
    _ncomp_coerce_internal<float>(from_ptr, from_type, from_missing,
                                  (float *)to_ptr, (float *)to_missing,
                                  has_missing, num);
    break;
  case NCOMP_BOOL:
    _ncomp_coerce_internal<char>(from_ptr, from_type, from_missing,
                                 (char *)to_ptr, (char *)to_missing,
                                 has_missing, num);
    break;
  case NCOMP_BYTE:
    _ncomp_coerce_internal<signed char>(
        from_ptr, from_type, from_missing, (signed char *)to_ptr,
        (signed char *)to_missing, has_missing, num);
    break;
  case NCOMP_UBYTE:
    _ncomp_coerce_internal<unsigned char>(
        from_ptr, from_type, from_missing, (unsigned char *)to_ptr,
        (unsigned char *)to_missing, has_missing, num);
    break;
  case NCOMP_SHORT:
    _ncomp_coerce_internal<short>(from_ptr, from_type, from_missing,
                                  (short *)to_ptr, (short *)to_missing,
                                  has_missing, num);
    break;
  case NCOMP_USHORT:
    _ncomp_coerce_internal<unsigned short>(
        from_ptr, from_type, from_missing, (unsigned short *)to_ptr,
        (unsigned short *)to_missing, has_missing, num);
    break;
  case NCOMP_INT:
    _ncomp_coerce_internal<int>(from_ptr, from_type, from_missing,
                                (int *)to_ptr, (int *)to_missing, has_missing,
                                num);
    break;
  case NCOMP_UINT:
    _ncomp_coerce_internal<unsigned int>(
        from_ptr, from_type, from_missing, (unsigned int *)to_ptr,
        (unsigned int *)to_missing, has_missing, num);
    break;
  case NCOMP_LONG:
    _ncomp_coerce_internal<long>(from_ptr, from_type, from_missing,
                                 (long *)to_ptr, (long *)to_missing,
                                 has_missing, num);
    break;
  case NCOMP_ULONG:
    _ncomp_coerce_internal<unsigned long>(
        from_ptr, from_type, from_missing, (unsigned long *)to_ptr,
        (unsigned long *)to_missing, has_missing, num);
    break;
  case NCOMP_LONGLONG:
    _ncomp_coerce_internal<long long>(
        from_ptr, from_type, from_missing, (long long *)to_ptr,
        (long long *)to_missing, has_missing, num);
    break;
  case NCOMP_ULONGLONG:
    _ncomp_coerce_internal<unsigned long long>(
        from_ptr, from_type, from_missing, (unsigned long long *)to_ptr,
        (unsigned long long *)to_missing, has_missing, num);
    break;
  case NCOMP_LONGDOUBLE:
    _ncomp_coerce_internal<long double>(
        from_ptr, from_type, from_missing, (long double *)to_ptr,
        (long double *)to_missing, has_missing, num);
    break;
  }
  return;
}

// Allocates a variable of a given size and initializes all the elements
template <typename T>
T* allocateAndInit(unsigned long size, T initValue) {
  T * tmpVar = new T[size];
  std::fill_n(tmpVar, size, initValue);
  return tmpVar;
}

extern "C" ncomp_single_attribute * create_ncomp_single_attribute_from_ncomp_array(
  char * name,
  ncomp_array * value)
{
  int size_name = strlen(name) + 1;
  char * copy_of_name = new char[size_name];
  std::copy(name, name+size_name, copy_of_name);

  ncomp_single_attribute * out_ncomp_single_attribute = new ncomp_single_attribute[1];
  out_ncomp_single_attribute->name = copy_of_name;
  out_ncomp_single_attribute->value = value;
  return out_ncomp_single_attribute;
}

// Creates a ncomp_single_attribute. Note that the data is not copied. So, make
// sure it is not scoped out if you need to use it later.
// Note that the name is copied. So, no need to manually copy it.
extern "C" ncomp_single_attribute * create_ncomp_single_attribute(
  char * name,
  void * data_ptr,
  int data_type,
  int data_ndim,
  size_t * data_shape)
{
  ncomp_array * value = ncomp_array_alloc(data_ptr, data_type, data_ndim, data_shape);
  return create_ncomp_single_attribute_from_ncomp_array(name, value);
}

extern "C" ncomp_single_attribute* create_ncomp_single_attribute_from_1DArray(
  char * name,
  void * data_ptr,
  int data_type,
  size_t nelem)
{
    int data_ndim = 1;
    size_t data_shape[1] = {nelem};

    return create_ncomp_single_attribute(name, data_ptr, data_type, data_ndim, data_shape);
}

extern "C" ncomp_single_attribute* create_ncomp_single_attribute_from_scalar(
  char * name,
  void * data_ptr,
  int data_type)
{
    return create_ncomp_single_attribute_from_1DArray(name, data_ptr, data_type, 1);
}

extern "C" ncomp_attributes* ncomp_attributes_allocate(int nAttribute) {
  ncomp_attributes * output = new ncomp_attributes;
  output->nAttribute = nAttribute;
  output->attribute_array = new ncomp_single_attribute*[nAttribute];
  return output;
}

ncomp_attributes * collectAttributeList(std::vector<ncomp_single_attribute *> attrVector) {
  ncomp_attributes * output_attribute_list = ncomp_attributes_allocate(attrVector.size());
  for (int i = 0; i < attrVector.size(); ++i) {
    output_attribute_list->attribute_array[i] = attrVector[i];
  }
  return output_attribute_list;
}

void collectAttributeList(std::vector<ncomp_single_attribute*> attrVector, ncomp_attributes * collectedAttributedList) {
  collectedAttributedList->nAttribute = attrVector.size();
  collectedAttributedList->attribute_array  = new ncomp_single_attribute*[attrVector.size()];
  for (int i = 0; i < attrVector.size(); ++i) {
    collectedAttributedList->attribute_array[i] = attrVector[i];
  }
}

// Searches for a given attribute. If it exists it returns 1 and sets the
// attributePosInList to the proper position in the array. Otherwise it returns
// zero (0) and attributePosInList is unchanged.
extern "C" int hasAttribute(
  const ncomp_attributes* attributeList,
  const char* attributeName,
  int* attributePosInList)
{
  for (int i=0; i < attributeList->nAttribute; ++i) {
    if (strcmp(attributeList->attribute_array[i]->name, attributeName) == 0) {
      *attributePosInList = i;
      return 1;
    }
  }
  return 0;
}

// Searches for an attribute and returns it. If the attribute doesn't exists, returns the default value provided by the user.
extern "C" ncomp_single_attribute* getAttributeOrDefault_ncomp_single_attribute(
  const ncomp_attributes* attributeList,
  const char* attributeName,
  const ncomp_single_attribute* defaultValue)
{
  int attributePosInList = -1;
  if (hasAttribute(attributeList, attributeName, &attributePosInList)==1) {
    return attributeList->attribute_array[attributePosInList];
  } else {
    return (ncomp_single_attribute*) defaultValue;
  }
}

// Searches for an attribute and returns it; If the attribute doesn't exists, it returns nullptr
extern "C" ncomp_single_attribute* getAttribute(
  const ncomp_attributes * attributeList,
  const char* attributeName)
{
  ncomp_single_attribute* defaultValue = nullptr;
  return getAttributeOrDefault_ncomp_single_attribute(attributeList,attributeName, defaultValue);
}

// a variant of getAttributeOrDefault with easier way of providing defaultValue.
extern "C" void* getAttributeOrDefault(
  const ncomp_attributes * attributeList,
  const char* attributeName,
  const void * defaultValue)
{
  ncomp_single_attribute * tmpOutput = getAttribute(attributeList, attributeName);
  if (tmpOutput != nullptr) {
    return tmpOutput->value->addr;  // It is assumed that addr could be casted to type T.
                                    // May be we should add a cast check here to make sure
                                    // we could indeed cast.
  } else {
    return (void *) defaultValue;
  }
}

//
size_t prod(const size_t* shape, int ndim) {
  size_t nElements = 1;
  for (int i = 0; i<ndim; ++i) {
    nElements *= *(shape+i);
  }
  return nElements;
}

void print_ncomp_array(const char * name, const ncomp_array * in) {
  std::cout << name << " type: " << in->type << "\n";

  std::cout << name << " ndim: " << in->ndim << "\n";

  std::cout << name << " dims: [ ";

  int nelem = 1;
  for (int i = 0; i < in->ndim; ++i) {
    std::cout << in->shape[i] << " ";
    nelem *= in->shape[i];
  }
  std::cout << "]\n";

  std::cout << name <<" has_missing: " << in->has_missing << "\n";

  switch(in->type) {
    case NCOMP_DOUBLE:
      std::cout << name << " msg: " << in->msg.msg_double << "\n";

      std::cout << name << " data: [ ";
      for (int i = 0; i < nelem; ++i) {
        std::cout << static_cast<double *>(in->addr)[i] << " ";
      }
      std::cout << "]\n";

      break;
    case NCOMP_FLOAT:
      std::cout << name << " msg: " << in->msg.msg_float << "\n";

      std::cout << name << " data: [ ";
      for (int i = 0; i < nelem; ++i) {
        ;
        std::cout << static_cast<float*>(in->addr)[i] << " ";
      }
      std::cout << "]\n";
      break;
    case NCOMP_INT:

      std::cout << name << " msg: " << in->msg.msg_int << "\n";

      std::cout << name << " data: [ ";
      for (int i = 0; i < nelem; ++i) {
        ;
        std::cout << static_cast<int *>(in->addr)[i] << " ";
      }
      std::cout<< "]\n";
      break;
    case NCOMP_CHAR:
      std::cout << name << " data: [ " << static_cast<char *>(in->addr) << " ]\n";
      break;
    default:
      std::cout << "ERROR: TYPE NOT FOUND OR NOT IMPLEMENTED YET\n";
  }
}

void print_ncomp_attributes(const ncomp_attributes * in) {
  std::cout << "nAttribtues: " << in->nAttribute << "\n";
  for (int i = 0; i < in->nAttribute; ++i) {
    char * attr_name = in->attribute_array[i]->name;
    std::cout << "attr Name:  "<< attr_name << "\n";
    print_ncomp_array(attr_name, in->attribute_array[i]->value);
    std::cout << std::endl;
  }
}

// explicit function instantiations
template void convert_to<double>(void *, size_t, size_t, int, double *);
template void convert_to<float>(void *, size_t, size_t, int, float *);
template double * allocateAndInit(size_t, double);
template float * allocateAndInit(size_t, float);
template double * convert_to_with_copy_avoiding(void *, size_t, size_t, int, NcompTypes);
