#include "ncomp/constants.h"
#include "ncomp/types.h"
#include "ncomp/util.h"
#include "ncomp_internal/util.hpp"
#include <stdlib.h>
#include <vector>


extern "C" ncomp_array *ncomp_array_alloc(void *array_ptr, int array_type, int array_ndim,
                               size_t *array_shape) {
  ncomp_array *new_array = (ncomp_array *)malloc(
      sizeof(ncomp_array) + (array_ndim - 1) * sizeof(size_t));
  new_array->type = array_type;
  new_array->addr = array_ptr;
  new_array->ndim = array_ndim;
  new_array->has_missing = 0;
  new_array->msg.msg_bool = 0;
  int i;
  for (i = 0; i < array_ndim; i++)
    new_array->shape[i] = array_shape[i];
  return new_array;

  // std::shared_ptr<size_t> local_array_shape(new size_t[array_ndim],
  //                                           std::default_delete<size_t[]>());

  // std::unique_ptr<ncomp_array> ncompStruct(new ncomp_array);
  // ncompStruct->type = array_type;
  // ncompStruct->addr = array_ptr;
  // ncompStruct->ndim = array_ndim;
  // ncompStruct->has_missing = 0;
  // ncompStruct->msg.msg_bool = 0;
  // ncompStruct->shape = local_array_shape.get();

  // return ncompStruct.get();
}

extern "C" void ncomp_array_free(ncomp_array* old_array, int keep_data_ptr) {
    /* free ptr by default */
    if (!keep_data_ptr)
        free(old_array->addr);

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
  // Check for missing value and coerce if neccesary.
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

    if (type_x != NCOMP_DOUBLE) {
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
  } else { // assign DEFAULT mising value just in case
    if (type_x != NCOMP_DOUBLE) {
      *missing_dx = static_cast<double>(NC_FILL_FLOAT);
      *missing_fx = NC_FILL_FLOAT;
    } else {
      *missing_dx = NC_FILL_DOUBLE;
    }
  }
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

// Sets a subset of the output data to missing.
void set_subset_output_missing(void *x, size_t index_x, int type_x,
                               size_t size_x, const double &missing_x) {
  if (type_x != NCOMP_DOUBLE) {
    for (auto i = 0; i < size_x; i++) {
      (static_cast<float *>(x))[index_x + i] = static_cast<float>(missing_x);
    }
  } else {
    for (auto i = 0; i < size_x; i++) {
      (static_cast<double *>(x))[index_x + i] = missing_x;
    }
  }
}
