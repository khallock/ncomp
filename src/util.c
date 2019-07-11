#include <stdlib.h>
#include "ncomp.h"
#include "util.h"

#define _ncomp_coerce_to_type(a, b, c, d, e, f) _ncomp_coerce_to_bool(a, b, c, d, e, f)
#define C_TYPE char
#define NCOMP_TYPE NCOMP_BOOL
#include "coerce.def"
#undef _ncomp_coerce_to_type
#undef C_TYPE
#undef NCOMP_TYPE

#define _ncomp_coerce_to_type(a, b, c, d, e, f) _ncomp_coerce_to_byte(a, b, c, d, e, f)
#define C_TYPE signed char
#define NCOMP_TYPE NCOMP_BYTE
#include "coerce.def"
#undef _ncomp_coerce_to_type
#undef C_TYPE
#undef NCOMP_TYPE

#define _ncomp_coerce_to_type(a, b, c, d, e, f) _ncomp_coerce_to_ubyte(a, b, c, d, e, f)
#define C_TYPE unsigned char
#define NCOMP_TYPE NCOMP_UBYTE
#include "coerce.def"
#undef _ncomp_coerce_to_type
#undef C_TYPE
#undef NCOMP_TYPE

#define _ncomp_coerce_to_type(a, b, c, d, e, f) _ncomp_coerce_to_short(a, b, c, d, e, f)
#define C_TYPE short
#define NCOMP_TYPE NCOMP_SHORT
#include "coerce.def"
#undef _ncomp_coerce_to_type
#undef C_TYPE
#undef NCOMP_TYPE

#define _ncomp_coerce_to_type(a, b, c, d, e, f) _ncomp_coerce_to_ushort(a, b, c, d, e, f)
#define C_TYPE unsigned short
#define NCOMP_TYPE NCOMP_USHORT
#include "coerce.def"
#undef _ncomp_coerce_to_type
#undef C_TYPE
#undef NCOMP_TYPE

#define _ncomp_coerce_to_type(a, b, c, d, e, f) _ncomp_coerce_to_int(a, b, c, d, e, f)
#define C_TYPE int
#define NCOMP_TYPE NCOMP_INT
#include "coerce.def"
#undef _ncomp_coerce_to_type
#undef C_TYPE
#undef NCOMP_TYPE

#define _ncomp_coerce_to_type(a, b, c, d, e, f) _ncomp_coerce_to_uint(a, b, c, d, e, f)
#define C_TYPE unsigned int
#define NCOMP_TYPE NCOMP_UINT
#include "coerce.def"
#undef _ncomp_coerce_to_type
#undef C_TYPE
#undef NCOMP_TYPE

#define _ncomp_coerce_to_type(a, b, c, d, e, f) _ncomp_coerce_to_long(a, b, c, d, e, f)
#define C_TYPE long
#define NCOMP_TYPE NCOMP_LONG
#include "coerce.def"
#undef _ncomp_coerce_to_type
#undef C_TYPE
#undef NCOMP_TYPE

#define _ncomp_coerce_to_type(a, b, c, d, e, f) _ncomp_coerce_to_ulong(a, b, c, d, e, f)
#define C_TYPE unsigned long
#define NCOMP_TYPE NCOMP_ULONG
#include "coerce.def"
#undef _ncomp_coerce_to_type
#undef C_TYPE
#undef NCOMP_TYPE

#define _ncomp_coerce_to_type(a, b, c, d, e, f) _ncomp_coerce_to_longlong(a, b, c, d, e, f)
#define C_TYPE long long
#define NCOMP_TYPE NCOMP_LONGLONG
#include "coerce.def"
#undef _ncomp_coerce_to_type
#undef C_TYPE
#undef NCOMP_TYPE

#define _ncomp_coerce_to_type(a, b, c, d, e, f) _ncomp_coerce_to_ulonglong(a, b, c, d, e, f)
#define C_TYPE unsigned long long
#define NCOMP_TYPE NCOMP_ULONGLONG
#include "coerce.def"
#undef _ncomp_coerce_to_type
#undef C_TYPE
#undef NCOMP_TYPE

#define _ncomp_coerce_to_type(a, b, c, d, e, f) _ncomp_coerce_to_float(a, b, c, d, e, f)
#define C_TYPE float
#define NCOMP_TYPE NCOMP_FLOAT
#include "coerce.def"
#undef _ncomp_coerce_to_type
#undef C_TYPE
#undef NCOMP_TYPE

#define _ncomp_coerce_to_type(a, b, c, d, e, f) _ncomp_coerce_to_double(a, b, c, d, e, f)
#define C_TYPE double
#define NCOMP_TYPE NCOMP_DOUBLE
#include "coerce.def"
#undef _ncomp_coerce_to_type
#undef C_TYPE
#undef NCOMP_TYPE

#define _ncomp_coerce_to_type(a, b, c, d, e, f) _ncomp_coerce_to_longdouble(a, b, c, d, e, f)
#define C_TYPE long double
#define NCOMP_TYPE NCOMP_LONGDOUBLE
#include "coerce.def"
#undef _ncomp_coerce_to_type
#undef C_TYPE
#undef NCOMP_TYPE

void _to_double(void* arr, size_t arr_index, size_t arr_size, int arr_type, double* out_arr) {
    long i;
    switch(arr_type) {
        case NCOMP_DOUBLE:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (double)((double*)arr)[i + arr_index];
            break;
        case NCOMP_FLOAT:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (double)((float*)arr)[i + arr_index];
            break;
        case NCOMP_BOOL:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (double)((char*)arr)[i + arr_index];
            break;
        case NCOMP_BYTE:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (double)((signed char*)arr)[i + arr_index];
            break;
        case NCOMP_UBYTE:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (double)((unsigned char*)arr)[i + arr_index];
            break;
        case NCOMP_SHORT:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (double)((short*)arr)[i + arr_index];
            break;
        case NCOMP_USHORT:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (double)((unsigned short*)arr)[i + arr_index];
            break;
        case NCOMP_INT:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (double)((int*)arr)[i + arr_index];
            break;
        case NCOMP_UINT:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (double)((unsigned int*)arr)[i + arr_index];
            break;
        case NCOMP_LONG:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (double)((long*)arr)[i + arr_index];
            break;
        case NCOMP_ULONG:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (double)((unsigned long*)arr)[i + arr_index];
            break;
        case NCOMP_LONGLONG:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (double)((long long*)arr)[i + arr_index];
            break;
        case NCOMP_ULONGLONG:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (double)((unsigned long long*)arr)[i + arr_index];
            break;
        case NCOMP_LONGDOUBLE:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (double)((long double*)arr)[i + arr_index];
            break;
    }

    return;
}

void _to_float(void* arr, size_t arr_index, size_t arr_size, int arr_type, float* out_arr) {
    long i;
    switch(arr_type) {
        case NCOMP_FLOAT:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (float)((float*)arr)[i + arr_index];
            break;
        case NCOMP_DOUBLE:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (float)((double*)arr)[i + arr_index];
            break;
        case NCOMP_BOOL:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (float)((char*)arr)[i + arr_index];
            break;
        case NCOMP_BYTE:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (float)((signed char*)arr)[i + arr_index];
            break;
        case NCOMP_UBYTE:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (float)((unsigned char*)arr)[i + arr_index];
            break;
        case NCOMP_SHORT:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (float)((short*)arr)[i + arr_index];
            break;
        case NCOMP_USHORT:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (float)((unsigned short*)arr)[i + arr_index];
            break;
        case NCOMP_INT:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (float)((int*)arr)[i + arr_index];
            break;
        case NCOMP_UINT:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (float)((unsigned int*)arr)[i + arr_index];
            break;
        case NCOMP_LONG:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (float)((long*)arr)[i + arr_index];
            break;
        case NCOMP_ULONG:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (float)((unsigned long*)arr)[i + arr_index];
            break;
        case NCOMP_LONGLONG:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (float)((long long*)arr)[i + arr_index];
            break;
        case NCOMP_ULONGLONG:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (float)((unsigned long long*)arr)[i + arr_index];
            break;
        case NCOMP_LONGDOUBLE:
            for (i = 0; i < arr_size; i++)
                out_arr[i] = (float)((long double*)arr)[i + arr_index];
            break;
    }

    return;
}

ncomp_array* ncomp_array_alloc(void* array_ptr, int array_type, int array_ndim, size_t* array_shape) {
    ncomp_array* new_array = (ncomp_array*)malloc(sizeof(ncomp_array) + (array_ndim - 1) * sizeof(size_t));
    new_array->type = array_type;
    new_array->addr = array_ptr;
    new_array->ndim = array_ndim;
    int i;
    for (i = 0; i < array_ndim; i++)
        new_array->shape[i] = array_shape[i];
    return new_array;
}

ncomp_array* ncomp_array_alloc_with_data(int array_type, int array_ndim, size_t* array_shape) {
    int i;
    long n = 1;
    for (i = 0; i < array_ndim; i++)
        n *= array_shape[i];
    void* array_ptr = (void*)malloc(sizeof_ncomp_array_data(array_type) * n);
    return ncomp_array_alloc(array_ptr, array_type, array_ndim, array_shape);
}

void ncomp_array_free(ncomp_array* old_array, int keep_data_ptr) {
    /* free ptr by default */
    if (!keep_data_ptr)
        free(old_array->addr);

    free(old_array);

    return;
}

int sizeof_ncomp_array_data(int array_type) {
    int array_type_size = 0;
    switch(array_type) {
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
void coerce_missing(
int     type_x,
int     has_missing_x,
void   *missing_x,
double *missing_dx,
float  *missing_rx)
{
/*
 * Check for missing value and coerce if neccesary.
 */
  if(has_missing_x) {
/*
 * Coerce missing value to double.
 */
    _ncomp_coerce_to_double(missing_x, type_x, NULL, missing_dx, NULL, 1);

    if(type_x != NCOMP_DOUBLE && missing_rx != NULL) {
        _to_float(missing_x, 0, 1, type_x, missing_rx);
        _ncomp_coerce_to_float(missing_x, type_x, NULL, missing_rx, NULL, 1);
    }
  }
  else {
    if(missing_dx != NULL) {
/*
 * Get the default missing value, just in case.
 */
      if(type_x != NCOMP_DOUBLE) {
        *missing_dx = NC_FILL_DOUBLE;
        if(missing_rx != NULL) {
          *missing_dx = NC_FILL_FLOAT;
        }
      }
      else {
        *missing_dx = NC_FILL_DOUBLE;
      }
    }
  }
}

/*
 * Coerce data to double, or just return a pointer to it if
 * it is already double.
 */
double *coerce_input_double(
void              *x,
int type_x,
size_t         size_x,
int               has_missing_x,
void         *missing_x,
double         *missing_dx)
{
  double *dx;
/*
 * Coerce x to double if necessary.
 */
  if(type_x != NCOMP_DOUBLE) {
    dx = (double*)calloc(size_x,sizeof(double));
    if( dx == NULL ) return(NULL);
    if(has_missing_x) {
      /*
      _Nclcoerce((NclTypeClass)nclTypedoubleClass,
                 (void*)dx,
                 x,
                 size_x,
                 missing_x,
                 missing_dx,
                 _NclTypeEnumToTypeClass(_NclBasicDataTypeToObjType(type_x)));
      */
        _ncomp_coerce(x, type_x, missing_x, dx, NCOMP_DOUBLE, (void*) missing_dx, size_x);
    }
    else {
      /*
      _Nclcoerce((NclTypeClass)nclTypedoubleClass,
                 (void*)dx,
                 x,
                 size_x,
                 NULL,
                 NULL,
                 _NclTypeEnumToTypeClass(_NclBasicDataTypeToObjType(type_x)));
      */
        _ncomp_coerce(x, type_x, NULL, dx, NCOMP_DOUBLE, NULL, size_x);
    }
  }
  else {
/*
 * x is already double.
 */
    dx = (double*)x;
  }
  return(dx);
}

void _ncomp_coerce(
    void   *from_ptr,
    int     from_type,
    void   *from_missing,
    void   *to_ptr,
    int     to_type,
    void   *to_missing,
    size_t  num
) {
    switch(to_type) {
        case NCOMP_DOUBLE:
            _ncomp_coerce_to_double(from_ptr, from_type, from_missing, (double*) to_ptr, (double*) to_missing, num);
            break;
        case NCOMP_FLOAT:
            _ncomp_coerce_to_float(from_ptr, from_type, from_missing, (float*) to_ptr, (float*) to_missing, num);
            break;
        case NCOMP_BOOL:
            _ncomp_coerce_to_bool(from_ptr, from_type, from_missing, (char*) to_ptr, (char*) to_missing, num);
            break;
        case NCOMP_BYTE:
            _ncomp_coerce_to_byte(from_ptr, from_type, from_missing, (signed char*) to_ptr, (signed char*) to_missing, num);
            break;
        case NCOMP_UBYTE:
            _ncomp_coerce_to_ubyte(from_ptr, from_type, from_missing, (unsigned char*) to_ptr, (unsigned char*) to_missing, num);
            break;
        case NCOMP_SHORT:
            _ncomp_coerce_to_short(from_ptr, from_type, from_missing, (short*) to_ptr, (short*) to_missing, num);
            break;
        case NCOMP_USHORT:
            _ncomp_coerce_to_ushort(from_ptr, from_type, from_missing, (unsigned short*) to_ptr, (unsigned short*) to_missing, num);
            break;
        case NCOMP_INT:
            _ncomp_coerce_to_int(from_ptr, from_type, from_missing, (int*) to_ptr, (int*) to_missing, num);
            break;
        case NCOMP_UINT:
            _ncomp_coerce_to_uint(from_ptr, from_type, from_missing, (unsigned int*) to_ptr, (unsigned int*) to_missing, num);
            break;
        case NCOMP_LONG:
            _ncomp_coerce_to_long(from_ptr, from_type, from_missing, (long*) to_ptr, (long*) to_missing, num);
            break;
        case NCOMP_ULONG:
            _ncomp_coerce_to_ulong(from_ptr, from_type, from_missing, (unsigned long*) to_ptr, (unsigned long*) to_missing, num);
            break;
        case NCOMP_LONGLONG:
            _ncomp_coerce_to_longlong(from_ptr, from_type, from_missing, (long long*) to_ptr, (long long*) to_missing, num);
            break;
        case NCOMP_ULONGLONG:
            _ncomp_coerce_to_ulonglong(from_ptr, from_type, from_missing, (unsigned long long*) to_ptr, (unsigned long long*) to_missing, num);
            break;
        case NCOMP_LONGDOUBLE:
            _ncomp_coerce_to_longdouble(from_ptr, from_type, from_missing, (long double*) to_ptr, (long double*) to_missing, num);
            break;
    }
}

/*
 * Coerce a non-contiguous subset of the data to double.
 */
void coerce_subset_input_double_step(
void   *x,
double *tmp_x,
size_t  index_x,
size_t  step_x,
int     type_x,
size_t  size_x,
int     has_missing_x,
void   *missing_x,
double *missing_dx
)
{ 
  size_t i, ii;
  int sizeof_x;
  
/*
 * typeclass_x is what allows us to get the size of the type of x.
 */
  sizeof_x = sizeof_ncomp_array_data(type_x);
/*
 * Coerce x to double.
 */
  if(has_missing_x) {
/*
 * Coerce subset to double, with missing values.
 */ 
    for(i = 0; i < size_x; i++ ) {
      ii = step_x*i;
/*
      _Nclcoerce((NclTypeClass)nclTypedoubleClass,
                 &tmp_x[i],
                 (void*)((char*)x+(index_x+ii)*(sizeof_x)),
                 1,
                 missing_x,
                 missing_dx,
                 typeclass_x);
*/
      _ncomp_coerce_to_double((void*)((char*)x+(index_x+ii)*(sizeof_x)), type_x, missing_x, &tmp_x[i], missing_dx, 1);
    }
  }
  else {
/*
 * Coerce subset to double, with no missing values.
 */ 
    for(i = 0; i < size_x; i++ ) {
      ii = step_x*i;
      _ncomp_coerce_to_double((void*)((char*)x+(index_x+ii)*(sizeof_x)), type_x, NULL, &tmp_x[i], NULL, 1);
    }
  }
}

/*
 * Coerce a non-contiguous subset of the data to float.
 */
void coerce_subset_input_float_step(
void   *x,
float  *tmp_x,
size_t  index_x,
size_t  step_x,
int     type_x,
size_t  size_x,
int     has_missing_x,
void   *missing_x,
float  *missing_fx
)
{ 
  size_t i, ii;
  int sizeof_x;
  
/*
 * typeclass_x is what allows us to get the size of the type of x.
 */
  sizeof_x = sizeof_ncomp_array_data(type_x);
/*
 * Coerce x to float.
 */
  if(has_missing_x) {
/*
 * Coerce subset to float, with missing values.
 */ 
    for(i = 0; i < size_x; i++ ) {
      ii = step_x*i;
      _ncomp_coerce_to_float((void*)((char*)x+(index_x+ii)*(sizeof_x)), type_x, missing_x, &tmp_x[i], missing_fx, 1);
    }
  }
  else {
/*
 * Coerce subset to float, with no missing values.
 */ 
    for(i = 0; i < size_x; i++ ) {
      ii = step_x*i;
      _ncomp_coerce_to_float((void*)((char*)x+(index_x+ii)*(sizeof_x)), type_x, NULL, &tmp_x[i], NULL, 1);
    }
  }
}

/*
 * Sets a subset of non-contiguous output data to missing.
 */
void set_subset_output_missing_step(
void   *x,
size_t  index_x,
size_t  step_x,
int     type_x,
size_t  size_x,
double  missing_x)
{
  size_t i;
  for(i = 0; i < size_x; i++) {
    if(type_x != NCOMP_DOUBLE) {
      ((float*)x)[index_x + i*step_x] = (float)missing_x;
    }
    else {
      ((double*)x)[index_x + i*step_x] = missing_x;
    }
  }
}

/*
 * Copy double data back to double or float non-contiguous array,
 * using a void array. 
 */
void coerce_output_float_or_double_step(
void   *x,
double *dx,
int     type_x,
size_t  size_x,
size_t  index_x,
size_t  step_x
)
{
  size_t i;

  if(type_x == NCOMP_DOUBLE) {
    for( i = 0; i < size_x; i++ ) ((double*)x)[index_x+(step_x*i)] = dx[i];
  }
  else {
    for( i = 0; i < size_x; i++ ) ((float*)x)[index_x+(step_x*i)]  = (float)dx[i];
  }
}
