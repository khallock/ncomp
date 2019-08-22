#ifndef NCOMP_UTIL_H
#define NCOMP_UTIL_H

#include <ncomp/types.h>

ncomp_array* ncomp_array_alloc(void*, int, int, size_t*);
size_t sizeof_ncomp_array_data(int);

template <typename T>
void convert_to(void *in_arr, size_t arr_size, size_t arr_offset, int arr_type,
                T *out_arr);

void coerce_missing(int, int, const ncomp_missing *, double *, float *);

template <typename T,
          typename = typename std::enable_if<
              std::is_same<T, double>::value || std::is_same<T, float>::value ||
                  std::is_same<T, int>::value ||
                  std::is_same<T, unsigned int>::value ||
                  std::is_same<T, long>::value,
              T>::type>
T *coerce_input_T(void *, int, size_t, int, void *, T *);

template <typename T>
void _ncomp_coerce_internal(void *from_ptr, int from_type, void *from_missing,
                            T *to_ptr, T *to_missing, bool has_missing,
                            size_t num);

void _ncomp_coerce(void *from_ptr, int from_type, void *from_missing,
                   void *to_ptr, int to_type, void *to_missing, size_t num);

#endif
