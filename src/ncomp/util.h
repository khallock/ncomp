#ifndef NCOMP_UTIL_H
#define NCOMP_UTIL_H

#include <ncomp/types.h>
#include <limits.h>

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
extern "C" {
#endif

ncomp_array* ncomp_array_alloc(void*, int, int, size_t*);
ncomp_array* ncomp_array_alloc_scalar(void* data_ptr, int data_type);
void ncomp_array_free(ncomp_array*, int);
size_t sizeof_ncomp_array_data(int);

ncomp_single_attribute* create_ncomp_single_attribute_from_ncomp_array(
  char * name,
  ncomp_array* value);

ncomp_single_attribute * create_ncomp_single_attribute_char(
    char * name,
    char * data);

ncomp_single_attribute* create_ncomp_single_attribute(
  char * name,
  void * data_ptr,
  int data_type,
  int data_ndim,
  size_t * data_shape);

ncomp_single_attribute* create_ncomp_single_attribute_from_1DArray(
    char * name,
    void * data_ptr,
    int data_type,
    size_t nelem);


ncomp_single_attribute* create_ncomp_single_attribute_from_scalar(
  char * name,
  void * data_ptr,
  int data_type);

ncomp_attributes* ncomp_attributes_allocate(int nAttribute);

int hasAttribute(
  const ncomp_attributes * attributeList,
  const char* attributeName,
  int* attributePosInList);

ncomp_single_attribute*  getAttributeOrDefault_ncomp_single_attribute(
  const ncomp_attributes * attributeList,
  const char* attributeName,
  const ncomp_single_attribute* defaultValue);

ncomp_single_attribute*  getAttribute(
  const ncomp_attributes * attributeList,
  const char* attributeName);

void* getAttributeOrDefault(
  const ncomp_attributes * attributeList,
  const char* attributeName,
  const void * defaultValue);

void print_ncomp_array(const char * name, const ncomp_array * in);
void print_ncomp_attributes(const ncomp_attributes * in);

#ifdef __cplusplus /* If this is a C++ compiler, end C linkage */
}
#endif

#endif
