#ifndef NCOMP_TYPES_H
#define NCOMP_TYPES_H

#include <stddef.h>

/* enumerated types, adapted from numpy's ndarraytypes.h */
enum NcompTypes {
  NCOMP_BOOL = 0,
  NCOMP_BYTE = 1,
  NCOMP_UBYTE = 2,
  NCOMP_SHORT = 3,
  NCOMP_USHORT = 4,
  NCOMP_INT = 5,
  NCOMP_UINT = 6,
  NCOMP_LONG = 7,
  NCOMP_ULONG = 8,
  NCOMP_LONGLONG = 9,
  NCOMP_ULONGLONG = 10,
  NCOMP_FLOAT = 11,
  NCOMP_DOUBLE = 12,
  NCOMP_LONGDOUBLE = 13
};

typedef union {
    char                msg_bool;
    signed char         msg_byte;
    unsigned char       msg_ubyte;
    short               msg_short;
    unsigned short      msg_ushort;
    int                 msg_int;
    unsigned int        msg_uint;
    long                msg_long;
    unsigned long       msg_ulong;
    long long           msg_longlong;
    unsigned long long  msg_ulonglong;
    float               msg_float;
    double              msg_double;
    long double         msg_longdouble;
} ncomp_missing;

typedef struct {
    /* ordered for efficient memory packing */
    int type;
    int ndim;
    void *addr;
    int has_missing; // = 0;
    ncomp_missing msg;
    size_t shape[1];
} ncomp_array;

typedef struct {
  char* name; // Name of the attribute
  ncomp_array* value; // the value for the attribute
} single_attribute;

typedef struct {
  int nAttribute; // number of attributes
  single_attribute * attribute_array; // the attributes.
} attributes;

#endif // NCOMP_TYPES_H
