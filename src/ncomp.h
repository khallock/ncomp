#include <stddef.h>

#ifndef NCOMP_H
#define NCOMP_H

/* enumerated types, adapted from numpy's ndarraytypes.h */
enum NCOMP_TYPES {NCOMP_BOOL=0,
#define           NCOMP_BOOL 0
                  NCOMP_BYTE,
#define           NCOMP_BYTE 1
                  NCOMP_UBYTE,
#define           NCOMP_UBYTE 2
                  NCOMP_SHORT,
#define           NCOMP_SHORT 3
                  NCOMP_USHORT,
#define           NCOMP_USHORT 4
                  NCOMP_INT,
#define           NCOMP_INT 5
                  NCOMP_UINT,
#define           NCOMP_UINT 6
                  NCOMP_LONG,
#define           NCOMP_LONG 7
                  NCOMP_ULONG,
#define           NCOMP_ULONG 8
                  NCOMP_LONGLONG,
#define           NCOMP_LONGLONG 9
                  NCOMP_ULONGLONG,
#define           NCOMP_ULONGLONG 10
                  NCOMP_FLOAT,
#define           NCOMP_FLOAT 11
                  NCOMP_DOUBLE,
#define           NCOMP_DOUBLE 12
                  NCOMP_LONGDOUBLE
#define           NCOMP_LONGDOUBLE 13
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
    int     type;
    int     ndim;
    void*   addr;
    int     has_missing;// = 0;
    ncomp_missing msg;
    size_t  shape[1]; /* use ndim to malloc extra space for dims */
} ncomp_array;

#include "constants.h"
#include "wrapper.h"

#endif
