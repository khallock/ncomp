#ifndef NCOMP_CONSTANTS_H
#define NCOMP_CONSTANTS_H

#define NC_FILL_BYTE    ((signed char)-127)
#define NC_FILL_CHAR    ((char)0)
#define NC_FILL_SHORT   ((short)-32767)
#define NC_FILL_INT     (-2147483647L)
#define NC_FILL_FLOAT   (9.9692099683868690e+36f) /* near 15 * 2^119 */
#define NC_FILL_DOUBLE  (9.9692099683868690e+36)
#define NC_FILL_UBYTE   (255)
#define NC_FILL_USHORT  (65535)
#define NC_FILL_UINT    (4294967295U)
#define NC_FILL_INT64   ((long long)-9223372036854775806LL)
#define NC_FILL_UINT64  ((unsigned long long)18446744073709551614ULL)
#define NC_FILL_STRING  ((char *)"")

#endif
