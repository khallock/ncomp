SUBDIRS = fortran wrappers
lib_LTLIBRARIES = libncomp.la
libncomp_la_SOURCES = util.cpp
libncomp_la_CFLAGS = $(AM_CFLAGS)
libncomp_la_LIBADD = fortran/libncomp_fort.la wrappers/libncomp_wrappers.la -lgfortran -llapack -lblas
include_HEADERS = ncomp.h
nobase_include_HEADERS = ncomp/constants.h ncomp/types.h ncomp/util.h ncomp/wrapper.h
EXTRA_DIST = ncomp_internal/
