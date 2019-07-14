SUBDIRS = fortran
lib_LTLIBRARIES = libncomp.la
libncomp_la_SOURCES = util.c ncomp.h linint2.c
libncomp_la_CFLAGS = $(AM_CFLAGS)
libncomp_la_LIBADD = fortran/libncomp_fort.la
include_HEADERS = ncomp.h
ncompincludedir = $(includedir)/ncomp
ncompinclude_HEADERS = util.h constants.h wrapper.h
