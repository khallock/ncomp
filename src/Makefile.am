SUBDIRS = fortran
lib_LTLIBRARIES = libncomp.la
libncomp_la_SOURCES = util.c ncomp.h
libncomp_la_CFLAGS = $(AM_CFLAGS)
if LINUX
libncomp_la_LDFLAGS = -Wl,--whole-archive,fortran/libncomp_fort.a -Wl,--no-whole-archive
endif
if OSX
libncomp_la_LDFLAGS = -Wl,-all_load,fortran/libncomp_fort.a
endif
