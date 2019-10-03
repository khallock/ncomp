check_PROGRAMS = 	test_linint2 \
								 	test_eofunc_01 test_eofunc_02 test_eofunc_03 test_eofunc_04 test_eofunc_05 test_eofunc_06 \
									test_eofunc_n_01 test_eofunc_n_02 test_eofunc_n_03 test_eofunc_n_04 test_eofunc_n_05
AM_CPPFLAGS = -I$(srcdir)/../src
AM_LDFLAGS = -lm -L$(srcdir)/../src/.libs -lncomp -Wl,-rpath,$(prefix)/lib
TESTS = $(check_PROGRAMS)
