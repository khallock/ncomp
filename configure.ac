AC_INIT([NComp], [0.1], [hallock@ucar.edu], [ncomp])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AX_CHECK_ENABLE_DEBUG([no], [DEBUG])
AC_PROG_CC
AC_PROG_FC
AC_PROG_CXX
AM_PROG_AR
AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_HEADERS([ncomp_config.h])
AC_CONFIG_FILES([
 Makefile
 src/Makefile
 src/fortran/Makefile
 src/wrappers/Makefile
 test/Makefile
])
LT_INIT
AX_BLAS
AX_LAPACK
AC_OUTPUT
