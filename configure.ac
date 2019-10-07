AC_INIT([NComp], [0.1], [hallock@ucar.edu], [ncomp])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CC
AC_PROG_F77
AC_PROG_FC
AC_PROG_CXX
AM_PROG_AR
AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
 src/Makefile
 src/fortran/Makefile
 src/wrappers/Makefile
 test/Makefile
])
LT_INIT
AC_ARG_WITH([blas],
            [--with-blas=<path to blas installation prefix>],
            [case "${withval}" in
                 no) blas_dir=${prefix}/lib ;;
                 *)  blas_dir="${withval}" ;;
             esac],
            [blas_dir="${prefix}/lib"])

AC_ARG_WITH([lapack],
            [--with-lapack=<path to lapack installation prefix>],
            [case "${withval}" in
                 no) lapack_dir=${prefix}/lib ;;
                 *)  lapack_dir="${withval}" ;;
             esac],
            [lapack_dir="${prefix}/lib"])

AC_SUBST(BLAS_DIR,[$blas_dir])
AC_SUBST(LAPACK_DIR,[$lapack_dir])
AC_OUTPUT
