#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.69])
AC_INIT(stap, 2.0)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_SRCDIR([src/seq2binding.cpp])
AM_INIT_AUTOMAKE

AC_ARG_ENABLE([debug],
[  --enable-debug	Turn on debugging],
[case "${enableval}" in
yes) debug=true ;;
no) debug=false ;;
*) AC_MSG_ERROR([bad value ${enablevar} for --enable-debug]) ;;
esac], [debug=false])
AM_CONDITIONAL([DEBUG], [test x${debug} = xtrue])

AC_ARG_WITH([StrongerSiteAttr],
[  --with-StrongerSiteAttr,	compile-time option present in original Makefile, but not in code.],
[ stronger_site_attr=$withval ],
[ stronger_site_attr=0])

AC_ARG_WITH([StrongerSiteAttr],
[  --with-RandomGen,     compile-time option present in original Makefile, but not in code.],
[ random_gen=$withval ],
[ random_gen=0])

# Checks for programs.
AC_PROG_CXX
AC_LANG(C++)

# Checks for libraries.
AX_PATH_GSL([1.8], [], [AC_MSG_ERROR(Could not find the required version of GSL)])

# Order matters here.
AC_CHECK_LIB([m], [cos])
AC_CHECK_LIB([gslcblas], [cblas_dgemm])
AC_CHECK_LIB([gsl], [gsl_blas_dgemm])

# Checks for header files.
AC_CHECK_HEADERS([gsl/gsl_linalg.h gsl/gsl_math.h gsl/gsl_multimin.h gsl/gsl_randist.h])


# Checks for typedefs, structures, and compiler characteristics.
#AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([floor pow sqrt strchr])

AC_CONFIG_FILES([Makefile
		 src/Makefile])

#TODO : add a set of regression tests
# Setup the test suite.
#AC_CONFIG_TESTDIR(tests)
#AC_CONFIG_FILES([tests/Makefile tests/atlocal])
#AC_CONFIG_FILES([tests/testexec], [chmod +x tests/testexec])
#AM_MISSING_PROG([AUTOM4TE], [autom4te])

AC_OUTPUT
