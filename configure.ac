dnl Process this file with autoconf to produce a configure script, e.g., by running ./bootstrap.
AC_INIT([cddlib], [0.94i])
dnl Version number of the library c:r:a.
dnl Set to c+1:0:0 on a non-backwards-compatible change.
dnl Set to c+1:0:a+1 on a backwards-compatible change.
dnl Set to c:r+1:a if the interface is unchanged.
AC_SUBST([libcdd_version_info], [0:0:0])

AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE

dnl Find C compiler.
AC_PROG_CC
dnl Set INSTALL or use install-sh.
AC_PROG_INSTALL
dnl Configure libtool.
LT_INIT

dnl Look for gmp. But do not prepend -lgmp to LIBS as we do not want to link everything against gmp.
AC_CHECK_LIB([gmp], [__gmpz_init], [true])
dnl We do not build the gmp enabled library if we don't have gmp.
AM_CONDITIONAL([GMP], [test "x$ac_cv_lib_gmp___gmpz_init" = "xyes" ])

AC_CONFIG_FILES([lib-src/Makefile src/Makefile lib-src-gmp/Makefile src-gmp/Makefile Makefile])
AC_OUTPUT
