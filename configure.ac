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
dnl Our "make distcheck" needs egrep.
AC_PROG_EGREP
dnl Configure libtool.
LT_INIT

dnl Look for gmp. But do not prepend -lgmp to LIBS as we do not want to link everything against gmp.
AC_CHECK_LIB([gmp], [__gmpz_init], [true])
dnl We do not build the gmp enabled library if we don't have gmp.
AM_CONDITIONAL([GMP], [test "x$ac_cv_lib_gmp___gmpz_init" = "xyes"])

dnl Check for latex to build the documentation dvi
AC_CHECK_PROGS([latex], [latex])
AM_CONDITIONAL([LATEX], [test "x$ac_cv_prog_latex" != "x"])

dnl Check for pdflatex to build the documentation pdf
AC_CHECK_PROGS([pdflatex], [pdflatex])
AM_CONDITIONAL([PDF], [test "x$ac_cv_prog_pdflatex" != "x"])

dnl Check for latex2html to build the documentation in HTML format
AC_CHECK_PROGS([latex2html], [latex2html])
AM_CONDITIONAL([HTML], [test "x$ac_cv_prog_latex2html" != "x"])

dnl Check for dvips to build the documentation ps
AC_CHECK_PROGS([dvips], [dvips])
AM_CONDITIONAL([PS], [test "x$ac_cv_prog_dvips" != "x"])

AC_CONFIG_FILES([doc/Makefile lib-src/Makefile src/Makefile lib-src-gmp/Makefile src-gmp/Makefile Makefile])
AC_OUTPUT
