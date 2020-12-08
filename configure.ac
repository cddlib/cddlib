dnl Process this file with autoconf to produce a configure script, e.g., by running ./bootstrap.

dnl We need at least autoconf 2.69 for this configure.ac to work.
AC_PREREQ([2.69])

dnl When changing the version, also change it in cddtypes.h's dd_DDVERSION
AC_INIT([cddlib], [0.94m])
dnl Version number of the library c:r:a.
dnl Set to c+1:0:0 on a non-backwards-compatible change.
dnl Set to c+1:0:a+1 on a backwards-compatible change.
dnl Set to c:r+1:a if the interface is unchanged.
AC_SUBST([libcdd_version_info], [1:2:1])

AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([1.14.1 foreign tar-ustar])

dnl Find C compiler.
AC_PROG_CC
dnl Set INSTALL or use install-sh.
AC_PROG_INSTALL
dnl Our "make distcheck" needs egrep.
AC_PROG_EGREP
dnl Configure libtool.
LT_INIT([win32-dll])

AC_CANONICAL_HOST
dnl libtool requires "-no-undefined" for win32 dll
AC_SUBST(CDD_LDFLAGS)
case $host_os in
    *cygwin*|*mingw*)
    if test x"$enable_shared" = "xyes"; then
        CDD_LDFLAGS="$CDD_LDFLAGS -no-undefined"
    fi
;;
esac

dnl Look for gmp. But do not prepend -lgmp to LIBS as we do not want to link everything against gmp.
dnl Check for GMP headers and library
AC_CHECK_HEADER(gmp.h, [
  AC_CHECK_LIB([gmp], [__gmpz_init], [true])
])

dnl We do not build the gmp enabled library if we don't have gmp.
AM_CONDITIONAL([GMP], [test "x$ac_cv_lib_gmp___gmpz_init" = "xyes"])
AS_IF([test "x$ac_cv_lib_gmp___gmpz_init" = "xyes"], [CDD_LIBS="-lcdd -lcddgmp"], [CDD_LIBS="-lcdd"])
AC_SUBST(CDD_LIBS)

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

AC_CONFIG_FILES([doc/Makefile lib-src/Makefile src/Makefile Makefile cddlib.pc])
AC_OUTPUT
