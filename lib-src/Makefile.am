lib_LTLIBRARIES = libcdd.la

libcdd_la_SOURCES = \
cddcore.c \
cddlp.c	\
cddmp.c \
cddio.c \
cddlib.c \
cddproj.c \
setoper.c
libcdd_la_CPPFLAGS = -UGMPRATIONAL

AM_LDFLAGS = -version-info $(libcdd_version_info) $(CDD_LDFLAGS)

pkginclude_HEADERS = \
cdd.h \
cddmp.h \
cddtypes.h \
setoper.h \
splitmix64.h

include ./Makefile.gmp.am
