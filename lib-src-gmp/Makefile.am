lib_LIBRARIES = libcddgmp.a

libcddgmp_a_SOURCES = \
cddcore.c \
cddlp.c	\
cddmp.c \
cddio.c \
cddlib.c \
cddproj.c \
cddcore_f.c \
cddlp_f.c	\
cddmp_f.c \
cddio_f.c \
cddlib_f.c \
cddproj_f.c \
setoper.c

include_HEADERS = \
cdd.h \
cddmp.h \
cddtypes.h \
cdd_f.h \
cddmp_f.h \
cddtypes_f.h \
setoper.h

CC       = $(CXX)
CPPFLAGS = -DGMPRATIONAL
CXXFLAGS = -O3
INCLUDES = -I$(gmpincludedir)

