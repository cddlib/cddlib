lib_LIBRARIES = libcddgmp.a

libcddgmp_a_SOURCES = \
cddcore.c \
cddlp.c	\
cddmp.c \
cddio.c \
cddlib.c \
setoper.c

include_HEADERS = \
cdd.h \
cddmp.h \
cddtypes.h \
setoper.h

CC       = $(CXX)
CPPFLAGS = -DGMPRATIONAL
CXXFLAGS = -O3
INCLUDES = -I$(gmpincludedir)

