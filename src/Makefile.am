bin_PROGRAMS = \
scdd \
lcdd \
redcheck \
fourier \
projection \
adjacency \
testcdd1 \
testshoot \
testcdd2 \
testlp1 \
testlp2 \
testlp3
#cddmathlink

scdd_SOURCES                   = simplecdd.c
lcdd_SOURCES                   = lcdd.c
redcheck_SOURCES               = redcheck.c
fourier_SOURCES                = fourier.c
projection_SOURCES             = projection.c
adjacency_SOURCES              = adjacency.c
testcdd1_SOURCES               = testcdd1.c
testshoot_SOURCES              = testshoot.c
testcdd2_SOURCES               = testcdd2.c
testlp1_SOURCES                = testlp1.c
testlp2_SOURCES                = testlp2.c
testlp3_SOURCES                = testlp3.c
# cddmathlink_SOURCES                = cddmathlink.c cddmlio.h cddmlio.c

CC           = $(CXX)
LDADD        = -L../lib-src -lcdd
INCLUDES     = -I../lib-src
CPPFLAGS     = -UGMPRATIONAL
CXXFLAGS     = -O2
