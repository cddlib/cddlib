bin_PROGRAMS = \
scdd \
lcdd \
redcheck \
fourier \
projection \
adjacency \
allfaces \
testcdd1 \
testshoot \
testcdd2 \
testlp1 \
testlp2 \
testlp3 \
cddexec

# Build these binaries but do not install them; these are sample applications
# that are mostly interesting to developers.
noinst_PROGRAMS = \
redexter \
redundancies \
redundancies_clarkson

LDADD = ../lib-src/libcdd.la
AM_CPPFLAGS = -I$(top_srcdir)/lib-src

include ./Makefile.gmp.am
