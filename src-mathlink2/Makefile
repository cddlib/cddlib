# Makefile for cddmathlink2
#
# see README
#
# make ml     to produce MathLink program cddmathlink2 with C double arithmetic
# make clean      to clean all objects and cddmathlink2 

# C-compiler
CC=/usr/local/bin/mcc
#CC=gcc
#CC=cc

# Archive maintainer ar
AR=ar

# Libtool ranlib
RL=ranlib

# ----------- cddlib Compilation Environment ---------- 
# Compiler options
CFLAGS = -g
#CFLAGS = -pg -Wall -O
#CFLAGS = -O3

# --------------- MathLink Environment --------------- 
#  Ask a system administrator if you are not familiar with 
#  the local Mathematica installation
#
# MathLink Top Directory  
MLINKDIR = /usr/local

# Subdirectory containing mprep 
MLBINDIR = ${MLINKDIR}/bin

# Subdirectory containing the ML header file mathlink.h 
MLINCDIR = ${MLINKDIR}/include

# Subdirectory containing the Library libML.a 
MLLIBDIR = ${MLINKDIR}/lib

# MathLink Compilar Flags (check mcc scripts for appropriate flags)
#MLFLAGS = -lML -lm -lsocket -lnsl   #Solaris
MLFLAGS = -lML # Generic


########## You shouldn't have to change anything after this point ##########
# Location of include files for cddlib compilation
IFLAGS = -I../lib-src -I../lib-src-gmp

# C-Library location  for cddlib compilation
LFLAGS = -L../lib-src -L../lib-src-gmp

# Library and Include FLAGS together
LIFLAGS = $(LFLAGS) $(IFLAGS) $(CFLAGS)

ml : cddmathlink2 cddmathlink2gmp

cddmathlink2 : cddmathlinktm.o cddmathlink.o cddmlio.o
	${CC} $(LIFLAGS) -I${MLINCDIR} cddmlio.o cddmathlinktm.o cddmathlink.o -L${MLLIBDIR} ${MLFLAGS} ../lib-src/libcdd.a -o cddmathlink2

cddmathlink2gmp : cddmathlinktm.o cddmathlink_gmp.o cddmlio_gmp.o
	${CC} $(LIFLAGS) -I${MLINCDIR} cddmlio_gmp.o cddmathlinktm.o cddmathlink_gmp.o -L${MLLIBDIR} ${MLFLAGS} -lgmp ../lib-src-gmp/libcddgmp.a -o cddmathlink2gmp

cddmathlinktm.o: cddmathlinktm.c $(HFILES)
	$(CC) $(LIFLAGS) -c cddmathlinktm.c

cddmathlink.o: cddmathlink.c $(HFILES)
	$(CC) $(LIFLAGS) -c cddmathlink.c

cddmathlink_gmp.o: cddmathlink.c $(HFILES)
	$(CC) $(LIFLAGS) -DGMPRATIONAL -c cddmathlink.c -o cddmathlink_gmp.o

cddmlio.o: cddmlio.c $(HFILES)
	$(CC) $(LIFLAGS) -c cddmlio.c

cddmlio_gmp.o: cddmlio.c $(HFILES)
	$(CC) $(LIFLAGS) -DGMPRATIONAL -c cddmlio.c -o cddmlio_gmp.o

cddmathlinktm.c : cddmathlink.tm
	${MLBINDIR}/mprep $? > $@

clean:
	rm -rf core *.out  *.o cddmathlinktm.* cddmathlink2 cddmathlink2gmp

all: ml
