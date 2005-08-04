# Makefile for cddmathlink
#
# see README
#
# make ml     to produce MathLink program cddmathlink with C double arithmetic
# make clean      to clean all objects and cddmathlink 

# C-compiler
CC=/usr/local/bin/mcc
#CC=cc

# Archive maintainer ar
AR=ar

# Libtool ranlib
RL=ranlib

# ----------- cddlib Compilation Environment ---------- 
# Compiler options
CFLAGS = -g -O3
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
MLFLAGS = -lML -lm  # Generic


########## You shouldn't have to change anything after this point ##########
# Location of include files for cddlib compilation
IFLAGS =  -I../lib-src

# C-Library location  for cddlib compilation
LFLAGS =  -L../lib-src

# Library and Include FLAGS together
LIFLAGS = $(LFLAGS) $(IFLAGS) $(CFLAGS)

ml : cddmathlink

cddmathlink : cddmathlinktm.o cddmathlink.o cddmlio.o
	${CC} -O3 -I${MLINCDIR} cddmlio.o cddmathlinktm.o cddmathlink.o -L${MLLIBDIR} ${MLFLAGS} $(LFLAGS) ../lib-src/libcdd.a -o cddmathlink

cddmathlinktm.o: cddmathlinktm.c $(HFILES)
	$(CC) $(LIFLAGS) -c cddmathlinktm.c

cddmathlink.o: cddmathlink.c $(HFILES)
	$(CC) $(LIFLAGS) -c cddmathlink.c

cddmlio.o: cddmlio.c $(HFILES)
	$(CC) $(LIFLAGS) -c cddmlio.c

cddmathlinktm.c : cddmathlink.tm
	${MLBINDIR}/mprep $? > $@

clean:
	rm -rf core *.out  *.o cddmathlinktm.* cddmathlink

all: ml
