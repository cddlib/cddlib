# Makefile for cddlib, cddmathlink and examples 
#
# see README
#
# make all    to produce cdd/dplex objects and test programs (cdd_test*, dplex_test* )
# make ml     to produce MathLink program cddmathlink
# make everything   to both of the above 
# make clean  to clean all objects, library and executables 

# C-compiler
CC=gcc
#CC=cc

# --------------- MathLink Environment --------------- 
#  Ignore this part if you are not using cddmathlink;  
#  go down to "cddlib Compilation Environment". 
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
MLFLAGS = -lML -lm                  #NEXTSTEP


ml : cddmathlink

cddmathlink : cddarith.o cddio.o cddlib.o cddmlio.o setoper.o cddmathlinktm.o cddmathlink.o
	${CC} -O3 -I${MLINCDIR} cddarith.o cddio.o cddlib.o cddmlio.o setoper.o cddmathlinktm.o cddmathlink.o -L${MLLIBDIR} ${MLFLAGS} -o cddmathlink

.c.o :
	${CC} -c -O3 -I${MLINCDIR} $<

cddmathlinktm.c : cddmathlink.tm
	${MLBINDIR}/mprep $? > $@

# ----------- cddlib Compilation Environment ---------- 
# Location of include files for cddlib compilation
INCLUDEDIR = .

# Compiler options
#CFLAGS = -g -O -I$(INCLUDEDIR)
#CFLAGS = -pg -Wall -O -I$(INCLUDEDIR)
CFLAGS = -O3 -I$(INCLUDEDIR)

# Standard C-Library Location
# LFLAGS = -L.
LFLAGS =

########## You shouldn't have to change anything after this point ##########

COFILES = cddio.o cddlib.o cddarith.o setoper.o
DOFILES = dplex.o setoper.o

cddio.o: cddio.c dplex.h cddlib.h
	$(CC) $(CFLAGS) -c cddio.c

cddmlio.o: cddmlio.c dplex.h cddlib.h
	$(CC) $(CFLAGS) -c cddmlio.c

cddarith.o: cddarith.c dplex.h cddlib.h
	$(CC) $(CFLAGS) -c cddarith.c

cddlib.o: cddlib.c dplex.h cddlib.h
	$(CC) $(CFLAGS) -c cddlib.c

dplex.o: dplex.c dplex.h setoper.h
	$(CC) $(CFLAGS) -c dplex.c

dplex_test1.o: dplex_test1.c dplex.h
	$(CC) $(CFLAGS) -c dplex_test1.c

dplex_test2.o: dplex_test2.c dplex.h
	$(CC) $(CFLAGS) -c dplex_test2.c

dplex_test3.o: dplex_test3.c dplex.h
	$(CC) $(CFLAGS) -c dplex_test3.c

cdd_test1.o: cdd_test1.c dplex.h cddlib.h
	$(CC) $(CFLAGS) -c cdd_test1.c

cdd_test2.o: cdd_test2.c dplex.h cddlib.h
	$(CC) $(CFLAGS) -c cdd_test2.c

setoper.o: setoper.c
	$(CC) $(CFLAGS) -c setoper.c

cdd_test1: cdd_test1.o $(COFILES)
	$(CC) $(CFLAGS) $(COFILES) cdd_test1.o -o cdd_test1
	./cdd_test1 < ../ine/reg24-5.ine

cdd_test2: cdd_test2.o $(COFILES)
	$(CC) $(CFLAGS) $(COFILES) cdd_test2.o -o cdd_test2
	./cdd_test2

dplex_test1: dplex_test1.o $(DOFILES)
	$(CC) $(CFLAGS) $(DOFILES) dplex_test1.o -o dplex_test1

dplex_test2: dplex_test2.o $(DOFILES)
	$(CC) $(CFLAGS) $(DOFILES) dplex_test2.o -o dplex_test2

dplex_test3: dplex_test3.o $(DOFILES)
	$(CC) $(CFLAGS) $(DOFILES) dplex_test3.o -o dplex_test3

clean:
	rm -rf core *.out cdd cdd_test1 cdd_test2 dplex_test1 dplex_test2 dplex_test3 *.o *~ cddmathlinktm.* cddmathlink

all: cdd_test1 cdd_test2 dplex_test1 dplex_test2 dplex_test3

everything: ml all
