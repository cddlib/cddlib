/* cddmlio.c: MathLink Basic Input and Output Procedures for cddlib.c
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.92dev, Sept. 22, 2001
*/

/* cdd.c : C-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#include "setoper.h"  /* set operation library header (Ver. March 16,1995 or later) */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "mathlink.h"
#include "cddmlio.h"

void dd_MLWriteAmatrix(dd_Amatrix A, long rowmax, long colmax)
{
  long i,j;
  double a;

  if (A==NULL){
    rowmax=0; colmax=0;
  }
  MLPutFunction(stdlink,"List",rowmax);
  for (i=0; i < rowmax; i++) {
    MLPutFunction(stdlink,"List",colmax);
    for (j=0; j < colmax; j++) {
      a=dd_get_d(A[i][j]);
      MLPutDouble(stdlink, a);
    }
  }
}

void dd_MLWriteMatrix(dd_MatrixPtr M)
{
  MLPutFunction(stdlink,"List",2);
  dd_MLWriteAmatrix(M->matrix, M->rowsize, M->colsize);
  dd_MLWriteSet(M->linset);
}

void dd_MLWriteSet(set_type S)
{
  long j;

  MLPutFunction(stdlink,"List",set_card(S));
  for (j=1; j <= S[0]; j++) {
    if (set_member(j, S)) MLPutLongInteger(stdlink, j);
  }
}

void dd_MLWriteSetFamily(dd_SetFamilyPtr F)
{
  long i,j;

  if (F!=NULL){
    MLPutFunction(stdlink,"List",F->famsize);
    for (i=0; i < F->famsize; i++) {
      MLPutFunction(stdlink,"List",set_card(F->set[i]));
      for (j=1; j <= F->setsize; j++) {
        if (set_member(j, F->set[i])) MLPutLongInteger(stdlink, j);
      }
    }
  }
}

void dd_MLWriteError(dd_PolyhedraPtr poly)
{
  MLPutFunction(stdlink,"List",3);
  MLPutFunction(stdlink,"List",0);
  MLPutFunction(stdlink,"List",1);
  MLPutString(stdlink,"Error occured: code");
  MLPutFunction(stdlink,"List",1);
  MLPutInteger(stdlink,poly->child->Error);
}

