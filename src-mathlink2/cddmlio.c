/* cddmlio.c: MathLink Basic Input and Output Procedures for cddlib.c
   written by Komei Fukuda, fukuda@cs.mcgill.ca     
   Version 0.93dev, Jan 15, 2003      
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
  char *str=NULL;

  if (A==NULL){
    rowmax=0; colmax=0;
  }
  MLPutFunction(stdlink,"List",rowmax);
  for (i=0; i < rowmax; i++) {
    MLPutFunction(stdlink,"List",colmax);
    for (j=0; j < colmax; j++) {
#if defined GMPRATIONAL
      str=dd_MLGetStrForNumber(A[i][j]);
      MLPutString(stdlink, str);
      if (str!=NULL) free(str);
#else
      a=dd_get_d(A[i][j]);
      MLPutDouble(stdlink, a);
#endif
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

char *dd_MLGetStrForNumber(mytype x)
{
  /* This is to make a string of rational expression for GMP rational number x.
     It does nothing (return NULL) if x is double.  */
     
  char *sd=NULL,*sn=NULL,*st=NULL;
#if defined GMPRATIONAL
  mpz_t zn,zd;
  int len;

  mpz_init(zn); mpz_init(zd);
  mpq_canonicalize(x);
  mpq_get_num(zn,x);
  mpq_get_den(zd,x);
  if (mpz_sgn(zn)==0){
    st=(char *)malloc(2);
    strcpy(st," 0");
  } else if (mpz_cmp_ui(zd,1U)==0){
    st=mpz_get_str(st,10,zn);
  } else {
    sn=mpz_get_str(sn,10,zn);sd=mpz_get_str(sd,10,zd);
    len=strlen(sn)+strlen(sd)+2;
    st=(char *)malloc(len);
    strcpy(st,sn);
    strcat(st,"/");strcat(st,sd);
  }
  mpz_clear(zn); mpz_clear(zd);
  if (sd!=NULL) free(sd);
  if (sn!=NULL) free(sn);
#else
  /* do nothing */
#endif
  /* printf("String for Number =%s\n",st);  */
  return st;
}

void dd_MLSetMatrixWithString(dd_rowrange m, dd_colrange d, char line[], dd_MatrixPtr M)
{
   dd_rowrange i=0;
   dd_colrange j=0;
   char *next,*copy;
   const char ct[]=", {}\n";  /* allows separators ","," ","{","}". */
   mytype value;
   double dval;

   dd_init(value);
   copy=(char *)malloc(strlen(line) + 4048); /* some extra space as buffer.  Somehow linux version needs this */
   strcpy(copy,line);
   next=strtok(copy,ct);
   i=0; j=0;
   do {
#if defined GMPRATIONAL
      dd_sread_rational_value(next, value);
#else
      dval=atof(next);
      dd_set_d(value,dval);
#endif
      dd_WriteNumber(stderr,value);
      dd_set(M->matrix[i][j],value);
      j++;
      if (j == d) {
         i++; j=0;
      }
  } while ((next=strtok(NULL,ct))!=NULL);
  free(copy);
  dd_clear(value);
  return;
}




