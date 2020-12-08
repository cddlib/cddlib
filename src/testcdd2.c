/* testcdd2.c: Main test program to call the cdd library cddlib
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
*/

/*  This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#include "setoper.h"
#include "cdd.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

int main(int argc, char *argv[])
{
  dd_PolyhedraPtr poly;
  dd_MatrixPtr A, B, G;
  dd_rowrange m; 
  dd_colrange d;
  dd_ErrorType err;

  dd_set_global_constants();  /* First, this must be called to use cddlib. */

  m=4; d=3;
  A=dd_CreateMatrix(m,d);
  dd_set_si(A->matrix[0][0],7); dd_set_si(A->matrix[0][1],-3); dd_set_si(A->matrix[0][2], 0);
  dd_set_si(A->matrix[1][0],7); dd_set_si(A->matrix[1][1], 0); dd_set_si(A->matrix[1][2],-3);
  dd_set_si(A->matrix[2][0],1); dd_set_si(A->matrix[2][1], 1); dd_set_si(A->matrix[2][2], 0);
  dd_set_si(A->matrix[3][0],1); dd_set_si(A->matrix[3][1], 0); dd_set_si(A->matrix[3][2], 1);
  /* 7 - 3 x1          >= 0
     7         - 3x2   >= 0
     1 +   x1          >= 0
     1         +  x2   >= 0
  */
  A->representation=dd_Inequality;
  poly=dd_DDMatrix2Poly(A, &err);  /* compute the second (generator) representation */
  if (err!=dd_NoError) goto _L99;
  printf("\nInput is H-representation:\n");
  G=dd_CopyGenerators(poly);
  dd_WriteMatrix(stdout,A);  printf("\n");
  dd_WriteMatrix(stdout,G);
  dd_FreeMatrix(A);
  dd_FreeMatrix(G);

  /* Add inequalities:
     7 +  x1   -3x2   = 0
     7 - 3x1   + x2   >= 0
  */
  m=2;
  B=dd_CreateMatrix(m,d);
  dd_set_d(B->matrix[0][0],7.0); dd_set_d(B->matrix[0][1], 1.0); dd_set_d(B->matrix[0][2],-3.0);
  dd_set_d(B->matrix[1][0],7.0); dd_set_d(B->matrix[1][1],-3.0); dd_set_d(B->matrix[1][2], 1.0);
  set_addelem(B->linset,1); /* setting the first to be equality */

/* Above dd_set_d is used instead of dd_set_si.  This might be useful when your input is float. Yet,
   0.33333 won't be converted to 1/3 when -DGMPRATIONAL flag is used.  Better alternative might be
   dd_set_si2 function to assign a rational number, e.g.
   dd_set_si2(B->matrix[0][0],1,3).  Use these three assignment functions according to your need.
   These functions are defined in cddmp.h and cddmp.c.
*/

  dd_DDInputAppend(&poly,B, &err); /* append the two inequalities and compute the generators */
  if (err!=dd_NoError) goto _L99;
  A=dd_CopyInequalities(poly);  /* get the inequalities (=input). */
  G=dd_CopyGenerators(poly);  /* get the generators (=output). */
  printf("\nNew H-representation with added inequalities:\n");
  dd_WriteMatrix(stdout,A);  printf("\n");
  dd_WriteMatrix(stdout,G);

  dd_FreeMatrix(A);
  dd_FreeMatrix(B);
  dd_FreeMatrix(G);
  dd_FreePolyhedra(poly);

_L99:
  if (err!=dd_NoError) dd_WriteErrorMessages(stdout,err);

  dd_free_global_constants();  /* At the end, this must be called. */
  return 0;
}


/* end of testcdd2.c */
