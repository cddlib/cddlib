/* cdd_test2.c: Main test program to call the cdd library cddlib
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.80, March 2, 1999
   Standard ftp site: ftp.ifor.math.ethz.ch, Directory: pub/fukuda/cdd
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
#include "cddlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

int main(int argc, char *argv[])
{
  dd_PolyhedraPtr poly;
  dd_MatrixPtr A, A1, Acopy, A1copy, G, G1;
  dd_rowrange m,m1; 
  dd_colrange d;

  m=5; d=4;
  A=dd_CreateMatrix(m,d);
  A->matrix[0][0]= 3; A->matrix[0][1]=-1; A->matrix[0][2]= 0; A->matrix[0][3]= 0;
  A->matrix[1][0]= 3; A->matrix[1][1]= 0; A->matrix[1][2]=-1; A->matrix[1][3]= 0;
  A->matrix[2][0]= 1; A->matrix[2][1]= 1; A->matrix[2][2]= 0; A->matrix[2][3]= 0;
  A->matrix[3][0]= 1; A->matrix[3][1]= 0; A->matrix[3][2]= 1; A->matrix[3][3]= 0;
  A->matrix[4][0]= 1; A->matrix[4][1]= 0; A->matrix[4][2]= 0; A->matrix[4][3]= 1;
  /* 3 - x1          >= 0
     3      -x2      >= 0
     1 + x1          >= 0
     1      +x2      >= 0
     1          +x3  >= 0
  */
  dd_PolyhedraLoadMatrix(&poly, Inequality, A);
  dd_DoubleDescription(poly);  /* compute the second (generator) representation */
  printf("\nInput is an H-representation:\n");
  Acopy=dd_CopyInequalities(poly);
  G=dd_CopyGenerators(poly);
  printf("\n\nH-representation\n");  dd_WriteMatrix(stdout,Acopy);
  printf("\n\nV-representation\n");  dd_WriteMatrix(stdout,G);

  printf("\nNow, we add three new inequalities, and recompute the generators:\n");
  m1=3; d=4;
  A1=dd_CreateMatrix(m1,d);  /* create a new matrix with etra inequalities */
  A1->matrix[0][0]= 2; A1->matrix[0][1]= 0; A1->matrix[0][2]= 0; A1->matrix[0][3]=-1;
  A1->matrix[1][0]= 8; A1->matrix[1][1]=-1; A1->matrix[1][2]=-1; A1->matrix[1][3]=-1;
  A1->matrix[2][0]=-2; A1->matrix[2][1]= 1; A1->matrix[2][2]= 1; A1->matrix[2][3]= 1;
  /* 2          -x3  >= 0
     8  -x1 -x2 -x3  >= 0
    -2  +x1 +x2 +x3  >= 0
  */

  dd_DDAddInequalities(poly, A1);
  A1copy=dd_CopyInequalities(poly);
  G1=dd_CopyGenerators(poly);
  printf("\n\nH-representation\n");  dd_WriteMatrix(stdout,A1copy);
  printf("\n\nV-representation\n");  dd_WriteMatrix(stdout,G1);

  dd_FreeMatrix(&A);
  dd_FreeMatrix(&Acopy);
  dd_FreeMatrix(&G);
  dd_FreeMatrix(&A1copy);
  dd_FreeMatrix(&G1);

  return 0;
}


/* end of cdd_test2.c */
