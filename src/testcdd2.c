/* testcdd2.c: Main test program to call the cdd library cddlib
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.90, May 28, 2000
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
#include "cdd.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

int main(int argc, char *argv[])
{
  dd_PolyhedraPtr poly;
  dd_MatrixPtr A, G;
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
  A->representation=Inequality;
  poly=dd_Matrix2Poly(A, &err);
  dd_DoubleDescription(poly, &err);  /* compute the second (generator) representation */
  printf("\nInput is an H-representation:\n\n");
  G=dd_CopyGenerators(poly);
  dd_WriteMatrix(stdout,A);  printf("\n");
  dd_WriteMatrix(stdout,G);

  dd_FreeMatrix(A);
  dd_FreeMatrix(G);

  return 0;
}


/* end of testcdd2.c */
