/* testlp2.c: Main test program to call the cdd lp library
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

FILE *reading, *writing;

int main(int argc, char *argv[])
{
  /* The original LP data m x n matrix 
     = | b   -A  |
       | c0  c^T |,
   
  where the LP to be solved is to
  maximize  c^T x  +   c0
  subj. to
            A   x  <=  b.
  */
        
  dd_ErrorType error=dd_NoError;
  dd_LPSolverType solver;  /* either DualSimplex or CrissCross */
  dd_LPPtr lp;   

  dd_rowrange m;
  dd_colrange n;
  dd_NumberType numb;
  dd_MatrixPtr A;
  dd_ErrorType err;


/* Define an LP  */
/*
  max  0 + 3 x1 + 4 x2 
  s.t.
       4/3 - 2 x1 -   x2  >= 0
       2/3        -   x2  >= 0
               x1         >= 0
                      x2  >= 0

  For this LP, we set up a matrix A as 4 x 3 matrix and a row vector:

      4/3  -2  -1   <- 1st constraint
      2/3   0  -1
      0     1   0
      0     0   1   <- last constraint

      0     3   4   <- objective row
*/

  dd_set_global_constants();
    
  numb=dd_Real;   /* set a number type */
  m=4;    /* number of rows  */
  n=3;    /* number of columns */
  A=dd_CreateMatrix(m,n);
  dd_set_si2(A->matrix[0][0],4,3); dd_set_si(A->matrix[0][1],-2); dd_set_si(A->matrix[0][2],-1);
  dd_set_si2(A->matrix[1][0],2,3); dd_set_si(A->matrix[1][1], 0); dd_set_si(A->matrix[1][2],-1);
  dd_set_si(A->matrix[2][0],0); dd_set_si(A->matrix[2][1], 1); dd_set_si(A->matrix[2][2], 0);
  dd_set_si(A->matrix[3][0],0); dd_set_si(A->matrix[3][1], 0); dd_set_si(A->matrix[3][2], 1);

  dd_set_si(A->rowvec[0],0);    dd_set_si(A->rowvec[1], 3); dd_set_si(A->rowvec[2], 4);

  A->objective=dd_LPmax;
  lp=dd_Matrix2LP(A, &err); /* load an LP */
  if (lp==NULL) goto _L99;

/* Print the LP. */
  printf("\n--- LP to be solved  ---\n");
  dd_WriteLP(stdout, lp);

/* Solve the LP by cdd LP solver. */
  printf("\n--- Running dd_LPSolve ---\n");
  solver=dd_DualSimplex;
  dd_LPSolve(lp, solver, &error);  /* Solve the LP */
  if (error!=dd_NoError) goto _L99;

/* Write the LP solutions by cdd LP reporter. */
  dd_WriteLPResult(stdout, lp, error);

/* Free allocated spaces. */
  dd_FreeLPData(lp);
  dd_FreeMatrix(A);

_L99:;
  if (error!=dd_NoError) dd_WriteErrorMessages(stdout, error);
  dd_free_global_constants();  /* At the end, this should be called. */
  return 0;
}

/* end of testlp2.c */

/* The dual LP is

  min  0 + 4 y1 + 2 y2
  s.t.
      -3 + 2 y1        >= 0
      -4     y1 +   y2 >= 0
             y1        >= 0
                    y2 >= 0
*/
