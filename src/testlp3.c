/* testlp3.c: Main test program to call the cdd lp library
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.90c, June 12, 2000
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

FILE *reading, *writing;


void SetWriteFile(FILE **f)
{
  char *fname;
  fname="testlp.out";
  *f = fopen(fname, "w");
  printf("file %s is open\n",fname);
}


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
        
  dd_ErrorType err=NoError;
  dd_LPSolverType solver=DualSimplex; 
     /* either DualSimplex or CrissCross */
  dd_LPPtr lp,lp1;   
    /* pointer to LP data structure that is not visible by user. */
  dd_LPSolutionPtr lps,lps1; 
    /* pointer to LP solution data that is visible by user. */

  dd_rowrange m;
  dd_colrange n;
  dd_NumberType numb;
  dd_MatrixPtr A;
  dd_colrange j;

/* Define an LP  */
/*
  max  0 + 3 x1 + 4 x2 
  s.t.
       4 - 2 x1 -   x2  >= 0
       2        -   x2  >= 0
             x1         >= 0
                    x2  >= 0

  For this LP, we set up a matrix A as 4 x 3 matrix and a row vector:

      4  -2  -1   <- 1st constraint
      2   0  -1
      0   1   0
      0   0   1   <- last constraint

      0   3   4   <- objective row
*/

  dd_set_global_constants();
    
  numb=Real;   /* set a number type */
  m=4;    /* number of rows  */
  n=3;    /* number of columns */
  A=dd_CreateMatrix(m,n);
  dd_set_si(A->matrix[0][0],4); dd_set_si(A->matrix[0][1],-2); dd_set_si(A->matrix[0][2],-1);
  dd_set_si(A->matrix[1][0],2); dd_set_si(A->matrix[1][1], 0); dd_set_si(A->matrix[1][2],-1);
  dd_set_si(A->matrix[2][0],0); dd_set_si(A->matrix[2][1], 1); dd_set_si(A->matrix[2][2], 0);
  dd_set_si(A->matrix[3][0],0); dd_set_si(A->matrix[3][1], 0); dd_set_si(A->matrix[3][2], 1);

  dd_set_si(A->rowvec[0],0);    dd_set_si(A->rowvec[1], 3); dd_set_si(A->rowvec[2], 4);

  A->objective=LPmax;
  lp=dd_Matrix2LP(A, &err); /* load an LP */
  if (lp==NULL) goto _L99;

/* Solve the LP by cdd LP solver. */
  printf("\n--- Running dd_LPSolve ---\n");
  solver=DualSimplex;
  dd_LPSolve(lp, solver, &err);  /* Solve the LP */
  if (err!=NoError) goto _L99;

/* Write the LP solutions by cdd LP reporter. */
/*  dd_WriteLPResult(stdout, lp, err); */
/*  dd_WriteLPResult(writing, lp, err); */

/* One can access the solutions by loading them.  See dd_WriteLPResult
   for outputing the results correctly. */
  lps=dd_CopyLPSolution(lp);
  if (lps->LPS==Optimal){
    printf("Optimal solution found:\n");
    printf("  primal_solution\n");
    for (j=1; j<lps->d; j++) {
      printf("  %3ld : ",j);
      dd_WriteNumber(stdout,lps->sol[j]);
      printf("\n");
    }
    printf("  dual_solution\n");
    for (j=1; j<lps->d; j++){
      if (lps->nbindex[j+1]>0) {
        printf("  %3ld : ",lps->nbindex[j+1]);
        dd_WriteNumber(stdout,lps->dsol[j]); printf("\n");
      }
    }
    printf("  optimal_value : "); dd_WriteNumber(stdout,lps->optvalue);
    printf("\n");
  }

/* Find an interior point with cdd LP library. */
    printf("\n--- Running dd_FindInteriorPoint ---\n");
    lp1=dd_MakeLPforInteriorFinding(lp);
    dd_LPSolve(lp1,solver,&err);
    if (err!=NoError) goto _L99;

    /* Write an interior point. */
    lps1=dd_CopyLPSolution(lp1);
    if (dd_Positive(lps1->optvalue)){
      printf("An interior point found: (");
      for (j=1; j <(lps1->d)-1; j++) {
        dd_WriteNumber(stdout,lps1->sol[j]);
      }
      printf(")\n");
    }
    if (dd_Negative(lps1->optvalue)) 
      printf("The feasible region is empty.\n");
    if (dd_EqualToZero(lps1->optvalue)) 
      printf("The feasible region is nonempty but has no interior point.\n");

/* Free allocated spaces. */
  dd_FreeLPSolution(lps);
  dd_FreeLPData(lp);
  dd_FreeLPSolution(lps1);
  dd_FreeLPData(lp1);

_L99:;
  if (err!=NoError) dd_WriteErrorMessages(stdout, err);
  return 0;
}

/* end of testlp3.c */

