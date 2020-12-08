/* testshoot.c: Main test program to call the cdd lp library
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
        
  dd_ErrorType err=dd_NoError;
  dd_LPSolverType solver=dd_DualSimplex; 
     /* either DualSimplex or CrissCross */
  dd_LPPtr lp,lp1;   
    /* pointer to LP data structure that is not visible by user. */
  dd_LPSolutionPtr lps1; 
    /* pointer to LP solution data that is visible by user. */

  dd_rowrange i,m;
  dd_colrange d;
  dd_NumberType numb;
  dd_MatrixPtr A;
  dd_Arow r;
  dd_colrange j;
  int iter;

/* Define an LP  */
/*
  max  0 + x1 + x2 
  s.t.
       1        -   x2  >= 0
       1 -   x1         >= 0
             x1         >= 0
                    x2  >= 0
       2 -   x1 -   x2  >= 0

  For this LP, we set up a matrix A as 4 x 3 matrix and a row vector:

      1   0  -1   <- 1st constraint
      1  -1   0
      0   1   0
      0   0   1 
      2  -1  -1   <- last constraint

      0   1   1   <- objective row
*/

  dd_set_global_constants();
    
  numb=dd_Real;   /* set a number type */
  m=5;    /* number of rows  */
  d=3;    /* number of columns */
  A=dd_CreateMatrix(m,d);
  dd_set_si(A->matrix[0][0],1); dd_set_si(A->matrix[0][1], 0); dd_set_si(A->matrix[0][2],-1);
  dd_set_si(A->matrix[1][0],1); dd_set_si(A->matrix[1][1],-1); dd_set_si(A->matrix[1][2], 0);
  dd_set_si(A->matrix[2][0],0); dd_set_si(A->matrix[2][1], 1); dd_set_si(A->matrix[2][2], 0);
  dd_set_si(A->matrix[3][0],0); dd_set_si(A->matrix[3][1], 0); dd_set_si(A->matrix[3][2], 1);
  dd_set_si(A->matrix[4][0],2); dd_set_si(A->matrix[4][1],-1); dd_set_si(A->matrix[4][2],-1);

  dd_set_si(A->rowvec[0],0);    dd_set_si(A->rowvec[1], 1); dd_set_si(A->rowvec[2], 1);

  A->objective=dd_LPmax;
  lp=dd_Matrix2LP(A, &err); /* load an LP */
  if (lp==NULL) goto _L99;

/* Find an interior point with cdd LP library. */
    printf("\n--- Running dd_FindInteriorPoint ---\n");
    lp1=dd_MakeLPforInteriorFinding(lp);
    dd_LPSolve(lp1,solver,&err);
    if (err!=dd_NoError) goto _L99;

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

/* Do shootings from the interior point. */

  dd_InitializeArow(d, &r);

  printf("Shooting test from the point:");
  dd_WriteArow(stdout,lps1->sol,d);  printf("\n");

  for (iter=1; iter<=3; iter++){
    dd_set_si(r[0],0); dd_set_si(r[1], 9998+iter); dd_set_si(r[2], 10000);
    printf("Shooting to the direction:");
    dd_WriteArow(stdout,r,d);  printf("\n");
    i=dd_RayShooting(A, lps1->sol, r);
    printf("The first hyperplane hit: %ld.\n\n", i);
  }

/* Free allocated spaces. */
  dd_FreeLPData(lp);
  dd_FreeLPSolution(lps1);
  dd_FreeLPData(lp1);
  dd_FreeArow(d, r);
  dd_FreeMatrix(A);

_L99:;
  if (err!=dd_NoError) dd_WriteErrorMessages(stdout, err);
  dd_free_global_constants();  /* At the end, this should be called. */
  return 0;
}

/* end of testlp3.c */

