/* dplex_test3.c: Main test program to call the dplex library
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
#include "dplex.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

FILE *reading, *writing;


void SetWriteFile(FILE **f)
{
  char *fname;
  fname="dplex_test.out";
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
        
  dp_ErrorType error=dp_None;
  dp_LPPtr lp;   
    /* pointer to LP data structure that is not visible by user. */
  dp_LPSolutionPtr lps; 
    /* pointer to LP solution data that is visible by user. */

  dp_rowrange M;
  dp_colrange N;
  dp_NumberType NUMB;
  dp_LPConversionType CONV;
  dp_Amatrix A;
  dp_rowrange OBJ; 
  dp_colrange RHS,j;


/* Define an LP  */
/*
  max  0 + 3 x1 + 4 x2 + 2 x3
  s.t.
       4 - 2 x1               >= 0
       8 -   x1        - 2 x3 >= 0
       6        - 3 x2 -   x3 >= 0
             x1               >= 0
                    x2        >= 0
                           x3 >= 0

  For this LP, we set up a matrix A as 7 x 4 matrix:

      4  -2   0   0   <- 1st constraint
      8  -1   0  -2
      6   0  -3  -1
      0   1   0   0
      0   0   1   0
      0   0   0   1   <- last constraint
      0   3   4   2   <- objective row
*/

    
  NUMB=dp_Real;   /* set a number type */
  M=7;    /* number of rows */
  N=4;    /* number of columns */
  dp_InitializeAmatrix(M,N,&A);  /* create a space for A */
  OBJ=M;   /* objective row number */
  RHS=1;  /* right hand side column number */
  CONV=dp_LPmax;
  A[0][0]= 4; A[0][1]=-2; A[0][2]= 0; A[0][3]= 0;
  A[1][0]= 8; A[1][1]=-1; A[1][2]= 0; A[1][3]=-2;
  A[2][0]= 6; A[2][1]= 0; A[2][2]=-3; A[2][3]=-1;
  A[3][0]= 0; A[3][1]= 1; A[3][2]= 0; A[3][3]= 0;
  A[4][0]= 0; A[4][1]= 0; A[4][2]= 1; A[4][3]= 0;
  A[5][0]= 0; A[5][1]= 0; A[5][2]= 0; A[5][3]= 1;
  A[6][0]= 0; A[6][1]= 3; A[6][2]= 4; A[6][3]= 2;

  lp=dp_LPLoad(CONV,NUMB,M,N,A,OBJ,RHS,&error); /* load an LP */
  if (error!=dp_None) goto _L99;

/* Solve the LP by dplex LP solver. */
  printf("\n--- Running dp_LPSolve ---\n");
  dp_LPSolve(lp, &error);  /* Solve the LP */
  if (error!=dp_None) goto _L99;

/* Write the LP solutions by dplex LP reporter. */
/*  dp_WriteLPResult(stdout, lp, error); */
/*  dp_WriteLPResult(writing, lp, error); */

/* One can access the solutions by loading them.  See dp_WriteLPResult
   for outputing the results correctly. */
  lps=dp_LPSolutionLoad(lp);
  if (lps->LPS==dp_Optimal){
    printf("Optimal solution found:\n");
    printf("  primal_solution\n");
    for (j=1; j<lps->d; j++) {
      printf("  %3ld : ",j);
      dp_WriteReal(stdout,lps->sol[j]);
      printf("\n");
    }
    printf("  dual_solution\n");
    for (j=1; j<lps->d; j++){
      if (lps->nbindex[j+1]>0) {
        printf("  %3ld : ",lps->nbindex[j+1]);
        dp_WriteReal(stdout,lps->dsol[j]); printf("\n");
      }
    }
    printf("  optimal_value : % .9E\n",lps->optvalue);
  }


/* Free allocated spaces. */
  dp_FreeLPSolution(&lps);
  dp_FreeLPData(&lp);

_L99:;
  if (error!=dp_None) dp_WriteErrorMessages(stdout, error);
  return 0;
}

/* end of dplex_test3.c */

/* The dual LP is

/*
  min  0 + 4 y1 + 8 y2 + 6 y3
  s.t.
      -3 + 2 x1 +   y2        >= 0
      -4               + 3 y3 >= 0
      -2        + 2 y2 +   y3 >= 0
             y1               >= 0
                    y2        >= 0
                           y3 >= 0
*/
