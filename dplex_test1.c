/* dplex_test1.c: Main test program to call the dplex library
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
  /* The original LP data  m x n matrix 
     = | b   -A  |
       | c0  c^T |,
   
  where the LP to be solved is to
  maximize  c^T x  +   c0
  subj. to
            A   x  <=  b.
  */
        
  dp_ErrorType error=dp_None;
  dp_LPPtr lp,lp1;   /* pointer to LP data structure that is not visible by user. */
  dp_LPSolutionPtr lps1; /* pointer to LP solution data that is visible by user. */
  dp_colrange j;

  while (error==dp_None) {

/* Input an LP using the dplex library  */
    lp=dp_LPInput(&reading, &error);
    if (error!=dp_None) goto _L99;

    SetWriteFile(&writing);

/* Solve the LP by dplex LP solver. */
    printf("\n--- Running dp_LPSolve ---\n");
    dp_LPSolve(lp, &error);
    if (error!=dp_None) goto _L99;

    /* Write the LP solutions by dplex LP reporter. */
    dp_WriteLPResult(stdout, lp, error);
    dp_WriteLPResult(writing, lp, error);

/* Find an interior point with dplex. */
    printf("\n--- Running dp_FindInteriorPoint ---\n");
    lp1=dp_MakeLPforInteriorFinding(lp);
    dp_LPSolve(lp1, &error);
    if (error!=dp_None) goto _L99;

    /* Write an interior point. */
    lps1=dp_LPSolutionLoad(lp1);
    if (dp_Positive(lps1->optvalue)){
      printf("An interior point found: (");
      for (j=1; j <(lps1->d)-1; j++) dp_WriteReal(stdout,lps1->sol[j]);
      printf(")\n");
    }
    if (dp_Negative(lps1->optvalue)) 
      printf("The feasible region is empty.\n");
    if (dp_Zero(lps1->optvalue)) 
      printf("The feasible region is nonempty but has no interior point.\n");


/* Free allocated spaces. */
    dp_FreeLPSolution(&lps1);
    dp_FreeLPData(&lp);
    dp_FreeLPData(&lp1);
  }
_L99:;
  if (error!=dp_None) dp_WriteErrorMessages(stdout, error);
  return 0;
}

/* end of dplex_test1.c */
