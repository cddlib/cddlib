/* dplex_test2.c: Main test program to call the dplex library
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
  dp_LPSolutionPtr lps; /* pointer to LP solution data that is visible by user. */
  dp_rowrange i;
  dp_Arow a;

  while (error==dp_None) {

/* Input an LP using the dplex library  */
    lp=dp_LPInput(&reading, &error);
    if (error!=dp_None) goto _L99;

    SetWriteFile(&writing);

/* Solve the LP by dplex LP solver. */
    printf("\n--- Running dp_LPSolve ---\n");
    dp_LPSolve(lp, &error);
    lps=dp_LPSolutionLoad(lp);
    if (error!=dp_None) goto _L99;

    /* Write the LP solutions by dplex LP reporter. */
    dp_WriteLPResult(stdout, lp, error);
    dp_WriteLPResult(writing, lp, error);

/* Modify the LP by reversing an inequality. */
    lp1=dp_LPCopy(lp);  /* create a copy of lp */
    i=5;
    if (dp_LPReverseRow(lp1, i)){
      printf("\n--- Reverse the %ld-th inequality ---\n", i);
      dp_LPSolve(lp1, &error);
      if (error!=dp_None) goto _L99;
      dp_WriteLPResult(stdout, lp1, error);
      dp_WriteLPResult(writing, lp1, error);
    } else {
      printf("\n--- Reversing the %ld-th inequality failed ---\n", i);
    }
    dp_FreeLPData(&lp1);  /* remove lp1 */

/* Modify the LP by replacing an inequality. */
    lp1=dp_LPCopy(lp);  /* create a copy of lp */
    i=5;
    a=dp_LPCopyRow(lp1,i); /* copy the i-th constraint to a */
    a[0]=a[0]+1;  /* modify the vector a  by relaxing the rhs by 1 */

    if (dp_LPReplaceRow(lp1, i, a)){  /* replace the i-th constraint with a */
      printf("\n--- Replace the %ld-th inequality ---\n", i);
      dp_LPSolve(lp1, &error);
      if (error!=dp_None) goto _L99;
      dp_WriteLPResult(stdout, lp1, error);
      dp_WriteLPResult(writing, lp1, error);
    } else {
      printf("\n--- Replacing the %ld-th inequality failed ---\n", i);
    }
    dp_FreeLPData(&lp1);      /* remove lp1 */

/* Free allocated spaces. */
    dp_FreeLPSolution(&lps);  /* remove lp solution lps */
    dp_FreeLPData(&lp);       /* remove lp */
  }
_L99:;
  if (error!=dp_None) dp_WriteErrorMessages(stdout, error);
  return 0;
}

/* end of dplex_test2.c */
