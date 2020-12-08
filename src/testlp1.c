/* testlp1.c: Main test program to call the cdd lp library
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
  /* The original LP data  m x n matrix 
     = | b   -A  |
       | c0  c^T |,
   
  where the LP to be solved is to
  maximize  c^T x  +   c0
  subj. to
            A   x  <=  b.
  */
        
  dd_ErrorType error=dd_NoError;
  dd_MatrixPtr M,G;
  dd_LPSolverType solver=dd_DualSimplex;  /* either DualSimplex or CrissCross */
  dd_LPPtr lp;   /* pointer to LP data structure that is not visible by user. */
  dd_LPSolutionPtr  lps1;
  dd_colrange j;
  dd_rowset ImL, Lbasis;

  dd_PolyhedraPtr poly;
  dd_DataFileType inputfile;
  int ans;

  dd_set_global_constants(); /* First, this must be called once to use cddlib. */

  printf("Welcome to cddlib %s\n",dd_DDVERSION);

  while (error==dd_NoError) {

/* Input an LP using the cdd library  */
    dd_SetInputFile(&reading,inputfile,&error);
    if (error!=dd_NoError) goto _L99;
    M=dd_PolyFile2Matrix(reading, &error);
    if (error!=dd_NoError) goto _L99;
    /* dd_WriteMatrix(stdout, M);  */
    lp=dd_Matrix2LP(M, &error);
    if (error!=dd_NoError) goto _L99;
	

/* Solve the LP by cdd LP solver. */
    printf("\n--- Running dd_LPSolve ---\n");
    dd_LPSolve(lp,solver,&error);
    if (error!=dd_NoError) goto _L99;

    /* Write the LP solutions by cdd LP reporter. */
    dd_WriteLPResult(stdout, lp, error);

/* Generate all vertices of the feasible reagion */
    printf("\nDo you want to compute the generator representation (y/n)? ");
    ans=getchar();
    if (ans=='y' || ans=='Y'){
      poly=dd_DDMatrix2Poly(M, &error);
      G=dd_CopyGenerators(poly);
      printf("\nGenerators (All the vertices of the feasible region if bounded.)\n");
      dd_WriteMatrix(stdout,G);

      /* Free allocated spaces. */
      dd_FreeMatrix(G);
      dd_FreePolyhedra(poly);
    }

/* Find an interior point with cdd LP library. */
    printf("\nDo you want to find a relative interior point (y/n)? ");
    ans=getchar(); ans=getchar();
    if (ans=='y' || ans=='Y'){
      printf("\n--- Running dd_FindRelativeInteriorPoint ---\n");
	  dd_FindRelativeInterior(M, &ImL, &Lbasis, &lps1, &error);
      if (error!=dd_NoError) goto _L99;

      /* Write an interior point. */
      if (dd_Positive(lps1->optvalue)){
        printf("A relative interior point found: (");
        for (j=1; j <(lps1->d)-1; j++) {
          dd_WriteNumber(stdout,lps1->sol[j]);
        }
        printf(")\nThe dimension of the region = ");
		printf("%ld\n",M->colsize-set_card(Lbasis)-1);
		if (set_card(ImL)>0) {
		  printf("Implicit equations: "); set_write(ImL); printf("\n");
		}
      } else {
        printf("The feasible region is empty.\n");
	  }
	  dd_FreeLPSolution(lps1);
	  set_free(ImL);
	  set_free(Lbasis);
    }

/* Free allocated spaces. */
    dd_FreeMatrix(M);
    dd_FreeLPData(lp);
  }
_L99:;
  fclose(reading);
  if (error!=dd_NoError) dd_WriteErrorMessages(stdout, error);
  dd_free_global_constants();  /* At the end, this should be called. */
  return 0;
}

/* end of testlp1.c */
