/* testlp1.c: Main test program to call the cdd lp library
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
  dd_PolyhedraPtr poly;
  dd_DataFileType inputfile;
  int ans;

  dd_set_global_constants(); /* First, this must be called once to use cddlib. */

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
    printf("\nShould I generate all extreme points (y/n)? ");
    ans=getchar();
    if (ans=='y' || ans=='Y'){
      poly=dd_DDMatrix2Poly(M, &error);
      G=dd_CopyGenerators(poly);
      printf("\nAll the vertices of the feasible region.\n");
      dd_WriteMatrix(stdout,G);

      /* Free allocated spaces. */
      dd_FreeMatrix(G);
      dd_FreePolyhedra(poly);
    }

/* Free allocated spaces. */
    dd_FreeMatrix(M);
    dd_FreeLPData(lp);
  }
_L99:;
  fclose(reading);
  if (error!=dd_NoError) dd_WriteErrorMessages(stdout, error);
  return 0;
}

/* end of testlp1.c */
