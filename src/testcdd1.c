/* testcdd1.c: Main test program to call the cdd library cddlib
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

dd_boolean SetInputFile(FILE **f, char *fname)
{
  dd_boolean success=dd_FALSE;

  if ( ( *f = fopen(fname, "r") )!= NULL) {
    printf("input file %s is open\n", fname);
    success=dd_TRUE;
  }
  else{
    printf("The input file %s not found\n",fname);
  }
  return success;
}

int main(int argc, char *argv[])
{
  dd_PolyhedraPtr poly;
  dd_MatrixPtr M;
  dd_ErrorType err;
  dd_DataFileType inputfile;
  FILE *reading=NULL;
  dd_MatrixPtr A, G;
  dd_SetFamilyPtr GI,GA;

  dd_set_global_constants();  /* First, this must be called. */

  dd_SetInputFile(&reading,inputfile, &err);
  if (err==dd_NoError) {
    M=dd_PolyFile2Matrix(reading, &err);
  }
  else {
    printf("Input file not found\n");
    goto _L99;
  }

  if (err==dd_NoError) {
    poly=dd_DDMatrix2Poly(M, &err); /* compute the second representation */
    if (err!=dd_NoError) {
      dd_WriteErrorMessages(stdout,err);  goto _L99;
    }
    A=dd_CopyInequalities(poly);
    G=dd_CopyGenerators(poly);
    GI=dd_CopyIncidence(poly);
    GA=dd_CopyAdjacency(poly);

    if (poly->representation==dd_Inequality) {
      printf("\nInput is an H-representation\n");
    } else {
      printf("\nInput is a V-representation\n");
    }
    dd_WriteMatrix(stdout,A); printf("\n");
    dd_WriteMatrix(stdout,G);

    printf("\nHere is the incidence list:\n");
    dd_WriteSetFamily(stdout,GI);

    printf("\nHere is the adjacency list:\n");
    dd_WriteSetFamily(stdout,GA);

    dd_FreePolyhedra(poly);
    /* This is to remove all the space allocated for poly. */
    dd_FreeMatrix(M);
    dd_FreeMatrix(A);
    dd_FreeMatrix(G);
    dd_FreeSetFamily(GI);
    dd_FreeSetFamily(GA);
  } else {
    dd_WriteErrorMessages(stdout,err);
  }
_L99:
  dd_free_global_constants();  /* At the end, this must be called. */
  return 0;
}


/* end of testcdd1.c */
