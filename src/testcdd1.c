/* testcdd1.c: Main test program to call the cdd library cddlib
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

boolean SetInputFile(FILE **f, char *fname)
{
  boolean success=FALSE;
  success=FALSE;

  if ( ( *f = fopen(fname, "r") )!= NULL) {
    printf("input file %s is open\n", fname);
    success=TRUE;
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
  boolean found=FALSE;
  dd_ErrorType err;
  dd_DataFileType inputfile;
  FILE *reading=NULL;
  dd_MatrixPtr A, G;
  dd_SetFamilyPtr GI,GA;

  dd_set_global_constants();  /* First, this must be called. */

  dd_SetInputFile(&reading,inputfile, &err);
  if (err==NoError) {
    M=dd_PolyFile2Matrix(reading, &err);
  }
  else {
    printf("Input file not found\n");
    goto _L99;
  }

  if (err==NoError) {
    poly=dd_Matrix2Poly(M, &err);
    found=dd_DoubleDescription(poly,&err);  /* compute the second representation */
    if (!found) {
      dd_WriteErrorMessages(stdout,err);  goto _L99;
    }
    A=dd_CopyInequalities(poly);
    G=dd_CopyGenerators(poly);
    GI=dd_CopyIncidence(poly);
    GA=dd_CopyAdjacency(poly);

    if (poly->representation==Inequality) {
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
    /* This is to remove the workspace created by dd_DoubleDescription;
       PolyhedraData poly won't be erased. */
    dd_FreeMatrix(A);
    dd_FreeMatrix(G);
    dd_FreeSetFamily(GI);
    dd_FreeSetFamily(GA);
  } else {
    dd_WriteErrorMessages(stdout,err);
  }
_L99:
  return 0;
}


/* end of testcdd1.c */
