/* cdd_test1.c: Main test program to call the cdd library cddlib
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.85, October 3, 1999
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

void SetWriteFile(FILE **f)
{
  char *fname;
  fname="cdd_test.out";
  *f = fopen(fname, "w");
  printf("file %s is open\n",fname);
}

int main(int argc, char *argv[])
{
  dd_PolyhedraPtr poly, poly1;
  dd_ErrorType err;
  FILE *writing=NULL;
  dd_MatrixPtr A, G, A1, G1;
  dd_SetFamilyPtr GI,GA;
  boolean found;

  dd_PolyhedraInput(&poly, &err);

  if (err==None) {
    SetWriteFile(&writing);
    found=dd_DoubleDescription(poly,&err);  /* compute the second representation */
    if (!found) {
      dd_WriteErrorMessages(stdout,err);  goto _L99;
    }
    A=dd_CopyInequalities(poly);
    G=dd_CopyGenerators(poly);
    GI=dd_CopyIncidence(poly);
    GA=dd_CopyAdjacency(poly);

    if (poly->Representation==Inequality) {
      printf("\nInput is an H-representation\n");
    } else {
      printf("\nInput is a V-representation\n");
    }
    printf("\nH-representation\n");  dd_WriteMatrix(stdout,A);
    printf("\nV-representation\n");  dd_WriteMatrix(stdout,G);

    fprintf(writing, "\nH-representation\n");  dd_WriteMatrix(writing,A);
    fprintf(writing, "\nV-representation\n");  dd_WriteMatrix(writing,G);

    printf("\nHere is the incidence list:\n");
    dd_WriteSetFamily(stdout,GI);

    printf("\nHere is the adjacency list:\n");
    dd_WriteSetFamilyWithNumbers(stdout,GA);

    printf("\nLet's try the reverse operation!\n");
    if (poly->Representation==Inequality){
      dd_PolyhedraLoadMatrix(&poly1, Generator, G);
      printf("Now, input is a V-representation (computed previously)\n");
    } else {
      dd_PolyhedraLoadMatrix(&poly1, Inequality, A);
      printf("Now, input is an H-representation (computed previously)\n");
    }
    found=dd_DoubleDescription(poly1,&err);  /* compute the second representation */
    if (!found) {
      dd_WriteErrorMessages(stdout,err);  goto _L99;
    }
    A1=dd_CopyInequalities(poly1);
    G1=dd_CopyGenerators(poly1);
    printf("\nH-representation\n");  dd_WriteMatrix(stdout,A1);
    printf("\nV-representation\n");  dd_WriteMatrix(stdout,G1);

    dd_FreeDDMemory(poly);  dd_FreeDDMemory(poly1);
    /* This is to remove the workspace created by dd_DoubleDescription;
       PolyhedraData poly won't be erased. */
    dd_FreeMatrix(&A);
    dd_FreeMatrix(&G);
    dd_FreeMatrix(&A1);
    dd_FreeMatrix(&G1);
    if (writing!=NULL) fclose(writing);
  } else {
    dd_WriteErrorMessages(stdout,err);
  }
_L99:
  return 0;
}


/* end of cdd_test1.c */
