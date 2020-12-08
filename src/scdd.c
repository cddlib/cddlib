/* scdd.c: Main test program to call the cdd library cddlib
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

dd_boolean SetInputFile(FILE **f, dd_DataFileType fname)
{
  dd_boolean success=dd_FALSE;
  success=dd_FALSE;

  if ( ( *f = fopen(fname, "r") )!= NULL) {
    printf("input file %s is open\n", fname);
    success=dd_TRUE;
  }
  else{
    printf("The input file %s not found\n",fname);
  }
  return success;
}

dd_boolean SetWriteFile(FILE **f, dd_DataFileType fname)
{
  dd_boolean success=dd_FALSE;

  if ( (*f = fopen(fname, "w")) != NULL){
    printf("output file %s is open\n",fname);
    success=dd_TRUE;
  }
  else{
    printf("The output file %s cannot be opened\n",fname);
  }
  return success;
}


int main(int argc, char *argv[])
{
  dd_PolyhedraPtr poly;
  dd_LPPtr lp;
  dd_MatrixPtr M,A;
  dd_ErrorType err=dd_NoError;
  dd_DataFileType inputfile,outputfile;
  FILE *reading=NULL, *writing;

  dd_set_global_constants();  /* First, this must be called. */
  dd_log=dd_TRUE;  /* output log */

  if (argc>1) strcpy(inputfile,argv[1]);
  if (argc<=1 || !SetInputFile(&reading,argv[1])){
    dd_WriteProgramDescription(stdout);
    dd_SetInputFile(&reading,inputfile, &err);
  }
  if (err==dd_NoError) {
    M=dd_PolyFile2Matrix(reading, &err);
  }
  else {
    printf("Input file not found\n");
    goto _L99;
  }

  if (err!=dd_NoError) goto _L99;

  if (M->objective==dd_LPnone){ /* do representation conversion */
    poly=dd_DDMatrix2Poly2(M, dd_LexMin, &err);
    /* equivalent to poly=dd_DDMatrix2Poly2(M, &err) when the second argument is set to dd_LexMin. */
    if (err!=dd_NoError) goto _L99;

    dd_SetWriteFileName(inputfile, outputfile, 'o', poly->representation);
    SetWriteFile(&writing, outputfile);
    dd_WriteProgramDescription(writing);
    dd_WriteRunningMode(writing, poly);
    switch (poly->representation) {
    case dd_Inequality:
      fprintf(writing, "ext_file: Generators\n");
      A=dd_CopyGenerators(poly);
      dd_WriteMatrix(writing,A);
      dd_FreeMatrix(A);
      break;

    case dd_Generator:
      fprintf(writing, "ine_file: Inequalities\n");
      A=dd_CopyInequalities(poly);
      dd_WriteMatrix(writing,A);
      dd_FreeMatrix(A);
      break;

    default:
      break;
    }
    dd_WriteDDTimes(writing,poly);
    fclose(writing);

    dd_SetWriteFileName(inputfile, outputfile, 'a', poly->representation);
    SetWriteFile(&writing, outputfile);
    dd_WriteAdjacency(writing,poly);
    fclose(writing);

    dd_SetWriteFileName(inputfile, outputfile, 'j', poly->representation);
    SetWriteFile(&writing, outputfile);
    dd_WriteInputAdjacency(writing,poly);
    fclose(writing);

    dd_SetWriteFileName(inputfile, outputfile, 'i', poly->representation);
    SetWriteFile(&writing, outputfile);
    dd_WriteIncidence(writing,poly);
    fclose(writing);

    dd_SetWriteFileName(inputfile, outputfile, 'n', poly->representation);
    SetWriteFile(&writing, outputfile);
    dd_WriteInputIncidence(writing,poly);
    fclose(writing);

    dd_FreeMatrix(M);
    dd_FreePolyhedra(poly);

  } else { /* solve the LP */
    lp=dd_Matrix2LP(M, &err);  if (err!=dd_NoError) goto _L99;
    dd_LPSolve(lp,dd_DualSimplex,&err);  if (err!=dd_NoError) goto _L99;

    dd_SetWriteFileName(inputfile, outputfile, 's', M->representation);
    SetWriteFile(&writing, outputfile);
    dd_WriteLPResult(writing, lp, err);
    fclose(writing);

    dd_FreeMatrix(M);
    dd_FreeLPData(lp);
  }
_L99:
  if (err!=dd_NoError) dd_WriteErrorMessages(stdout,err);
  return 0;
}


/* end of simplecdd.c */
