/* projection.c: Test program to call the cdd library cddlib
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.93, July 18,  2003
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
  dd_MatrixPtr M=NULL,M1=NULL,M2=NULL;
  dd_colrange j,s,t,d;
  dd_ErrorType err=dd_NoError;
  dd_rowset redrows,linrows;
  dd_colset delset;
  mytype val;
  dd_DataFileType inputfile;
  FILE *reading=NULL;

  dd_set_global_constants();  /* First, this must be called. */

  dd_init(val);
  if (argc>1) strcpy(inputfile,argv[1]);
  if (argc<=1 || !SetInputFile(&reading,argv[1])){
    dd_WriteProgramDescription(stdout);
    fprintf(stdout,"\ncddlib test program to apply the Block Elimination to an H-polyhedron.\n");
    dd_SetInputFile(&reading,inputfile, &err);
  }
  if (err==dd_NoError) {
    M=dd_PolyFile2Matrix(reading, &err);
  }
  else {
    fprintf(stderr,"Input file not found\n");
    goto _L99;
  }

  if (err!=dd_NoError) goto _L99;

  d=M->colsize;
  set_initialize(&delset, d);

  printf("How many variables to elminate? (max %ld): ",d-1);
  scanf("%ld",&s);

  for (j=1; j<=s; j++){
    printf("\n%ld th deletion variable): ",j);
    scanf("%ld",&t);
    set_addelem(delset, t+1);
  }
  
  M1=dd_BlockElimination(M, delset, &err);

  dd_WriteMatrix(stdout, M1);

  fprintf(stdout, "redundant rows: ");
  redrows=dd_RedundantRowsViaShooting(M1, &err);
  set_fwrite(stdout, redrows);

  M2=dd_MatrixSubmatrix(M1, redrows);

  fprintf(stdout, "Implicit linearity (after removal of redundant rows): ");
  linrows=dd_ImplicitLinearityRows(M2, &err);

  if (M->representation==dd_Generator)
    fprintf(stdout," %ld  ", set_card(linrows));
  else 
    fprintf(stdout," %ld  ", set_card(linrows));
  set_fwrite(stdout,linrows);
  set_uni(M2->linset, M2->linset, linrows); 

  printf("\nNonredundant representation (except for the linearity part):\n");
  dd_WriteMatrix(stdout, M2);

  dd_FreeMatrix(M);
  dd_FreeMatrix(M1);
  dd_FreeMatrix(M2);
  set_free(delset);
  set_free(redrows);
  set_free(linrows);

_L99:;
  /* if (err!=dd_NoError) dd_WriteErrorMessages(stderr,err); */
  return 0;
}


/* end of projection.c */
