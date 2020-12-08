/* adjacency.c: Test program to call the cdd library cddlib
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
  dd_MatrixPtr M=NULL,M2=NULL,M3=NULL;
  dd_SetFamilyPtr A=NULL;
  dd_colrange d;
  dd_ErrorType err=dd_NoError;
  dd_rowset redrows,linrows,ignoredrows, basisrows;
  dd_colset ignoredcols, basiscols;
  long rank;
  mytype val;
  time_t starttime, endtime;
  dd_DataFileType inputfile;
  FILE *reading=NULL;

  dd_set_global_constants();  /* First, this must be called. */

  dd_init(val);
  if (argc>1) strcpy(inputfile,argv[1]);
  if (argc<=1 || !SetInputFile(&reading,argv[1])){
    dd_WriteProgramDescription(stdout);
    fprintf(stdout,"\ncddlib test program to remove redundancy and compute adjacency of the resulting representation.\n");
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

  if (M->representation==dd_Generator) d=M->colsize+1; else d=M->colsize;

  fprintf(stdout, "redundant rows:\n");
  time(&starttime);
  redrows=dd_RedundantRows(M, &err);
  time(&endtime);
  set_fwrite(stdout, redrows);
  dd_WriteTimes(stdout,starttime,endtime);

  M2=dd_MatrixSubmatrix(M, redrows);

  fprintf(stdout, "Implicit linearity (after removal of redundant rows): ");
  linrows=dd_ImplicitLinearityRows(M2, &err);

  if (M->representation==dd_Generator)
    fprintf(stdout," %ld  ", set_card(linrows));
  else 
    fprintf(stdout," %ld  ", set_card(linrows));
  set_fwrite(stdout,linrows);
  set_uni(M2->linset, M2->linset, linrows); 
      /* add the implicit linrows to the explicit linearity rows */

  printf("\nNonredundant representation (except possibly for the linearity part):\n");
  dd_WriteMatrix(stdout, M2);

  /* To remove redundancy of the linearity part, 
     we need to compute the rank and a basis of the linearity part. */
  set_initialize(&ignoredrows, M2->rowsize);
  set_initialize(&ignoredcols, M2->colsize);
  set_compl(ignoredrows, M2->linset);
  rank=dd_MatrixRank(M2,ignoredrows,ignoredcols, &basisrows, &basiscols);
  set_diff(ignoredrows, M2->linset, basisrows);
  M3=dd_MatrixSubmatrix(M2, ignoredrows);
  if (rank>0){
    if (set_card(ignoredrows)) {
       fprintf(stdout,"\nThe following %ld linearity rows are dependent and unnecessary:", set_card(ignoredrows));
       set_fwrite(stdout,ignoredrows);
    }
  } else
    fprintf(stdout,"\nThe linearity rows are independent and thus minimal\n");  

  printf("Nonredundant representation (= minimal representation):\n");
  dd_WriteMatrix(stdout, M3);

  printf("\nAdjacency of the minimal representation:\n");
  A=dd_Matrix2Adjacency(M3, &err);
  dd_WriteSetFamily(stdout, A);
  
  dd_clear(val);
  set_free(linrows);
  set_free(basisrows);
  set_free(basiscols);
  set_free(ignoredrows);
  set_free(ignoredcols);
  set_free(redrows);
  
  if (A!=NULL) dd_FreeSetFamily(A);
  dd_FreeMatrix(M);
  dd_FreeMatrix(M2);
  dd_FreeMatrix(M3);

_L99:;
  if (err!=dd_NoError) dd_WriteErrorMessages(stderr,err);
  dd_free_global_constants();  /* At the end, this should be called. */
  return 0;
}


/* end of adjacency.c */
