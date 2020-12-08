/* redcheck.c: Test program to call the cdd library cddlib
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
  dd_MatrixPtr M=NULL,M2=NULL;
  dd_ErrorType err=dd_NoError;
  dd_rowset redrows,linrows;
  mytype val;
  dd_DataFileType inputfile;
  FILE *reading=NULL;
  time_t starttime,endtime;

  dd_set_global_constants();  /* First, this must be called. */

  dd_init(val);
  if (argc>1) strcpy(inputfile,argv[1]);
  if (argc<=1 || !SetInputFile(&reading,argv[1])){
    dd_WriteProgramDescription(stdout);
    fprintf(stdout,"\ncddlib test program to check redundancy of an H/V-representation\nby Clarkson's algorithm.");
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

  time(&starttime);
  fprintf(stdout, "redundant rows: ");
  redrows=dd_RedundantRowsViaShooting(M, &err);
  set_fwrite(stdout, redrows);

  M2=dd_MatrixSubmatrix(M, redrows);

  fprintf(stdout, "Implicit linearity (after removal of redundant rows): ");
  linrows=dd_ImplicitLinearityRows(M2, &err);

  if (M->representation==dd_Generator)
    fprintf(stdout," %ld  ", set_card(linrows));
  else 
    fprintf(stdout," %ld  ", set_card(linrows));
  set_fwrite(stdout,linrows);
  set_uni(M2->linset, M2->linset, linrows); 
      /* add the implicit linrows to the given linearity rows */

  time(&endtime);
  printf("\nNonredundant representation (except for the linearity part):\n");
  dd_WriteMatrix(stdout, M2);
  dd_WriteTimes(stdout,starttime,endtime);

  dd_FreeMatrix(M);
  dd_FreeMatrix(M2);
  set_free(linrows);
  set_free(redrows);
_L99:;
  if (err!=dd_NoError) dd_WriteErrorMessages(stderr,err);
  return 0;
}


/* end of redcheck.c */
