/* projection.c: Test program to call the cdd library cddlib
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
  dd_MatrixPtr M=NULL,M1=NULL;
  dd_colrange j,s,t,d;
  dd_ErrorType err=dd_NoError;
  dd_rowset redset,impl_linset;
  dd_colset delset;
  dd_rowindex newpos;
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

  dd_MatrixCanonicalize(&M1,&impl_linset,&redset,&newpos,&err);

  if (err!=dd_NoError) goto _L99;

  fprintf(stdout, "\nRedundant rows: ");
  set_fwrite(stdout, redset);
  fprintf(stdout, "\n");

  dd_WriteMatrix(stdout, M1);

  dd_FreeMatrix(M);
  dd_FreeMatrix(M1);
  set_free(delset);
  set_free(redset);
  set_free(impl_linset);
  free(newpos);

_L99:;
  /* if (err!=dd_NoError) dd_WriteErrorMessages(stderr,err); */
  dd_free_global_constants();  /* At the end, this should be called. */
  return 0;
}


/* end of projection.c */
