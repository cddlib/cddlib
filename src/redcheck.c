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
  dd_MatrixPtr M=NULL;
  dd_rowrange i,m;
  dd_ErrorType err=dd_NoError;
  dd_rowindex newpos;
  dd_rowset impl_linset,redset;
  time_t starttime, endtime;
  dd_DataFileType inputfile;
  FILE *reading=NULL;

  dd_set_global_constants();  /* First, this must be called. */

  if (argc>1) strcpy(inputfile,argv[1]);
  if (argc<=1 || !SetInputFile(&reading,argv[1])){
    dd_WriteProgramDescription(stdout);
    fprintf(stdout,"\ncddlib test program to check redundancy of an H/V-representation.\n");
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

  m=M->rowsize;
  fprintf(stdout, "Canonicalize the matrix.\n");
    
  time(&starttime);
  dd_MatrixCanonicalize(&M, &impl_linset, &redset, &newpos, &err);
  time(&endtime);
  
  if (err!=dd_NoError) goto _L99;

  fprintf(stdout, "Implicit linearity rows are: "); set_fwrite(stdout, impl_linset);

  fprintf(stdout, "\nRedundant rows are: "); set_fwrite(stdout, redset);
  fprintf(stdout, "\n");
  
  fprintf(stdout, "Nonredundant representation:\n");
  fprintf(stdout, "The new row positions are as follows (orig:new).\nEach redundant row has the new number 0.\nEach deleted duplicated row has a number nagative of the row that\nrepresents its equivalence class.\n");
  
  for (i=1; i<=m; i++){
   fprintf(stdout, " %ld:%ld",i, newpos[i]); 
  }
  fprintf(stdout, "\n");
  dd_WriteMatrix(stdout, M);
  
  dd_WriteTimes(stdout,starttime,endtime);

  set_free(redset);
  set_free(impl_linset);
  dd_FreeMatrix(M);
  free(newpos);

_L99:;
  if (err!=dd_NoError) dd_WriteErrorMessages(stderr,err);
  return 0;
}


/* end of redcheck.c */
