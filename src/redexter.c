/* redexter.c: Test program to call the cdd library cddlib
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
  dd_MatrixPtr M1=NULL,M2=NULL,M2row=NULL,M1plus=NULL;
  dd_colrange d1;
  dd_rowrange i,m1,m2,m1plus;
  dd_ErrorType err=dd_NoError,err1=dd_NoError,err2=dd_NoError;
  dd_rowset delset,rowset2;
  dd_Arow cvec; /* certificate */  

  time_t starttime, endtime;
  dd_DataFileType inputfile1,inputfile2;
  FILE *reading1=NULL,*reading2=NULL;

  dd_set_global_constants();  /* First, this must be called. */

  dd_WriteProgramDescription(stdout);
  fprintf(stdout,"\ncddlib test program to check redundancy of additional data.\n");
  if (argc>2){
    strcpy(inputfile1,argv[1]);
    strcpy(inputfile2,argv[2]);
  }
  /* 
  if (argc<=2){
    fprintf(stdout,"\nUsage:\n   redexter file1 file2\n");
	goto _L99;
  }
  */
  if (!SetInputFile(&reading1,argv[1])){
    fprintf(stdout,"\nSpecify file1.\n");
    dd_SetInputFile(&reading1,inputfile1, &err1);
  }
  if (!SetInputFile(&reading2,argv[2])){
    fprintf(stdout,"\nSpecify the secondary file.\n");
    dd_SetInputFile(&reading2,inputfile2, &err2);
  }
  if ((err1==dd_NoError) && (err2==dd_NoError)) {
    M1=dd_PolyFile2Matrix(reading1, &err1);
    M2=dd_PolyFile2Matrix(reading2, &err2);
  }
  else {
    fprintf(stderr,"Input file(s) not found\n");
    goto _L99;
  }

  if ((err1!=dd_NoError) || (err2!=dd_NoError)) goto _L99;

  m1=M1->rowsize;
  m2=M2->rowsize;
  set_initialize(&delset,m2);
  m1plus=m1+1;
  if (M1->representation==dd_Generator){
    d1=(M1->colsize)+1;
  } else {
    d1=M1->colsize;
  }
  dd_InitializeArow(d1,&cvec);

  fprintf(stdout, "\nThe first matrix\n");
  dd_WriteMatrix(stdout, M1);
  fprintf(stdout, "\nThe second matrix\n");
  dd_WriteMatrix(stdout, M2);
  
  printf("\nChecking whether each row of the second matrix is redundant w.r.t. the first.\n");

  time(&starttime);

  for (i=1; i<=m2; i++){
    set_initialize(&rowset2,m2);
	set_addelem(rowset2, i);
    set_compl(delset, rowset2);
    M2row=dd_MatrixSubmatrix(M2, delset);
	M1plus=dd_MatrixAppend(M1,M2row); 
	
    if (dd_Redundant(M1plus, m1plus, cvec, &err)) {
	  printf("%ld-th row: redundant\n", i);
	} else {
	  printf("%ld-th row: non-redundant\n A certificate:", i);
	  dd_WriteArow(stdout, cvec, d1);
	}

    dd_FreeMatrix(M1plus);
	dd_FreeMatrix(M2row);
    set_free(rowset2);
  }

  time(&endtime);

  dd_WriteTimes(stdout,starttime,endtime);

  set_free(delset);
  dd_FreeMatrix(M1);
  dd_FreeMatrix(M2);

_L99:;
  if (err1!=dd_NoError) dd_WriteErrorMessages(stderr,err1);
  if (err2!=dd_NoError) dd_WriteErrorMessages(stderr,err2);
  return 0;
}


/* end of redexter.c */
