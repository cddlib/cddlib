/* cddmathlink.c: Main test program to call the cdd library cddlib
   from Mathematica using MathLink.
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.92dev, September 22, 2001
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
#include "mathlink.h"
#include "cddmlio.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

void allvertices(int m_input, int d_input, double *a_input)
/* output vertices and incidences */
{
  dd_PolyhedraPtr poly;
  dd_MatrixPtr A=NULL,G=NULL;
  dd_SetFamilyPtr GI=NULL;
  dd_rowrange i,m; 
  dd_colrange j,d;
  dd_ErrorType err;

  m=(dd_rowrange)m_input; d=(dd_colrange)d_input;
  A=dd_CreateMatrix(m,d);
  for (i=0; i<m; i++){
    for (j=0; j<d; j++) dd_set_d(A->matrix[i][j],a_input[i*d+j]);
  }
  A->representation=dd_Inequality;
  poly=dd_DDMatrix2Poly(A, &err);
    /* compute the second (generator) representation */
  if (err==dd_NoError) {
    G=dd_CopyGenerators(poly);
    GI=dd_CopyIncidence(poly);

    MLPutFunction(stdlink,"List",2);
    dd_MLWriteMatrix(G);
    dd_MLWriteSetFamily(GI);
  } else {
    dd_MLWriteError(poly);
  }

  dd_FreeMatrix(A);
  dd_FreeMatrix(G);
  dd_FreeSetFamily(GI);
}

void allvertices2(int m_input, int d_input, double *a_input)
/* output vertices, incidences and adjacency */
{
  dd_PolyhedraPtr poly;
  dd_MatrixPtr A=NULL,G=NULL;
  dd_SetFamilyPtr GI=NULL,GA=NULL;
  dd_SetFamilyPtr AI=NULL,AA=NULL;
  dd_rowrange i,m; 
  dd_colrange j,d;
  dd_ErrorType err;

  m=(dd_rowrange)m_input; d=(dd_colrange)d_input;
  A=dd_CreateMatrix(m,d);
  for (i=0; i<m; i++){
    for (j=0; j<d; j++) dd_set_d(A->matrix[i][j],a_input[i*d+j]);
  }
  A->representation=dd_Inequality;
  poly=dd_DDMatrix2Poly(A, &err);
    /* compute the second (generator) representation */
  if (err==dd_NoError){
    G=dd_CopyGenerators(poly);
    GI=dd_CopyIncidence(poly);
    GA=dd_CopyAdjacency(poly);
    AI=dd_CopyInputIncidence(poly);
    AA=dd_CopyInputAdjacency(poly);

    MLPutFunction(stdlink,"List",5);
    dd_MLWriteMatrix(G);
    dd_MLWriteSetFamily(GI);
    dd_MLWriteSetFamily(GA);
    dd_MLWriteSetFamily(AI);
    dd_MLWriteSetFamily(AA);
  } else {
    dd_MLWriteError(poly);
  }
  dd_FreeMatrix(A);
  dd_FreeMatrix(G);
  dd_FreeSetFamily(GI);
  dd_FreeSetFamily(GA);
  dd_FreeSetFamily(AI);
  dd_FreeSetFamily(AA);
}

void allfacets(int n_input, int d_input, double *g_input)
/* output facets and incidences */
{
  dd_PolyhedraPtr poly;
  dd_MatrixPtr A=NULL,G=NULL;
  dd_SetFamilyPtr AI=NULL;
  dd_rowrange i,n; 
  dd_colrange j,d;
  dd_ErrorType err;

  n=(dd_rowrange)n_input; d=(dd_colrange)d_input;
  G=dd_CreateMatrix(n,d);
  for (i=0; i<n; i++){
    for (j=0; j<d; j++) dd_set_d(G->matrix[i][j],g_input[i*d+j]);
  }
  G->representation=dd_Generator;
  poly=dd_DDMatrix2Poly(G, &err);
    /* compute the second (inequality) representation */
  if (err==dd_NoError){
    A=dd_CopyInequalities(poly);
    AI=dd_CopyIncidence(poly);

    MLPutFunction(stdlink,"List",2);
    dd_MLWriteMatrix(A);
    dd_MLWriteSetFamily(AI);
  } else {
    dd_MLWriteError(poly);
  }

  dd_FreeMatrix(A);
  dd_FreeMatrix(G);
  dd_FreeSetFamily(AI);
}


void allfacets2(int n_input, int d_input, double *g_input)
/* output facets, incidences and adjacency */
{
  dd_PolyhedraPtr poly;
  dd_MatrixPtr A=NULL,G=NULL;
  dd_SetFamilyPtr AI=NULL, AA=NULL;
  dd_SetFamilyPtr GI=NULL, GA=NULL;
  dd_rowrange i,n; 
  dd_colrange j,d;
  dd_ErrorType err;

  n=(dd_rowrange)n_input; d=(dd_colrange)d_input;
  G=dd_CreateMatrix(n,d);
  for (i=0; i<n; i++){
    for (j=0; j<d; j++) dd_set_d(G->matrix[i][j],g_input[i*d+j]);
  }
  G->representation=dd_Generator;
  poly=dd_DDMatrix2Poly(G, &err);
    /* compute the second (inequality) representation */
  if (err==dd_NoError){
    A=dd_CopyInequalities(poly);
    AI=dd_CopyIncidence(poly);
    AA=dd_CopyAdjacency(poly);
    GI=dd_CopyInputIncidence(poly);
    GA=dd_CopyInputAdjacency(poly);

    MLPutFunction(stdlink,"List",5);
    dd_MLWriteMatrix(A);
    dd_MLWriteSetFamily(AI);
    dd_MLWriteSetFamily(AA);
    dd_MLWriteSetFamily(GI);
    dd_MLWriteSetFamily(GA);
  } else {
    dd_MLWriteError(poly);
  }

  dd_FreeMatrix(A);
  dd_FreeMatrix(G);
  dd_FreeSetFamily(AI);
  dd_FreeSetFamily(AA);
  dd_FreeSetFamily(GI);
  dd_FreeSetFamily(GA);
}

#if MACINTOSH_MATHLINK

int main( int argc, char* argv[])
{
	/* Due to a bug in some standard C libraries that have shipped with
	 * MPW, zero is passed to MLMain below.  (If you build this program
	 * as an MPW tool, you can change the zero to argc.)
	 */
    dd_set_global_constants();  /* First, this must be called to use cddlib. */

	argc = argc; /* suppress warning */
	return MLMain( 0, argv);
}

#elif WINDOWS_MATHLINK

#if __BORLANDC__
#pragma argsused
#endif

int PASCAL WinMain( HINSTANCE hinstCurrent, HINSTANCE hinstPrevious, LPSTR lpszCmdLine, int nCmdShow)
{
	char  buff[512];
	char FAR * buff_start = buff;
	char FAR * argv[32];
	char FAR * FAR * argv_end = argv + 32;

    dd_set_global_constants();  /* First, this must be called to use cddlib. */

	hinstPrevious = hinstPrevious; /* suppress warning */

	if( !MLInitializeIcon( hinstCurrent, nCmdShow)) return 1;
	MLScanString( argv, &argv_end, &lpszCmdLine, &buff_start);
	return MLMain( argv_end - argv, argv);
}

#else

int main(argc, argv)
	int argc; char* argv[];
{
    dd_set_global_constants();  /* First, this must be called to use cddlib. */

	return MLMain(argc, argv);
}

#endif

/* end of cddmathlink.c */
