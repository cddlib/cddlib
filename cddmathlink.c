/* cddmathlink.c: Main test program to call the cdd library cddlib
   from Mathematica using MathLink.
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   March 13, 1999
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
#include "cddlib.h"
#include "cddmlio.h"
#include "mathlink.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

void allvertices(int m_input, int d_input, double *a_input)
/* output vertices and incidences */
{
  dd_PolyhedraPtr poly;
  dd_MatrixPtr A,G;
  dd_SetFamilyPtr GI;
  dd_rowrange i,m; 
  dd_colrange j,d;
  boolean found;

  m=(dd_rowrange)m_input; d=(dd_colrange)d_input;
  A=dd_CreateMatrix(m,d);
  for (i=0; i<m; i++){
    for (j=0; j<d; j++) A->matrix[i][j]=a_input[i*d+j];
  }
  dd_PolyhedraLoadMatrix(&poly, Inequality, A);
  found=dd_DoubleDescription(poly);
    /* compute the second (generator) representation */
  if (found) {
    G=dd_CopyGenerators(poly);
    GI=dd_CopyIncidence(poly);

    MLPutFunction(stdlink,"List",2);
    dd_MLWriteMatrix(G);
    dd_MLWriteSetFamily(GI);
  } else {
    dd_MLWriteError(poly);
  }

  dd_FreeMatrix(&A);
  dd_FreeMatrix(&G);
  dd_FreeSetFamily(&GI);
}

void allvertices2(int m_input, int d_input, double *a_input)
/* output vertices, incidences and adjacency */
{
  dd_PolyhedraPtr poly;
  dd_MatrixPtr A,G;
  dd_SetFamilyPtr GI,GA;
  dd_rowrange i,m; 
  dd_colrange j,d;
  boolean found;

  m=(dd_rowrange)m_input; d=(dd_colrange)d_input;
  printf("m=%d   d=%d\n",m, d);
  A=dd_CreateMatrix(m,d);
  for (i=0; i<m; i++){
    for (j=0; j<d; j++) A->matrix[i][j]=a_input[i*d+j];
  }
  dd_PolyhedraLoadMatrix(&poly, Inequality, A);
  found=dd_DoubleDescription(poly);
    /* compute the second (generator) representation */
  if (found){
    G=dd_CopyGenerators(poly);
    GI=dd_CopyIncidence(poly);
    GA=dd_CopyAdjacency(poly);

    MLPutFunction(stdlink,"List",3);
    dd_MLWriteMatrix(G);
    dd_MLWriteSetFamily(GI);
    dd_MLWriteSetFamily(GA);
  } else {
    dd_MLWriteError(poly);
  }
  dd_FreeMatrix(&A);
  dd_FreeMatrix(&G);
  dd_FreeSetFamily(&GI);
  dd_FreeSetFamily(&GA);
}

void allfacets(int n_input, int d_input, double *g_input)
/* output facets and incidences */
{
  dd_PolyhedraPtr poly;
  dd_MatrixPtr A,G;
  dd_SetFamilyPtr AI;
  dd_rowrange i,n; 
  dd_colrange j,d;
  boolean found;

  n=(dd_rowrange)n_input; d=(dd_colrange)d_input;
  G=dd_CreateMatrix(n,d);
  for (i=0; i<n; i++){
    for (j=0; j<d; j++) G->matrix[i][j]=g_input[i*d+j];
  }
  dd_PolyhedraLoadMatrix(&poly, Generator, G);
  found=dd_DoubleDescription(poly);
    /* compute the second (inequality) representation */
  if (found){
    A=dd_CopyInequalities(poly);
    AI=dd_CopyIncidence(poly);

    MLPutFunction(stdlink,"List",2);
    dd_MLWriteMatrix(A);
    dd_MLWriteSetFamily(AI);
  } else {
    dd_MLWriteError(poly);
  }

  dd_FreeMatrix(&A);
  dd_FreeMatrix(&G);
  dd_FreeSetFamily(&AI);
}


void allfacets2(int n_input, int d_input, double *g_input)
/* output facets, incidences and adjacency */
{
  dd_PolyhedraPtr poly;
  dd_MatrixPtr A,G;
  dd_SetFamilyPtr AI, AA;
  dd_rowrange i,n; 
  dd_colrange j,d;
  boolean found;

  n=(dd_rowrange)n_input; d=(dd_colrange)d_input;
  G=dd_CreateMatrix(n,d);
  for (i=0; i<n; i++){
    for (j=0; j<d; j++) G->matrix[i][j]=g_input[i*d+j];
  }
  dd_PolyhedraLoadMatrix(&poly, Generator, G);
  found=dd_DoubleDescription(poly);  
    /* compute the second (inequality) representation */
  if (found){
    A=dd_CopyInequalities(poly);
    AI=dd_CopyIncidence(poly);
    AA=dd_CopyAdjacency(poly);

    MLPutFunction(stdlink,"List",3);
    dd_MLWriteMatrix(A);
    dd_MLWriteSetFamily(AI);
    dd_MLWriteSetFamily(AA);
  } else {
    dd_MLWriteError(poly);
  }

  dd_FreeMatrix(&A);
  dd_FreeMatrix(&G);
  dd_FreeSetFamily(&AI);
  dd_FreeSetFamily(&AA);
}


int main(argc, argv)
        int argc; char* argv[];
{
        return MLMain(argc, argv);
}

/* end of cddmathlink.c */
