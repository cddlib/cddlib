/* cddlib.c: cdd library  (library version of cdd)
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.90d, June 25, 2000
   Standard ftp site: ftp.ifor.math.ethz.ch, Directory: pub/fukuda/cdd
*/

/* cdd : C-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.
   Please read COPYING (GNU General Public Licence) and
   the manual cddlibman.tex for detail.
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

/* The first version C0.21 was created on November 10,1993 
   with Dave Gillespie's p2c translator 
   from the Pascal program pdd.p written by Komei Fukuda. 
*/

#include "setoper.h" 
  /* set operation library header (Ver. June 1, 2000 or later) */
#include "cdd.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

/* Global Variables */
boolean debug               =FALSE;
/* GLOBAL CONSTANTS (to be set by dd_set_global_constants() */
mytype dd_zero;
mytype dd_one;
mytype dd_purezero;
mytype dd_minuszero;

/* #include <profile.h>    THINK C PROFILER */
/* #include <console.h>    THINK C PROFILER */

void DDInit(dd_ConePtr cone)
{
  cone->Error=NoError;
  cone->CompStatus=InProgress;
  cone->RayCount = 0;
  cone->TotalRayCount = 0;
  cone->FeasibleRayCount = 0;
  cone->WeaklyFeasibleRayCount = 0;
  cone->EdgeCount=0; /* active edge count */
  cone->TotalEdgeCount=0; /* active edge count */
  SetInequalitySets(cone);
  ComputeRowOrderVector(cone);
  cone->RecomputeRowOrder=FALSE;
}

void DDMain(dd_ConePtr cone)
{
  dd_rowrange hh, itemp, otemp;
  boolean localdebug=TRUE;

  if (cone->d<=0){
    cone->Iteration=cone->m;
    cone->FeasibleRayCount=0;
    cone->CompStatus=AllFound;
    goto _L99;
  }
  if (localdebug) {
     fprintf(stdout,"(Initially added rows ) = ");
     set_fwrite(stdout,cone->InitialHalfspaces);
  }
  while (cone->Iteration <= cone->m) {
    SelectNextHalfspace(cone, cone->WeaklyAddedHalfspaces, &hh);
    if (set_member(hh,cone->NonequalitySet)){  /* Skip the row hh */
      if (debug) {
        printf("*The row # %3ld should be inactive and thus skipped.\n", hh);
      }
      set_addelem(cone->WeaklyAddedHalfspaces, hh);
    } else {
      if (cone->PreOrderedRun)
        AddNewHalfspace2(cone, hh);
      else{
        AddNewHalfspace1(cone, hh);
      }
      set_addelem(cone->AddedHalfspaces, hh);
      set_addelem(cone->WeaklyAddedHalfspaces, hh);
    }
    if (!cone->PreOrderedRun){
      for (itemp=1; cone->OrderVector[itemp]!=hh; itemp++);
        otemp=cone->OrderVector[cone->Iteration];
      cone->OrderVector[cone->Iteration]=hh;
        /* store the dynamic ordering in ordervec */
      cone->OrderVector[itemp]=otemp;
        /* store the dynamic ordering in ordervec */
    }
    if (localdebug){
      printf("(Iter, Row, #Total, #Curr, #Feas)= %5ld %5ld %9ld %6ld %6ld\n",
        cone->Iteration, hh, cone->TotalRayCount, cone->RayCount,
        cone->FeasibleRayCount);
    }
    if (cone->CompStatus==AllFound||cone->CompStatus==RegionEmpty) {
      set_addelem(cone->AddedHalfspaces, hh);
      goto _L99;
    }
    (cone->Iteration)++;
  }
  _L99:;
  if (cone->d<=0 || cone->newcol[1]==0){ /* fixing the number of output */
     cone->parent->n=cone->LinearityDim + cone->FeasibleRayCount -1;
     cone->parent->ldim=cone->LinearityDim - 1;
  } else {
    cone->parent->n=cone->LinearityDim + cone->FeasibleRayCount;
    cone->parent->ldim=cone->LinearityDim;
  }
}


void dd_InitialDataSetup(dd_ConePtr cone)
{
  long j, r;
  dd_rowset ZSet;
  static dd_Arow Vector1,Vector2;
  static dd_colrange last_d=0;

  if (last_d < cone->d){
    if (last_d>0) {
    for (j=0; j<last_d; j++){
      dd_init(Vector1[j]);
      dd_init(Vector2[j]);
    }
    free(Vector1); free(Vector2);
    }
    Vector1=(mytype*)calloc(cone->d,sizeof(mytype));
    Vector2=(mytype*)calloc(cone->d,sizeof(mytype));
    for (j=0; j<cone->d; j++){
      dd_init(Vector1[j]);
      dd_init(Vector2[j]);
    }
    last_d=cone->d;
  }

  cone->RecomputeRowOrder=FALSE;
  cone->ArtificialRay = NULL;
  cone->FirstRay = NULL;
  cone->LastRay = NULL;
  set_initialize(&ZSet,cone->m);
  AddArtificialRay(cone);
  set_copy(cone->AddedHalfspaces, cone->InitialHalfspaces);
  set_copy(cone->WeaklyAddedHalfspaces, cone->InitialHalfspaces);
  UpdateRowOrderVector(cone, cone->InitialHalfspaces);
  for (r = 1; r <= cone->d; r++) {
    for (j = 0; j < cone->d; j++){
      dd_set(Vector1[j], cone->B[j][r-1]);
      dd_neg(Vector2[j], cone->B[j][r-1]);
    }
    Normalize(cone->d, Vector1);
    Normalize(cone->d, Vector2);
    ZeroIndexSet(cone->m, cone->d, cone->A, Vector1, ZSet);
    if (set_subset(cone->EqualitySet, ZSet)){
      if (debug) {
        printf("add an initial ray with zero set:");
        set_write(ZSet);
      }
      AddRay(cone, Vector1);
      if (cone->InitialRayIndex[r]==0) {
        AddRay(cone, Vector2);
        if (debug) {
          printf("and add its negative also.\n");
        }
      }
    }
  }
  CreateInitialEdges(cone);
  cone->Iteration = cone->d + 1;
  if (cone->Iteration >= cone->m) cone->CompStatus=AllFound;
  set_free(ZSet);
}

boolean DoubleDescription(dd_PolyhedraPtr poly, dd_ErrorType *err)
{
  dd_ConePtr cone=NULL;
  boolean found=FALSE;

  *err=NoError;
  if (poly!=NULL && (poly->child==NULL || poly->child->CompStatus!=AllFound)){
    cone=ConeDataLoad(poly);
    /* create a cone associated with poly by homogenization */
    time(&cone->starttime);
    DDInit(cone);
    if (poly->representation==Generator && poly->m<=0){
       *err=EmptyVrepresentation;
       cone->Error=*err;
    } else {
      FindInitialRays(cone, &found);
      if (found) {
        dd_InitialDataSetup(cone);
        DDMain(cone);
      }
    }
    time(&cone->endtime);
  }
  return found;
}

boolean dd_DDInputAppend(dd_PolyhedraPtr *poly, dd_MatrixPtr M,
  dd_ErrorType *err)
{
  /* This is imcomplete.  It simply solves the problem from scratch.  */
  boolean found;
  dd_ErrorType error;

  if ((*poly)->child!=NULL) dd_FreeDDMemory(*poly);
  AppendMatrix2Poly(poly, M);
  (*poly)->representation=Inequality;
  found=DoubleDescription(*poly, &error);
  *err=error;
  return found;
}

boolean dd_DDFile2File(char *ifile, char *ofile, dd_ErrorType *err)
{
  /* The representation conversion from an input file to an outfile.  */
  boolean found;
  dd_ErrorType error;
  FILE *reading=NULL,*writing=NULL;
  dd_PolyhedraPtr poly;
  dd_MatrixPtr M, A, G;

  if ( ( reading = fopen(ifile, "r") )!= NULL) {
    printf("input file %s is open\n", ifile);
    found=TRUE;
  }
  else{
    printf("The input file %s not found\n",ifile);
    found=FALSE;
    error=IFileNotFound;
    goto _L99;
  }

  if (found && (writing = fopen(ofile, "w")) != NULL){
    printf("output file %s is open\n",ofile);
    found=TRUE;
  }
  else{
    printf("The output file %s cannot be opened\n",ofile);
    found=FALSE;
    error=OFileNotOpen;
    goto _L99;
  }

  M=dd_PolyFile2Matrix(reading, &error);
  poly=dd_DDMatrix2Poly(M, &error);  /* compute the second representation */
  dd_FreeMatrix(M);

  if (error==NoError) {
    A=dd_CopyInequalities(poly);
    G=dd_CopyGenerators(poly);

    if (poly->representation==Inequality) {
      dd_WriteMatrix(writing,G);
    } else {
      dd_WriteMatrix(writing,A);
    }

    dd_FreePolyhedra(poly);
    dd_FreeMatrix(A);
    dd_FreeMatrix(G);
  } 

_L99:  *err=error;
  if (reading!=NULL) fclose(reading);
  if (writing!=NULL) fclose(writing);
  return found;
}


/* end of cddlib.c */
