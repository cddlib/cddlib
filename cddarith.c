/* cddarith.c:  Floating Point Arithmetic Procedures for cddlib.c
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.85, October 3, 1999
*/

/* cdd.c : C-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#include "setoper.h"  /* set operation library header (Ver. March 16,1995 or later) */
#include "cdd.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

void CheckAdjacency(dd_ConePtr cone,
    dd_RayPtr *RP1, dd_RayPtr *RP2, boolean *adjacent)
{
  dd_RayPtr TempRay;
  boolean localdebug=FALSE;
  static dd_rowset Face, Face1;
  static dd_rowrange last_m=0;
  
  if (last_m!=cone->m) {
    if (last_m>0){
      set_free(Face); set_free(Face1);
    }
    set_initialize(&Face, cone->m); 
    set_initialize(&Face1, cone->m); 
    last_m=cone->m;
  }

  if (debug) localdebug=TRUE;
  *adjacent = TRUE;
  set_int(Face1, (*RP1)->ZeroSet, (*RP2)->ZeroSet);
  set_int(Face, Face1, cone->AddedHalfspaces);
  if (set_card(Face)< cone->d - 2) {
    *adjacent = FALSE;
    if (localdebug) {
      printf("non adjacent: set_card(face) %ld < %ld = cone->d.\n",
        set_card(Face),cone->d);
    }
    return;
  }
  else if (cone->parent->NondegAssumed) {
  	*adjacent = TRUE;
  	return;
  }
  TempRay = cone->FirstRay;
  while (TempRay != NULL && *adjacent) {
    if (TempRay != *RP1 && TempRay != *RP2) {
    	set_int(Face1, TempRay->ZeroSet, cone->AddedHalfspaces);
      	if (set_subset(Face, Face1)) *adjacent = FALSE;
    }
    TempRay = TempRay->Next;
  }
}

void Eliminate(dd_ConePtr cone, dd_RayPtr*Ptr)
{
  /*eliminate the record pointed by Ptr^.Next*/
  dd_RayPtr TempPtr;

  TempPtr = (*Ptr)->Next;
  (*Ptr)->Next = (*Ptr)->Next->Next;
  if (TempPtr == cone->FirstRay)   /*Update the first pointer*/
    cone->FirstRay = (*Ptr)->Next;
  if (TempPtr == cone->LastRay)   /*Update the last pointer*/
    cone->LastRay = *Ptr;
  free(TempPtr->Ray);          /* free the ray vector memory */
  set_free(TempPtr->ZeroSet);  /* free the ZeroSet memory */
  free(TempPtr);   /* free the dd_Ray structure memory */
  cone->RayCount--; 
}

void SetInequalitySets(dd_ConePtr cone)
{
  dd_rowrange i;
  
  set_emptyset(cone->GroundSet);
  set_emptyset(cone->EqualitySet);
  set_emptyset(cone->NonequalitySet);  
  for (i = 1; i <= (cone->m); i++){
    set_addelem(cone->GroundSet, i);
    if (cone->parent->EqualityIndex[i]==1) set_addelem(cone->EqualitySet,i);
    if (cone->parent->EqualityIndex[i]==-1) set_addelem(cone->NonequalitySet,i);
  }
}


double AValue(dd_colrange d_size, dd_Amatrix A, double *p, dd_rowrange i)
{
  /*return the ith component of the vector  A x p */
  dd_colrange j;
  double temp;

  temp = 0.0;
  for (j = 0; j < d_size; j++)
    temp += A[i - 1][j] * p[j];
  return temp;
}

void StoreRay1(dd_ConePtr cone, double *p, boolean *feasible)
{  /* Original ray storing routine when RelaxedEnumeration is FALSE */
  dd_rowrange i,k,fii=cone->m+1;
  dd_colrange j;
  double temp;
  dd_RayPtr RR;
  boolean localdebug=FALSE;

  RR=cone->LastRay;
  if (debug) localdebug=TRUE;
  *feasible = TRUE;
  set_initialize(&(RR->ZeroSet),cone->m);
  RR->ARay = 0.0;
  for (j = 0; j < cone->d; j++)
    RR->Ray[j] = p[j];
  for (i = 1; i <= cone->m; i++) {
    k=cone->OrderVector[i];
    temp = AValue(cone->d, cone->A, p, k);
    if (dd_Zero(temp)) 
      set_addelem(RR->ZeroSet, k);
    if (dd_Negative(temp)){
      *feasible = FALSE;
      if (fii>cone->m) fii=i;  /* the first violating inequality index */
    }
  }
  RR->FirstInfeasIndex=fii;
  RR->feasible = *feasible;
}

void StoreRay2(dd_ConePtr cone, double *p, 
    boolean *feasible, boolean *weaklyfeasible)
   /* Ray storing routine when RelaxedEnumeration is TRUE.
       weaklyfeasible is true iff it is feasible with
       the strict_inequality conditions deleted. */
{
  dd_RayPtr RR;
  dd_rowrange i,k,fii=cone->m+1;
  dd_colrange j;
  double temp;
  boolean localdebug=FALSE;

  RR=cone->LastRay;
  if (debug) localdebug=TRUE;
  *feasible = TRUE;
  *weaklyfeasible = TRUE;
  set_initialize(&(RR->ZeroSet),cone->m);
  RR->ARay = 0.0;
  for (j = 0; j < cone->d; j++)
    RR->Ray[j] = p[j];
  for (i = 1; i <= cone->m; i++) {
    k=cone->OrderVector[i];
    temp = AValue(cone->d, cone->A, p, k);
/*    if (fabs(temp) < zero){ */
    if (dd_Zero(temp)){
      set_addelem(RR->ZeroSet, k);
      if (cone->parent->EqualityIndex[k]==-1) 
        *feasible=FALSE;  /* strict inequality required */
    }
/*    if (temp < -zero){ */
    if (dd_Negative(temp)){
      *feasible = FALSE;
      if (fii>cone->m && cone->parent->EqualityIndex[k]>=0) {
        fii=i;  /* the first violating inequality index */
        *weaklyfeasible=FALSE;
      }
    }
  }
  RR->FirstInfeasIndex=fii;
  RR->feasible = *feasible;
}


void AddRay(dd_ConePtr cone, double *p)
{  
  boolean feasible, weaklyfeasible;

  if (cone->FirstRay == NULL) {
    cone->FirstRay = (dd_RayPtr) malloc(sizeof(dd_Ray));
    cone->FirstRay->Ray = (double *) calloc(cone->d, sizeof(double));
    if (debug)
      printf("Create the first ray pointer\n");
    cone->LastRay = cone->FirstRay;
    cone->ArtificialRay->Next = cone->FirstRay;
  } else {
    cone->LastRay->Next = (dd_RayPtr) malloc(sizeof(dd_Ray));
    cone->LastRay->Next->Ray = (double *) calloc(cone->d, sizeof(double));
    if (debug) printf("Create a new ray pointer\n");
    cone->LastRay = cone->LastRay->Next;
  }
  cone->LastRay->Next = NULL;
  cone->RayCount++;
  cone->TotalRayCount++;
  if (debug) {
    if (cone->TotalRayCount % 100 == 0) {
      printf("*Rays (Total, Currently Active, Feasible) =%8ld%8ld%8ld\n",
	 cone->TotalRayCount, cone->RayCount, cone->FeasibleRayCount);
    }
  }
  if (cone->parent->RelaxedEnumeration){
    StoreRay2(cone, p, &feasible, &weaklyfeasible);
    if (weaklyfeasible) (cone->WeaklyFeasibleRayCount)++;
  } else {
    StoreRay1(cone, p, &feasible);
    if (feasible) (cone->WeaklyFeasibleRayCount)++;
    /* weaklyfeasible is equiv. to feasible in this case. */
  }
  if (!feasible) return;
  else {
    (cone->FeasibleRayCount)++;
  }
}

void AddArtificialRay(dd_ConePtr cone)
{  
  dd_Arow zerovector;
  long j;
  boolean feasible;

  dd_InitializeArow(cone->d, &zerovector);
  if (cone->ArtificialRay != NULL) {
    printf("Warning !!!  FirstRay in not nil.  Illegal Call\n");
    return;
  }
  cone->ArtificialRay = (dd_RayPtr) malloc(sizeof(dd_Ray));
  cone->ArtificialRay->Ray = (double *) calloc(cone->d, sizeof(double));
  if (debug)
    printf("Create the artificial ray pointer\n");
  for (j = 0; j < cone->d; j++)
    zerovector[j] = 0.0;
  cone->LastRay=cone->ArtificialRay;
  StoreRay1(cone, zerovector, &feasible);  
    /* This stores a vector to the record pointed by cone->LastRay */
  cone->ArtificialRay->Next = NULL;
}

void ConditionalAddEdge(dd_ConePtr cone, 
    dd_RayPtr Ray1, dd_RayPtr Ray2, dd_RayPtr ValidFirstRay)
{
  long it,it_row,fii1,fii2,fmin,fmax;
  boolean adjacent,lastchance;
  dd_RayPtr TempRay,Rmin,Rmax;
  dd_Adjacency *NewEdge;
  boolean localdebug=FALSE;
  dd_rowset ZSmin, ZSmax;
  static dd_rowset Face, Face1;
  static dd_rowrange last_m=0;
  
  if (last_m!=cone->m) {
    if (last_m>0){
      set_free(Face); set_free(Face1);
    }
    set_initialize(&Face, cone->m);
    set_initialize(&Face1, cone->m);
    last_m=cone->m;
  }
  
  fii1=Ray1->FirstInfeasIndex;
  fii2=Ray2->FirstInfeasIndex;
  if (fii1<fii2){
    fmin=fii1; fmax=fii2;
    Rmin=Ray1;
    Rmax=Ray2;
  }
  else{
    fmin=fii2; fmax=fii1;
    Rmin=Ray2;
    Rmax=Ray1;
  }
  ZSmin = Rmin->ZeroSet;
  ZSmax = Rmax->ZeroSet;
  if (localdebug) {
    printf("ConditionalAddEdge: FMIN = %ld (row%ld)   FMAX=%ld\n",
      fmin, cone->OrderVector[fmin], fmax);
  }
  if (fmin==fmax){
    if (localdebug) printf("ConditionalAddEdge: equal FII value-> No edge added\n");
  }
  else if (set_member(cone->OrderVector[fmin],ZSmax)){
    if (localdebug) printf("ConditionalAddEdge: No strong separation -> No edge added\n");
  }
  else {  /* the pair will be separated at the iteration fmin */
    lastchance=TRUE;
    /* flag to check it will be the last chance to store the edge candidate */
    set_int(Face1, ZSmax, ZSmin);
    (cone->count_int)++;
    if (localdebug){
      printf("Face: ");
      for (it=1; it<=cone->m; it++) {
        it_row=cone->OrderVector[it];
        if (set_member(it_row, Face1)) printf("%ld ",it_row);
      }
      printf("\n");
    }
    for (it=cone->Iteration+1; it < fmin && lastchance; it++){
      it_row=cone->OrderVector[it];
      if (cone->parent->EqualityIndex[it_row]>=0 && set_member(it_row, Face1)){
        lastchance=FALSE;
        (cone->count_int_bad)++;
        if (localdebug){
          printf("There will be another chance iteration %ld (row %ld) to store the pair\n", it, it_row);
        }
      }
    }
    if (lastchance){
      adjacent = TRUE;
      (cone->count_int_good)++;
      /* adjacent checking */
      set_int(Face, Face1, cone->AddedHalfspaces);
      if (localdebug){
        printf("Check adjacency\n");
        printf("AddedHalfspaces: "); set_write(cone->AddedHalfspaces);
        printf("Face: ");
        for (it=1; it<=cone->m; it++) {
          it_row=cone->OrderVector[it];
          if (set_member(it_row, Face)) printf("%ld ",it_row);
        }
        printf("\n");
      }
      if (set_card(Face)< cone->d - 2) {
        adjacent = FALSE;
      }
      else if (cone->parent->NondegAssumed) {
    	adjacent = TRUE;
      }
      else{
        TempRay = ValidFirstRay;  /* the first ray for adjacency checking */
        while (TempRay != NULL && adjacent) {
          if (TempRay != Ray1 && TempRay != Ray2) {
            set_int(Face1, TempRay->ZeroSet, cone->AddedHalfspaces);
            if (set_subset(Face, Face1)) {
              if (localdebug) set_write(Face1);
              adjacent = FALSE;
            }
          }
          TempRay = TempRay->Next;
        }
      }
      if (adjacent){
        if (localdebug) printf("The pair is adjacent and the pair must be stored for iteration %ld (row%ld)\n",
          fmin, cone->OrderVector[fmin]);
        NewEdge=(dd_AdjacencyPtr) malloc(sizeof *NewEdge);
        NewEdge->Ray1=Rmax;  /* save the one remains in iteration fmin in the first */
        NewEdge->Ray2=Rmin;  /* save the one deleted in iteration fmin in the second */
        NewEdge->Next=NULL;
        (cone->EdgeCount)++; 
        (cone->TotalEdgeCount)++;
        if (cone->Edges[fmin]==NULL){
          cone->Edges[fmin]=NewEdge;
          if (localdebug) printf("Create a new edge list of %ld\n", fmin);
        }else{
          NewEdge->Next=cone->Edges[fmin];
          cone->Edges[fmin]=NewEdge;
        }
      }
    }
  }
}

void CreateInitialEdges(dd_ConePtr cone)
{
  dd_RayPtr Ptr1, Ptr2;
  dd_rowrange fii1,fii2;
  long count=0;
  boolean adj,localdebug=FALSE;

  cone->Iteration=cone->d;  /* CHECK */
  if (cone->FirstRay ==NULL || cone->LastRay==NULL){
    printf("Error found: CreateInitialEdges called with NULL pointer(s)\n");
    goto _L99;
  }
  Ptr1=cone->FirstRay;
  while(Ptr1!=cone->LastRay && Ptr1!=NULL){
    fii1=Ptr1->FirstInfeasIndex;
    Ptr2=Ptr1->Next;
    while(Ptr2!=NULL){
      fii2=Ptr2->FirstInfeasIndex;
      count++;
      if (localdebug) printf("CreateInitialEdges: edge %ld \n",count);
      CheckAdjacency(cone, &Ptr1, &Ptr2, &adj);
      if (fii1!=fii2 && adj) 
        ConditionalAddEdge(cone, Ptr1, Ptr2, cone->FirstRay);
      Ptr2=Ptr2->Next;
    }
    Ptr1=Ptr1->Next;
  }
_L99:;  
}


void UpdateEdges(dd_ConePtr cone, dd_RayPtr RRbegin, dd_RayPtr RRend)
/* This procedure must be called after the ray list is sorted
   by EvaluateARay2 so that FirstInfeasIndex's are monotonically
   increasing.
*/
{
  dd_RayPtr Ptr1, Ptr2begin, Ptr2;
  dd_rowrange fii1;
  boolean ptr2found,quit,localdebug=FALSE;
  long count=0,pos1, pos2;
  float workleft,prevworkleft=110.0,totalpairs;

  totalpairs=(cone->ZeroRayCount-1.0)*(cone->ZeroRayCount-2.0)+1.0;
  Ptr2begin = NULL; 
  if (RRbegin ==NULL || RRend==NULL){
    if (1) printf("Warning: UpdateEdges called with NULL pointer(s)\n");
    goto _L99;
  }
  Ptr1=RRbegin;
  pos1=1;
  do{
    ptr2found=FALSE;
    quit=FALSE;
    fii1=Ptr1->FirstInfeasIndex;
    pos2=2;
    for (Ptr2=Ptr1->Next; !ptr2found && !quit; Ptr2=Ptr2->Next,pos2++){
      if  (Ptr2->FirstInfeasIndex > fii1){
        Ptr2begin=Ptr2;
        ptr2found=TRUE;
      }
      else if (Ptr2==RRend) quit=TRUE;
    }
    if (ptr2found){
      quit=FALSE;
      for (Ptr2=Ptr2begin; !quit ; Ptr2=Ptr2->Next){
        count++;
        if (localdebug) printf("UpdateEdges: edge %ld \n",count);
        ConditionalAddEdge(cone, Ptr1,Ptr2,RRbegin);
        if (Ptr2==RRend || Ptr2->Next==NULL) quit=TRUE;
      }
    }
    Ptr1=Ptr1->Next;
    pos1++;
    workleft = 100.0 * (cone->ZeroRayCount-pos1) * (cone->ZeroRayCount - pos1-1.0) / totalpairs;
    if (cone->ZeroRayCount>=500 && debug && pos1%10 ==0 && prevworkleft-workleft>=10 ) {
      printf("*Work of iteration %5ld(/%ld): %4ld/%4ld => %4.1f%% left\n",
	     cone->Iteration, cone->m, pos1, cone->ZeroRayCount, workleft);
      prevworkleft=workleft;
    }    
  }while(Ptr1!=RRend && Ptr1!=NULL);
_L99:;  
}

void FreeDDMemory0(dd_ConePtr *cone)
{
  dd_RayPtr Ptr, PrevPtr;
  long count;
  dd_rowrange i;
  boolean localdebug=FALSE;
  
  /* THIS SHOULD BE REWRITTEN carefully */
  PrevPtr=(*cone)->ArtificialRay;
  if (PrevPtr!=NULL){
    count=0;
    for (Ptr=(*cone)->ArtificialRay->Next; Ptr!=NULL; Ptr=Ptr->Next){
      free(PrevPtr->Ray);
      free(PrevPtr->ZeroSet);
      free(PrevPtr);
      count++;
      PrevPtr=Ptr;
    };
    (*cone)->LastRay=NULL;
    (*cone)->FirstRay=NULL;
    (*cone)->ArtificialRay=NULL;
    if (localdebug) printf("%ld ray storage spaces freed\n",count);
  }
  

  set_free((*cone)->GroundSet); 
  set_free((*cone)->EqualitySet); 
  set_free((*cone)->NonequalitySet); 
  set_free((*cone)->AddedHalfspaces); 
  set_free((*cone)->WeaklyAddedHalfspaces); 
  set_free((*cone)->InitialHalfspaces);
  free((*cone)->InitialRayIndex);
  free((*cone)->OrderVector);
  free((*cone)->newcol);
  dd_FreeBmatrix((*cone)->d,&((*cone)->B));
  dd_FreeBmatrix((*cone)->d,&((*cone)->Bsave));

  for (i=0; i<(*cone)->m; i++)
    free((*cone)->A[i]);

  free((*cone)->A);

  free(*cone);
}

void dd_FreeDDMemory(dd_PolyhedraPtr poly)
{
  FreeDDMemory0(&(poly->child));
}

void Normalize(dd_colrange d_size, double *V)
{
  long j;
  double min, temp;

  min = 1.0e+20;
  for (j = 0; j < d_size; j++) {
    temp = fabs(V[j]);
/*    if (temp > zero && temp < min) */
    if (dd_Positive(temp) && temp < min)
      min = temp;
  }
  for (j = 0; j < d_size; j++)
    V[j] /= min;
}


void ZeroIndexSet(dd_rowrange m_size, dd_colrange d_size, dd_Amatrix A, double *x, dd_rowset ZS)
{
  dd_rowrange i;
  double temp;

  set_emptyset(ZS);
  for (i = 1; i <= m_size; i++) {
    temp = AValue(d_size, A, x, i);
/*    if (fabs(temp) < zero)*/
    if (dd_Zero(temp))
      set_addelem(ZS, i);
  }
}

void CopyBmatrix(dd_colrange d_size, dd_Bmatrix T, dd_Bmatrix TCOPY)
{
  dd_rowrange i;
  dd_colrange j;

  for (i=0; i < d_size; i++) {
    for (j=0; j < d_size; j++) {
      TCOPY[i][j] = T[i][j];
    }
  }
}

double TableauEntry(dd_rowrange m_size, dd_colrange d_size, dd_Amatrix X, dd_Bmatrix T,
				dd_rowrange r, dd_colrange s)
/* Compute the (r,s) entry of X.T   */
{
  dd_colrange j;
  double temp;
  
  temp=0;
  for (j=0; j< d_size; j++) {
    temp = temp + X[r-1][j] * T[j][s-1];
  }
  return temp;
}


void SelectPivot(dd_ConePtr cone, dd_rowrange rowmax, dd_rowset NopivotRow,
            dd_colset NopivotCol, dd_rowrange *r, dd_colrange *s,
            boolean *selected)
/* Select a position (*r,*s) in the cone matrix A.B such that (A.B)[*r][*s] is nonzero
   The choice is feasible, i.e., not on NopivotRow and NopivotCol, and
   best with respect to the specified roworder 
 */
{
  boolean stop;
  dd_rowrange i,rtemp;
  dd_rowset rowexcluded;
  double Xtemp;
  boolean localdebug=FALSE;

  if (debug) localdebug=TRUE;
  stop = FALSE;
  set_initialize(&rowexcluded, cone->m);
  set_copy(rowexcluded,NopivotRow);
  if (localdebug) {
    switch (cone->HalfspaceOrder) {

    case MinIndex:
      fprintf(stdout, "*SelectPivot: MinIndex\n");
      break;

    case MaxIndex:
      fprintf(stdout, "*SelectPivot: MaxIndex\n");
      break;

    case MinCutoff:
      fprintf(stdout, "*SelectPivot: MinCutoff\n");
      break;

    case MaxCutoff:
      fprintf(stdout, "*SelectPivot: MaxCutoff\n");
    break;

    case MixCutoff:
      fprintf(stdout, "*SelectPivot: MixCutoff\n");
      break;

    case LexMin:
      fprintf(stdout, "*SelectPivot: LexMin\n");
      break;

    case LexMax:
      fprintf(stdout, "*SelectPivot: LexMax\n");
      break;

    case RandomRow:
      fprintf(stdout, "*SelectPivot: Random\n");
      break;

    case LineShelling:
      fprintf(stdout, "*SelectPivot: LineShelling\n");
      break;
    }
    printf("select pivot: rowexcluded=");
    set_fwrite(stdout,rowexcluded);
  }
  for (rtemp=rowmax+1;rtemp<=cone->m;rtemp++) {
    set_addelem(rowexcluded,rtemp);   /* cannot pivot on any row > rmax */
  }
  *selected = FALSE;
  do {
    i=1;rtemp=0;
    while (i<=cone->m && rtemp==0) {  /* EqualitySet vars have highest priorities */
      if (set_member(i,cone->EqualitySet) && !set_member(i,rowexcluded)){
        if (localdebug) printf("marked set %ld chosen as a candidate\n",i);
        rtemp=i;
      }
      i++;
    }
    if (rtemp==0) SelectNextHalfspace(cone, rowexcluded, &rtemp);
    if (rtemp>=1) {
      *r=rtemp;
      *s=1;
      while (*s <= cone->d && !*selected) {
        Xtemp=TableauEntry(cone->m,cone->d, cone->A, cone->B,*r,*s);
/*        if (!set_member(*s,NopivotCol) && fabs(Xtemp) > zero) {*/
        if (!set_member(*s,NopivotCol) && dd_Nonzero(Xtemp)) {
          *selected = TRUE;
          stop = TRUE;
        } else {
          (*s)++;
        }
      }
      if (!*selected) {
        set_addelem(rowexcluded, rtemp);
      }
    }
    else {
      *r = 0;
      *s = 0;
      stop = TRUE;
    }
  } while (!stop);
  set_free(rowexcluded);
}

void GaussianColumnPivot(dd_rowrange m_size, dd_colrange d_size, 
    dd_Amatrix X, dd_Bmatrix T, dd_rowrange r, dd_colrange s)
/* Update the Transformation matrix T with the pivot operation on (r,s) 
   This procedure performs a implicit pivot operation on the matrix X by
   updating the dual basis inverse  T.
 */
{
  long j, j1;
  double Xtemp0, Xtemp;
  static dd_Arow Rtemp;
  static dd_colrange last_d=0;

  if (last_d!=d_size){
    if (last_d>0) free(Rtemp);
    Rtemp=(double*)calloc(d_size,sizeof(double));
    last_d=d_size;
  }

  for (j=1; j<=d_size; j++) Rtemp[j-1]=TableauEntry(m_size, d_size, X, T, r,j);
  Xtemp0 = Rtemp[s-1];
  for (j = 1; j <= d_size; j++) {
    if (j != s) {
      Xtemp = Rtemp[j-1];
      for (j1 = 1; j1 <= d_size; j1++)
        T[j1-1][j-1] -= T[j1-1][s - 1] * Xtemp / Xtemp0;
    }
  }
  for (j = 1; j <= d_size; j++)
    T[j-1][s - 1] /= Xtemp0;
}


void CopyArow(double *acopy, double *a, dd_colrange d)
{
  dd_colrange j;

  for (j = 0; j < d; j++) {
    acopy[j]=a[j];
  }
}

void CopyAmatrix(double **Acopy, double **A, dd_rowrange m, dd_colrange d)
{
  dd_rowrange i;

  for (i = 0; i< m; i++) {
    CopyArow(Acopy[i],A[i],d);
  }
}


void dd_InitializeArow(dd_colrange d,dd_Arow *a)
{
  if (d>0) *a=(double*) calloc(d,sizeof(double));
}

void dd_InitializeAmatrix(dd_rowrange m,dd_colrange d,dd_Amatrix *A)
{
  dd_rowrange i;

  (*A)=(double**) calloc(m,sizeof(double*));
  for (i = 0; i < m; i++) {
    dd_InitializeArow(d,&((*A)[i]));
  }
}

void dd_FreeAmatrix(dd_rowrange m,dd_Amatrix *A)
{
  dd_rowrange i;

  if (*A!=NULL) {
    for (i = 0; i < m; i++) {
      free((*A)[i]);
    }
    free(*A);
  }
}

void dd_InitializeBmatrix(dd_colrange d,dd_Bmatrix *B)
{
  dd_colrange j;

  (*B)=(double**) calloc(d,sizeof(double*));
  for (j = 0; j < d; j++) {
    (*B)[j]=(double*) calloc(d,sizeof(double));
  }
}

void dd_FreeBmatrix(dd_colrange d,dd_Bmatrix *B)
{
  dd_colrange j;

  if (*B!=NULL) {
    for (j = 0; j < d; j++) {
      free((*B)[j]);
    }
    free(*B);
  }
}

dd_SetFamilyPtr dd_CreateSetFamily(dd_bigrange fsize, dd_bigrange ssize)
{
  dd_SetFamilyPtr F;
  dd_bigrange i;

  F=(dd_SetFamilyPtr) malloc (sizeof(dd_SetFamily));
  F->set=(set_type*) calloc(fsize,sizeof(set_type));
  for (i=0; i<fsize; i++) {
    set_initialize(&(F->set[i]), ssize);
  }
  F->famsize=fsize;
  F->setsize=ssize;
  return F;
}


void dd_FreeSetFamily(dd_SetFamilyPtr *F)
{
  dd_bigrange i;

  if (*F!=NULL){
    for (i=0; i<(*F)->famsize; i++) {
      set_free((*F)->set[i]);
    }
    free((*F)->set);
    free(*F);
  }
}

dd_MatrixPtr dd_CreateMatrix(dd_rowrange m_size,dd_colrange d_size)
{
  dd_MatrixPtr M;
  dd_rowrange m1;

  if (m_size<=0) m1=1; else m1=m_size;
  M=(dd_MatrixPtr) malloc (sizeof(dd_Matrix));
  dd_InitializeAmatrix(m_size,d_size,&(M->matrix));
  M->rowsize=m_size;
  set_initialize(&(M->linset), m1);
  M->colsize=d_size;
  return M;
}

void dd_FreeMatrix(dd_MatrixPtr *M)
{
  if (*M!=NULL) {
    dd_FreeAmatrix((*M)->rowsize,&((*M)->matrix));
    set_free((*M)->linset);
    free(*M);
  }
}

void dd_SetToIdentity(dd_colrange d_size, dd_Bmatrix T)
{
  dd_colrange j1, j2;

  for (j1 = 1; j1 <= d_size; j1++) {
    for (j2 = 1; j2 <= d_size; j2++) {
      if (j1 == j2)
        T[j1 - 1][j2 - 1] = 1.0;
      else
        T[j1 - 1][j2 - 1] = 0.0;
    }
  }
}

void ColumnReduce(dd_ConePtr cone)
{
  dd_colrange j,j1=0;
  dd_rowrange i;
  boolean localdebug=FALSE;

  for (j=1;j<=cone->d;j++) {
    if (cone->InitialRayIndex[j]>0){
      j1=j1+1;
      if (j1<j) {
        for (i=1; i<=cone->m; i++) cone->A[i-1][j1-1]=cone->A[i-1][j-1];
        cone->newcol[j]=j1;
        if (localdebug){
          printf("shifting the column %ld to column %ld\n", j, j1);
        }
          /* shifting forward */
      }
    } else {
      cone->newcol[j]=0;
      if (localdebug) {
        printf("a generator (or an equation) of the linearity space: ");
        for (i=1; i<=cone->d; i++) dd_WriteReal(stdout, cone->B[i-1][j-1]);
        printf("\n");
      }
    }
  }
  cone->d=j1;  /* update the dimension. cone->d_orig remembers the old. */
  CopyBmatrix(cone->d_orig, cone->B, cone->Bsave);  
    /* save the dual basis inverse as Bsave.  This matrix contains the linearity space generators. */
  cone->ColReduced=TRUE;
}

void FindBasis(dd_ConePtr cone, long *rank)
{
  boolean stop, chosen, localdebug=FALSE;
  dd_rowset NopivotRow;
  dd_colset ColSelected;
  dd_rowrange r;
  dd_colrange j,s;

  *rank = 0;
  stop = FALSE;
  for (j=0;j<=cone->d;j++) cone->InitialRayIndex[j]=0;
  set_emptyset(cone->InitialHalfspaces);
  set_initialize(&ColSelected, cone->d);
  set_initialize(&NopivotRow, cone->m);
  set_copy(NopivotRow,cone->NonequalitySet);
  dd_SetToIdentity(cone->d, cone->B);
  if (localdebug){
    printf("*Initial set of rows:");
  }
  do {   /* Find a set of rows for a basis */
      SelectPivot(cone, cone->m, NopivotRow, ColSelected, &r, &s, &chosen);
      if (debug && chosen) 
        printf("Procedure FindBasis: pivot on (r,s) =(%ld, %ld).\n", r, s);
      if (chosen) {
        set_addelem(cone->InitialHalfspaces, r);
        set_addelem(NopivotRow, r);
        set_addelem(ColSelected, s);
        cone->InitialRayIndex[s]=r;    /* cone->InitialRayIndex[s] stores the corr. row index */
        (*rank)++;
        GaussianColumnPivot(cone->m, cone->d, cone->A, cone->B, r, s);
	if (localdebug){
	  printf(" %ld",  r);
	}
      } else {
        stop=TRUE;
      }
      if (*rank==cone->d) stop = TRUE;
  } while (!stop);
  set_free(ColSelected);
}


void FindInitialRays(dd_ConePtr cone, boolean *found)
{
  dd_rowset CandidateRows;
  dd_rowrange i;
  long rank;
  dd_HalfspaceOrderType roworder_save=LexMin;

  *found = FALSE;
  set_initialize(&CandidateRows, cone->m);
  if (cone->parent->InitBasisAtBottom==TRUE) {
    roworder_save=cone->HalfspaceOrder;
    cone->HalfspaceOrder=MaxIndex;
    cone->PreOrderedRun=FALSE;
  }
  else cone->PreOrderedRun=TRUE;
  for (i = 1; i <= cone->m; i++)
    if (!set_member(i,cone->NonequalitySet)) set_addelem(CandidateRows, i);
    /*all rows not in NonequalitySet are candidates for initial cone*/
  FindBasis(cone, &rank);
  if (debug) dd_WriteBmatrix(stdout, cone->d, cone->B);
  if (debug) printf("FindInitialRays: rank of Amatrix = %ld\n", rank);
  cone->LinearityDim=cone->d - rank;
  if (debug) printf("Linearity Dimension = %ld\n", cone->LinearityDim);
  if (cone->LinearityDim > 0) {
     cone->Error = LowColumnRank;
     ColumnReduce(cone);
     FindBasis(cone, &rank);
  }
  if (!set_subset(cone->EqualitySet,cone->InitialHalfspaces)) {
    if (debug) {
      printf("Equality set is dependent. Equality Set and an initial basis:\n");
      set_write(cone->EqualitySet);
      set_write(cone->InitialHalfspaces);
    };
  }
  *found = TRUE;
  set_free(CandidateRows);
  if (cone->parent->InitBasisAtBottom==TRUE) {
    cone->HalfspaceOrder=roworder_save;
  }
  if (cone->HalfspaceOrder==MaxCutoff||
      cone->HalfspaceOrder==MinCutoff||
      cone->HalfspaceOrder==MixCutoff){
    cone->PreOrderedRun=FALSE;
  } else cone->PreOrderedRun=TRUE;
}

void CheckEquality(dd_colrange d_size, dd_RayPtr*RP1, dd_RayPtr*RP2, boolean *equal)
{
  long j;

  if (debug)
    printf("Check equality of two rays\n");
  *equal = TRUE;
  j = 1;
  while (j <= d_size && *equal) {
/*    if (fabs((*RP1)->Ray[j - 1] - (*RP2)->Ray[j - 1]) > 2.0 * zero)*/
    if (dd_Positive((*RP1)->Ray[j - 1] - (*RP2)->Ray[j - 1]))
      *equal = FALSE;
    j++;
  }
  if (*equal)
    printf("Equal records found !!!!\n");
}

void CreateNewRay(dd_ConePtr cone, 
    dd_RayPtr Ptr1, dd_RayPtr Ptr2, dd_rowrange ii)
{
  /*Create a new ray by taking a linear combination of two rays*/
  dd_colrange j;
  double v1, v2;
  static dd_Arow NewRay;
  static dd_colrange last_d=0;

  if (last_d!=cone->d){
    if (last_d>0) free(NewRay);
    NewRay=(double*)calloc(cone->d,sizeof(double));
    last_d=cone->d;
  }

  v1 = fabs(AValue(cone->d, cone->A, Ptr1->Ray, ii));
  v2 = fabs(AValue(cone->d, cone->A, Ptr2->Ray, ii));
  for (j = 0; j < cone->d; j++)
    NewRay[j] = Ptr1->Ray[j] * v2 + Ptr2->Ray[j] * v1;
  Normalize(cone->d, NewRay);
  AddRay(cone, NewRay);
}

void EvaluateARay1(dd_rowrange i, dd_ConePtr cone)
/* Evaluate the ith component of the vector  A x RD.Ray 
    and rearrange the linked list so that
    the infeasible rays with respect to  i  will be
    placed consecutively from First 
 */
{
  dd_colrange j;
  double temp;
  dd_RayPtr Ptr, PrevPtr, TempPtr;

  Ptr = cone->FirstRay;
  PrevPtr = cone->ArtificialRay;
  if (PrevPtr->Next != Ptr) {
    printf("Error.  Artificial Ray does not point to FirstRay!!!\n");
  }
  while (Ptr != NULL) {
    temp = 0.0;
    for (j = 0; j < cone->d; j++)
      temp += cone->A[i - 1][j] * Ptr->Ray[j];
    Ptr->ARay = temp;
/*    if ( temp <= -zero && Ptr != cone->FirstRay) {*/
    if ( dd_Negative(temp) && Ptr != cone->FirstRay) {
      /* printf("Moving an infeasible record w.r.t. %ld to FirstRay\n",i); */
      if (Ptr==cone->LastRay) cone->LastRay=PrevPtr;
      TempPtr=Ptr;
      Ptr = Ptr->Next;
      PrevPtr->Next = Ptr;
      cone->ArtificialRay->Next = TempPtr;
      TempPtr->Next = cone->FirstRay;
      cone->FirstRay = TempPtr;
    }
    else {
      PrevPtr = Ptr;
      Ptr = Ptr->Next;
    }
  }
}

void EvaluateARay2(dd_rowrange i, dd_ConePtr cone)
/* Evaluate the ith component of the vector  A x RD.Ray 
   and rearrange the linked list so that
   the infeasible rays with respect to  i  will be
   placed consecutively from First. Also for all feasible rays,
   "positive" rays and "zero" rays will be placed consecutively.
 */
{
  dd_colrange j;
  double temp;
  dd_RayPtr Ptr, NextPtr;
  boolean zerofound=FALSE,negfound=FALSE,posfound=FALSE;

  cone->PosHead=NULL;cone->ZeroHead=NULL;cone->NegHead=NULL;
  cone->PosLast=NULL;cone->ZeroLast=NULL;cone->NegLast=NULL;
  Ptr = cone->FirstRay;
  while (Ptr != NULL) {
    NextPtr=Ptr->Next;  /* remember the Next record */
    Ptr->Next=NULL;     /* then clear the Next pointer */
    temp = 0.0;
    for (j = 0; j < cone->d; j++)
      temp += cone->A[i - 1][j] * Ptr->Ray[j];
    Ptr->ARay = temp;
/*    if ( temp < -zero) {*/
    if ( dd_Negative(temp)) {
      if (!negfound){
        negfound=TRUE;
        cone->NegHead=Ptr;
        cone->NegLast=Ptr;
      }
      else{
        Ptr->Next=cone->NegHead;
        cone->NegHead=Ptr;
      }
    }
/*    else if (temp > zero){*/
    else if (dd_Positive(temp)){
      if (!posfound){
        posfound=TRUE;
        cone->PosHead=Ptr;
        cone->PosLast=Ptr;
      }
      else{  
        Ptr->Next=cone->PosHead;
        cone->PosHead=Ptr;
       }
    }
    else {
      if (!zerofound){
        zerofound=TRUE;
        cone->ZeroHead=Ptr;
        cone->ZeroLast=Ptr;
      }
      else{
        Ptr->Next=cone->ZeroHead;
        cone->ZeroHead=Ptr;
      }
    }
    Ptr=NextPtr;
  }
  /* joining three neg, pos and zero lists */
  if (negfound){                 /* -list nonempty */
    cone->FirstRay=cone->NegHead;
    if (posfound){               /* -list & +list nonempty */
      cone->NegLast->Next=cone->PosHead;
      if (zerofound){            /* -list, +list, 0list all nonempty */
        cone->PosLast->Next=cone->ZeroHead;
        cone->LastRay=cone->ZeroLast;
      } 
      else{                      /* -list, +list nonempty but  0list empty */
        cone->LastRay=cone->PosLast;      
      }
    }
    else{                        /* -list nonempty & +list empty */
      if (zerofound){            /* -list,0list nonempty & +list empty */
        cone->NegLast->Next=cone->ZeroHead;
        cone->LastRay=cone->ZeroLast;
      } 
      else {                      /* -list nonempty & +list,0list empty */
        cone->LastRay=cone->NegLast;
      }
    }
  }
  else if (posfound){            /* -list empty & +list nonempty */
    cone->FirstRay=cone->PosHead;
    if (zerofound){              /* -list empty & +list,0list nonempty */
      cone->PosLast->Next=cone->ZeroHead;
      cone->LastRay=cone->ZeroLast;
    } 
    else{                        /* -list,0list empty & +list nonempty */
      cone->LastRay=cone->PosLast;
    }
  }
  else{                          /* -list,+list empty & 0list nonempty */
    cone->FirstRay=cone->ZeroHead;
    cone->LastRay=cone->ZeroLast;
  }
  cone->ArtificialRay->Next=cone->FirstRay;
  cone->LastRay->Next=NULL;
}

void DeleteNegativeRays(dd_ConePtr cone)
/* Eliminate the infeasible rays with respect to  i  which
   are supposed to be consecutive from the head of the dd_Ray list,
   and sort the zero list assumed to be consecutive at the
   end of the list.
 */
{
  dd_rowrange fii,fiitest;
  double temp;
  dd_RayPtr Ptr, PrevPtr,NextPtr,ZeroPtr1,ZeroPtr0;
  boolean found, completed, zerofound=FALSE,negfound=FALSE,posfound=FALSE;
  boolean localdebug=FALSE;
  
  cone->PosHead=NULL;cone->ZeroHead=NULL;cone->NegHead=NULL;
  cone->PosLast=NULL;cone->ZeroLast=NULL;cone->NegLast=NULL;

  /* Delete the infeasible rays  */
  PrevPtr= cone->ArtificialRay;
  Ptr = cone->FirstRay;
  if (PrevPtr->Next != Ptr) 
    printf("Error at DeleteNegativeRays: ArtificialRay does not point the FirstRay.\n");
  completed=FALSE;
  while (Ptr != NULL && !completed){
/*    if ( (Ptr->ARay) < -zero ){ */
    if ( dd_Negative(Ptr->ARay)){
      Eliminate(cone, &PrevPtr);
      Ptr=PrevPtr->Next;
    }
    else{
      completed=TRUE;
    }
  }
  
  /* Sort the zero rays */
  Ptr = cone->FirstRay;
  cone->ZeroRayCount=0;
  while (Ptr != NULL) {
    NextPtr=Ptr->Next;  /* remember the Next record */
    temp = Ptr->ARay;
    if (localdebug) printf("Ptr->ARay : %5.3f \n", temp);
/*    if ( temp < -zero) {*/
    if ( dd_Negative(temp)) {
      if (!negfound){
        printf("Error: An infeasible ray found after their removal\n");
        negfound=TRUE;
      }
    }
/*    else if (temp > zero){*/
    else if (dd_Positive(temp)){
      if (!posfound){
        posfound=TRUE;
        cone->PosHead=Ptr;
        cone->PosLast=Ptr;
      }
      else{  
        cone->PosLast=Ptr;
       }
    }
    else {
      (cone->ZeroRayCount)++;
      if (!zerofound){
        zerofound=TRUE;
        cone->ZeroHead=Ptr;
        cone->ZeroLast=Ptr;
        cone->ZeroLast->Next=NULL;
      }
      else{/* Find a right position to store the record sorted w.r.t. FirstInfeasIndex */
        fii=Ptr->FirstInfeasIndex; 
        found=FALSE;
        ZeroPtr1=NULL;
        for (ZeroPtr0=cone->ZeroHead; !found && ZeroPtr0!=NULL ; ZeroPtr0=ZeroPtr0->Next){
          fiitest=ZeroPtr0->FirstInfeasIndex;
          if (fiitest >= fii){
            found=TRUE;
          }
          else ZeroPtr1=ZeroPtr0;
        }
        /* printf("insert position found \n %d  index %ld\n",found, fiitest); */
        if (!found){           /* the new record must be stored at the end of list */
          cone->ZeroLast->Next=Ptr;
          cone->ZeroLast=Ptr;
          cone->ZeroLast->Next=NULL;
        }
        else{
          if (ZeroPtr1==NULL){ /* store the new one at the head, and update the head ptr */
            /* printf("Insert at the head\n"); */
            Ptr->Next=cone->ZeroHead;
            cone->ZeroHead=Ptr;
          }
          else{                /* store the new one inbetween ZeroPtr1 and 0 */
            /* printf("Insert inbetween\n");  */
            Ptr->Next=ZeroPtr1->Next;
            ZeroPtr1->Next=Ptr;
          }
        }
        /*
        Ptr->Next=cone->ZeroHead;
        cone->ZeroHead=Ptr;
        */
      }
    }
    Ptr=NextPtr;
  }
  /* joining the pos and zero lists */
  if (posfound){            /* -list empty & +list nonempty */
    cone->FirstRay=cone->PosHead;
    if (zerofound){              /* +list,0list nonempty */
      cone->PosLast->Next=cone->ZeroHead;
      cone->LastRay=cone->ZeroLast;
    } 
    else{                        /* 0list empty & +list nonempty */
      cone->LastRay=cone->PosLast;
    }
  }
  else{                          /* +list empty & 0list nonempty */
    cone->FirstRay=cone->ZeroHead;
    cone->LastRay=cone->ZeroLast;
  }
  cone->ArtificialRay->Next=cone->FirstRay;
  cone->LastRay->Next=NULL;
}

void FeasibilityIndices(long *fnum, long *infnum, dd_rowrange i, dd_ConePtr cone)
{
  /*Evaluate the number of feasible rays and infeasible rays*/
  /*  w.r.t the hyperplane  i*/
  dd_colrange j;
  double temp;
  dd_RayPtr Ptr;

  *fnum = 0;
  *infnum = 0;
  Ptr = cone->FirstRay;
  while (Ptr != NULL) {
    temp = 0.0;
    for (j = 0; j < cone->d; j++)
      temp += cone->A[i - 1][j] * Ptr->Ray[j];
    if (temp >= 0)
      (*fnum)++;
    else
      (*infnum)++;
    Ptr = Ptr->Next;
  }
}

boolean LexSmaller(double *v1, double *v2, long dmax)
{ /* dmax is the size of vectors v1,v2 */
  boolean determined, smaller;
  dd_colrange j;

  smaller = FALSE;
  determined = FALSE;
  j = 1;
  do {
/*    if (fabs(v1[j - 1]-v2[j - 1])>zero) {*/
    if (dd_Positive(v1[j - 1]-v2[j - 1])) {
      if (v1[j - 1] < v2[j - 1]) {
	    smaller = TRUE;
	  }
      determined = TRUE;
    } else
      j++;
  } while (!(determined) && (j <= dmax));
  return smaller;
}


boolean LexLarger(double *v1, double *v2, long dmax)
{
  dd_colrange j;
  static dd_Arow u1,u2;
  static dd_colrange last_d=0;

  if (last_d < dmax){
    if (last_d>0) {free(u1); free(u2);}
    u1=(double*)calloc(dmax,sizeof(double));
    u2=(double*)calloc(dmax,sizeof(double));
    last_d=dmax;
  }

  for (j = 1; j <= dmax; j++) {
    u1[j-1] = -v1[j-1];
    u2[j-1] = -v2[j-1];
  }
  return (LexSmaller(u1, u2, dmax));
}

void Copydd_Arow(double *vcopy, double *v, long dmax)
{
 dd_colrange j;

  for (j = 1; j <= dmax; j++) {
    vcopy[j-1] = v[j-1];
  }
}

void AddNewHalfspace1(dd_ConePtr cone, dd_rowrange hnew)
/* This procedure 1 must be used with PreorderedRun=FALSE 
   This procedure is the most elementary implementation of
   DD and can be used with any type of ordering, including
   dynamic ordering of rows, e.g. MaxCutoff, MinCutoff.
   The memory requirement is minimum because it does not
   store any adjacency among the rays.
*/
{
  dd_RayPtr RayPtr0,RayPtr1,RayPtr2,RayPtr2s,RayPtr3;
  long pos1, pos2;
  double prevprogress, progress, value1, value2;
  boolean adj, equal, completed;

  EvaluateARay1(hnew, cone);  
   /*Check feasibility of rays w.r.t. hnew 
     and put all infeasible ones consecutively */

  RayPtr0 = cone->ArtificialRay;   /*Pointer pointing RayPrt1*/
  RayPtr1 = cone->FirstRay;        /*1st hnew-infeasible ray to scan and compare with feasible rays*/
  value1 = cone->FirstRay->ARay;
/*  if (value1 > -zero) {*/
  if (dd_Nonnegative(value1)) {
    if (cone->RayCount==cone->WeaklyFeasibleRayCount) cone->CompStatus=AllFound;
    goto _L99;        /* Sicne there is no hnew-infeasible ray and nothing to do */
  }
  else {
    RayPtr2s = RayPtr1->Next;/* RayPtr2s must point the first feasible ray */
    pos2=1;
/*    while (RayPtr2s!=NULL && RayPtr2s->ARay <= -zero) {*/
    while (RayPtr2s!=NULL && dd_Negative(RayPtr2s->ARay)) {
      RayPtr2s = RayPtr2s->Next;
      pos2++;
    }
  }
  if (RayPtr2s==NULL) {
    cone->FirstRay=NULL;
    cone->ArtificialRay->Next=cone->FirstRay;
    cone->RayCount=0;
    cone->CompStatus=AllFound;
    goto _L99;   /* All rays are infeasible, and the computation must stop */
  }
  RayPtr2 = RayPtr2s;   /*2nd feasible ray to scan and compare with 1st*/
  RayPtr3 = cone->LastRay;    /*Last feasible for scanning*/
  prevprogress=-10.0;
  pos1 = 1;
  completed=FALSE;
  while ((RayPtr1 != RayPtr2s) && !completed) {
    value1 = RayPtr1->ARay;
    value2 = RayPtr2->ARay;
    CheckEquality(cone->d, &RayPtr1, &RayPtr2, &equal);
/*    if ((value1 >= zero && value2 <= -zero) || (value2 >= zero && value1 <= -zero)) {*/
    if ((dd_Positive(value1) && dd_Negative(value2)) || (dd_Negative(value1) && dd_Positive(value2))){
      CheckAdjacency(cone, &RayPtr1, &RayPtr2, &adj);
      if (adj) CreateNewRay(cone, RayPtr1, RayPtr2, hnew);
    }
    if (RayPtr2 != RayPtr3) {
      RayPtr2 = RayPtr2->Next;
      continue;
    }
/*    if (value1 <= -zero || equal) {*/
    if (dd_Negative(value1) || equal) {
      Eliminate(cone, &RayPtr0);
      RayPtr1 = RayPtr0->Next;
      RayPtr2 = RayPtr2s;
    } else {
      completed=TRUE;
    }
    pos1++;
    progress = 100.0 * ((double)pos1 / pos2) * (2.0 * pos2 - pos1) / pos2;
    if (progress-prevprogress>=10 && pos1%10==0 && debug) {
      printf("*Progress of iteration %5ld(/%ld):   %4ld/%4ld => %4.1f%% done\n",
	     cone->Iteration, cone->m, pos1, pos2, progress);
      prevprogress=progress;
    }
  }
  if (cone->RayCount==cone->WeaklyFeasibleRayCount) cone->CompStatus=AllFound;
  _L99:;
}

void AddNewHalfspace2(dd_ConePtr cone, dd_rowrange hnew)
/* This procedure must be used under PreOrderedRun mode */
{
  dd_RayPtr RayPtr0,RayPtr1,RayPtr2;
  dd_Adjacency *EdgePtr, *EdgePtr0;
  long pos1;
  dd_rowrange fii1, fii2;
  boolean localdebug=FALSE;

  EvaluateARay2(hnew, cone);
   /* Check feasibility of rays w.r.t. hnew 
      and sort them. ( -rays, +rays, 0rays)*/

  if (cone->PosHead==NULL && cone->ZeroHead==NULL) {
    cone->FirstRay=NULL;
    cone->ArtificialRay->Next=cone->FirstRay;
    cone->RayCount=0;
    cone->CompStatus=AllFound;
    goto _L99;   /* All rays are infeasible, and the computation must stop */
  }
  
  if (localdebug){
    pos1=0;
    printf("(pos, FirstInfeasIndex, A Ray)=\n");
    for (RayPtr0=cone->FirstRay; RayPtr0!=NULL; RayPtr0=RayPtr0->Next){
      pos1++;
      printf("(%ld,%ld,",pos1,RayPtr0->FirstInfeasIndex);
      dd_WriteReal(stdout,RayPtr0->ARay); 
      printf(") ");
   }
    printf("\n");
  }
  
  if (cone->ZeroHead==NULL) cone->ZeroHead=cone->LastRay;

  EdgePtr=cone->Edges[cone->Iteration];
  while (EdgePtr!=NULL){
    RayPtr1=EdgePtr->Ray1;
    RayPtr2=EdgePtr->Ray2;
    fii1=RayPtr1->FirstInfeasIndex;   
    CreateNewRay(cone, RayPtr1, RayPtr2, hnew);
    fii2=cone->LastRay->FirstInfeasIndex;
    if (fii1 != fii2) 
      ConditionalAddEdge(cone,RayPtr1,cone->LastRay,cone->PosHead);
    EdgePtr0=EdgePtr;
    EdgePtr=EdgePtr->Next;
    free(EdgePtr0);
    (cone->EdgeCount)--;
  }
  cone->Edges[cone->Iteration]=NULL;
  
  DeleteNegativeRays(cone);
    
  set_addelem(cone->AddedHalfspaces, hnew);

  if (cone->Iteration<cone->m){
    if (cone->ZeroHead!=NULL && cone->ZeroHead!=cone->LastRay){
      if (cone->ZeroRayCount>200 && debug) printf("*New edges being scanned...\n");
      UpdateEdges(cone, cone->ZeroHead, cone->LastRay);
    }
  }

  if (cone->RayCount==cone->WeaklyFeasibleRayCount) cone->CompStatus=AllFound;
_L99:;
}


void SelectNextHalfspace0(dd_ConePtr cone, dd_rowset excluded, dd_rowrange *hnext)
{
  /*A natural way to choose the next hyperplane.  Simply the largest index*/
  long i;
  boolean determined;

  i = cone->m;
  determined = FALSE;
  do {
    if (set_member(i, excluded))
      i--;
    else
      determined = TRUE;
  } while (!determined && i>=1);
  if (determined) 
    *hnext = i;
  else
    *hnext = 0;
}

void SelectNextHalfspace1(dd_ConePtr cone, dd_rowset excluded, dd_rowrange *hnext)
{
  /*Natural way to choose the next hyperplane.  Simply the least index*/
  long i;
  boolean determined;

  i = 1;
  determined = FALSE;
  do {
    if (set_member(i, excluded))
      i++;
    else
      determined = TRUE;
  } while (!determined && i<=cone->m);
  if (determined) 
    *hnext = i;
  else 
    *hnext=0;
}

void SelectNextHalfspace2(dd_ConePtr cone, dd_rowset excluded, dd_rowrange *hnext)
{
  /*Choose the next hyperplane with maximum infeasibility*/
  long i, fea, inf, infmin, fi=0;   /*feasibility and infeasibility numbers*/

  infmin = cone->RayCount + 1;
  for (i = 1; i <= cone->m; i++) {
    if (!set_member(i, excluded)) {
      FeasibilityIndices(&fea, &inf, i, cone);
      if (inf < infmin) {
	infmin = inf;
	fi = fea;
	*hnext = i;
      }
    }
  }
  if (debug) {
    printf("*infeasible rays (min) =%5ld, #feas rays =%5ld\n", infmin, fi);
  }
}

void SelectNextHalfspace3(dd_ConePtr cone, dd_rowset excluded, dd_rowrange *hnext)
{
  /*Choose the next hyperplane with maximum infeasibility*/
  long i, fea, inf, infmax, fi=0;   /*feasibility and infeasibility numbers*/

  infmax = -1;
  for (i = 1; i <= cone->m; i++) {
    if (!set_member(i, excluded)) {
      FeasibilityIndices(&fea, &inf, i, cone);
      if (inf > infmax) {
	infmax = inf;
	fi = fea;
	*hnext = i;
      }
    }
  }
  if (debug) {
    printf("*infeasible rays (max) =%5ld, #feas rays =%5ld\n", infmax, fi);
  }
}

void SelectNextHalfspace4(dd_ConePtr cone, dd_rowset excluded, dd_rowrange *hnext)
{
  /*Choose the next hyperplane with the most unbalanced cut*/
  long i, fea, inf, max, tmax, fi=0, infi=0;
      /*feasibility and infeasibility numbers*/

  max = -1;
  for (i = 1; i <= cone->m; i++) {
    if (!set_member(i, excluded)) {
      FeasibilityIndices(&fea, &inf, i, cone);
      if (fea <= inf)
        tmax = inf;
      else
        tmax = fea;
      if (tmax > max) {
        max = tmax;
        fi = fea;
        infi = inf;
        *hnext = i;
      }
    }
  }
  if (!debug)
    return;
  if (max == fi) {
    printf("*infeasible rays (min) =%5ld, #feas rays =%5ld\n", infi, fi);
  } else {
    printf("*infeasible rays (max) =%5ld, #feas rays =%5ld\n", infi, fi);
  }
}

void SelectNextHalfspace5(dd_ConePtr cone, dd_rowset excluded, dd_rowrange *hnext)
{
  /*Choose the next hyperplane which is lexico-min*/
  long i, minindex;
  double *v1, *v2;

  minindex = 0;
  v1 = NULL;
  for (i = 1; i <= cone->m; i++) {
    if (!set_member(i, excluded)) {
	  v2 = cone->A[i - 1];
      if (minindex == 0) {
	    minindex = i;
	    v1=v2;
      } else if (LexSmaller(v2,v1,cone->d)) {
        minindex = i;
	    v1=v2;
      }
    }
  }
  *hnext = minindex;
}


void SelectNextHalfspace6(dd_ConePtr cone, dd_rowset excluded, dd_rowrange *hnext)
{
  /*Choose the next hyperplane which is lexico-max*/
  long i, maxindex;
  double *v1, *v2;

  maxindex = 0;
  v1 = NULL;
  for (i = 1; i <= cone->m; i++) {
    if (!set_member(i, excluded)) {
      v2= cone->A[i - 1];
      if (maxindex == 0) {
        maxindex = i;
        v1=v2;
      } else if (LexLarger(v2, v1, cone->d)) {
        maxindex = i;
        v1=v2;
     }
    }
  }
  *hnext = maxindex;
}

long Partition(dd_rowindex OV, long p, long r, dd_Amatrix A, long dmax)
{
  double *x;
  long i,j,ovi;
  
  x=A[OV[p]-1];
  i=p-1;
  j=r+1;
  while (TRUE){
    do{
      j--;
    } while (LexLarger(A[OV[j]-1],x,dmax));
    do{
      i++;
    } while (LexSmaller(A[OV[i]-1],x,dmax));
    if (i<j){
      ovi=OV[i];
      OV[i]=OV[j];
      OV[j]=ovi;
    }
    else{
      return j;
    }
  }
}

void QuickSort(dd_rowindex OV, long p, long r, dd_Amatrix A, long dmax)
{
  long q;
  
  if (p < r){
    q = Partition(OV, p, r, A, dmax);
    QuickSort(OV, p, q, A, dmax);
    QuickSort(OV, q+1, r, A, dmax);
  }
}

void LineShellingOrder(dd_rowrange m_size, dd_colrange d_size, dd_Amatrix A, dd_rowindex OV, double *z, double *d)
/* find the shelling ordering induced by a point 
   z (interior point, i.e. A z > 0) and a direction vector  d */
{
  long i,j;
  double temp1,temp2,infinity=10.0e+20;
  static double *beta;
  static long mlast=0;
  boolean localdebug=FALSE;
  
  if ( mlast<m_size ){
    if (beta!=NULL) free(beta);
    beta=(double *)calloc(m_size, sizeof *beta);
    /* initialize only for the first time or when last m_size is smaller */
    mlast=m_size;
  }
  for (i=1; i<= m_size; i++) beta[i-1]=A[i-1][0]; /* store the first column in beta */
  for (i=1; i<= m_size; i++){
    temp1 = 0.0;
    temp2 = 0.0;
    for (j = 1; j <= d_size; j++){
      temp1 += A[i - 1][j-1] * z[j-1];
      temp2 += A[i - 1][j-1] * d[j-1];
    }
/*    if (abs(temp1)>zero) A[i-1][0]=temp2/temp1;  */
    if (dd_Nonzero(temp1)) A[i-1][0]=temp2/temp1;  
    else if (temp1*temp2 > 0) A[i-1][0]= infinity;
    else A[i-1][0]= -infinity;
     /* use the first column of A tentatively */
  }
  if (localdebug) 
    for (i=1; i<= m_size; i++){
      printf("set A[%ld] = %g\n", i, A[i-1][0]);
    }
  QuickSort(OV, 1, m_size, A, 1);
  for (i=1; i<= m_size; i++) {
    A[i-1][0]=beta[i-1]; 
     /* restore the first column of A */ 
    if (localdebug) printf("restore A[%ld] with %g\n", i, A[i-1][0]);
  }
}


#ifndef RAND_MAX 
#define RAND_MAX 32767 
#endif

void RandomPermutation(dd_rowindex OV, long t, unsigned int seed)
{
  long k,j,ovj;
  double u,xk,r,rand_max=(double) RAND_MAX;
  boolean localdebug=FALSE;

  srand(seed);
  for (j=t; j>1 ; j--) {
    r=rand();
    u=r/rand_max;
    xk=j*u +1;
    k=xk;
    if (localdebug) printf("u=%g, k=%ld, r=%g, randmax= %g\n",u,k,r,rand_max);
    ovj=OV[j];
    OV[j]=OV[k];
    OV[k]=ovj;
    if (localdebug) printf("row %ld is exchanged with %ld\n",j,k); 
  }
}

void ComputeRowOrderVector(dd_ConePtr cone)
{
  long i,itemp,j;
  static dd_Arow zvec,dvec;
  static dd_colrange last_d=0;

  if (last_d < cone->d){
    if (last_d>0) {free(zvec); free(dvec);}
    zvec=(double*)calloc(cone->d,sizeof(double));
    dvec=(double*)calloc(cone->d,sizeof(double));
    last_d=cone->d;
  }
  
  cone->OrderVector[0]=0;
  switch (cone->HalfspaceOrder){
  case MaxIndex:
    for(i=1; i<=cone->m; i++) cone->OrderVector[i]=cone->m-i+1;
    break;

  case MinIndex: 
    for(i=1; i<=cone->m; i++) cone->OrderVector[i]=i;
    break;

  case LexMin: case MinCutoff: case MixCutoff: case MaxCutoff:
    for(i=1; i<=cone->m; i++) cone->OrderVector[i]=i;
    QuickSort(cone->OrderVector, 1, cone->m, cone->A, cone->d);
    break;

  case LexMax:
    for(i=1; i<=cone->m; i++) cone->OrderVector[i]=i;
    QuickSort(cone->OrderVector, 1, cone->m, cone->A, cone->d);
    for(i=1; i<=cone->m/2;i++){   /* just reverse the order */
      itemp=cone->OrderVector[i];
      cone->OrderVector[i]=cone->OrderVector[cone->m-i+1];
      cone->OrderVector[cone->m-i+1]=itemp;
    }
    break;

  case RandomRow:
    for(i=1; i<=cone->m; i++) cone->OrderVector[i]=i;
    RandomPermutation(cone->OrderVector, cone->m, cone->rseed);
    break;

  case LineShelling:
    for(i=1; i<=cone->m; i++) cone->OrderVector[i]=i;
    zvec[0]=1;
    dvec[0]=0;
    srand(cone->rseed);
    for(j=2; j<=cone->d; j++){
      zvec[j-1]=0;
      dvec[j-1]=cone->d-j+1;
      /* dvec[j-1]=rand(); */
    }
    LineShellingOrder(cone->m, cone->d, cone->A, cone->OrderVector, zvec, dvec);
    break;
  }
}

void UpdateRowOrderVector(dd_ConePtr cone, dd_rowset PriorityRows)
/* Update the RowOrder vector to shift selected rows
in highest order.
*/
{
  dd_rowrange i,j,k,j1=0,oj=0;
  long rr;
  boolean found, localdebug=FALSE;
  
  if (debug) localdebug=TRUE;
  found=TRUE;
  rr=set_card(PriorityRows);
  if (localdebug) set_write(PriorityRows);
  for (i=1; i<=rr; i++){
    found=FALSE;
    for (j=i; j<=cone->m && !found; j++){
      oj=cone->OrderVector[j];
      if (set_member(oj, PriorityRows)){
        found=TRUE;
        if (localdebug) printf("%ldth in sorted list (row %ld) is in PriorityRows\n", j, oj);
        j1=j;
      }
    }
    if (found){
      if (j1>i) {
        /* shift everything lower: ov[i]->cone->ov[i+1]..ov[j1-1]->cone->ov[j1] */
        for (k=j1; k>=i; k--) cone->OrderVector[k]=cone->OrderVector[k-1];
        cone->OrderVector[i]=oj;
        if (localdebug){
          printf("OrderVector updated to:\n");
          for (j = 1; j <= cone->m; j++) printf(" %2ld", cone->OrderVector[j]);
          printf("\n");
        }
      }
    } else {
      printf("UpdateRowOrder: Error.\n");
      goto _L99;
    }
  }
_L99:;
}

void SelectPreorderedNext(dd_ConePtr cone, dd_rowset excluded, dd_rowrange *hh)
{
  dd_rowrange i,k;
  
  *hh=0;
  for (i=1; i<=cone->m && *hh==0; i++){
    k=cone->OrderVector[i];
    if (!set_member(k, excluded)) *hh=k ;
  }
}

void SelectNextHalfspace(dd_ConePtr cone, dd_rowset excluded, dd_rowrange *hh)
{
  if (cone->PreOrderedRun){
    if (debug) {
      printf("debug SelectNextHalfspace: Use PreorderNext\n");
    }
    SelectPreorderedNext(cone, excluded, hh);
  }
  else {
    if (debug) {
      printf("debug SelectNextHalfspace: Use DynamicOrderedNext\n");
    }

    switch (cone->HalfspaceOrder) {

    case MaxIndex:
      SelectNextHalfspace0(cone, excluded, hh);
      break;

    case MinIndex:
      SelectNextHalfspace1(cone, excluded, hh);
      break;

    case MinCutoff:
      SelectNextHalfspace2(cone, excluded, hh);
      break;

    case MaxCutoff:
      SelectNextHalfspace3(cone, excluded, hh);
      break;

    case MixCutoff:
      SelectNextHalfspace4(cone, excluded, hh);
      break;

    default:
      SelectNextHalfspace0(cone, excluded, hh);
      break;
    }
  }
}

boolean dd_Nonnegative(double val)
{
  if (val>=-dd_zero) return TRUE;
  else return FALSE;
}

boolean dd_Nonpositive(double val)
{
  if (val<=dd_zero) return TRUE;
  else return FALSE;
}

boolean dd_Positive(double val)
{
  return !dd_Nonpositive(val);
}

boolean dd_Negative(double val)
{
  return !dd_Nonnegative(val);
}

boolean dd_Zero(double val)
{
  return (dd_Nonnegative(val) && dd_Nonpositive(val));
}

boolean dd_Nonzero(double val)
{
  return (dd_Positive(val) || dd_Negative(val));
}



/* end of cddarith.c */


