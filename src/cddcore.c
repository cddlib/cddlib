/* cddcore.c:  Core Procedures for cddlib
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.90e, July 12, 2000
*/

/* cddlib : C-library of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddlibman.tex for detail.
*/

#include "setoper.h"  /* set operation library header (Ver. June 1, 2000 or later) */
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
  for (i = 1; i <= (cone->parent->m); i++){
    set_addelem(cone->GroundSet, i);
    if (cone->parent->EqualityIndex[i]==1) set_addelem(cone->EqualitySet,i);
    if (cone->parent->EqualityIndex[i]==-1) set_addelem(cone->NonequalitySet,i);
  }
}


void AValue(mytype *val, dd_colrange d_size, dd_Amatrix A, mytype *p, dd_rowrange i)
{
  /*return the ith component of the vector  A x p */
  dd_colrange j;
  mytype x;

  dd_init(x); 
  dd_init(*val);
  for (j = 0; j < d_size; j++){
    dd_mul(x,A[i - 1][j], p[j]);
    dd_add(*val, *val, x);
  }
  dd_clear(x);
}

void StoreRay1(dd_ConePtr cone, mytype *p, boolean *feasible)
{  /* Original ray storing routine when RelaxedEnumeration is FALSE */
  dd_rowrange i,k,fii=cone->m+1;
  dd_colrange j;
  mytype temp;
  dd_RayPtr RR;
  boolean localdebug=debug;

  dd_init(temp);
  RR=cone->LastRay;
  *feasible = TRUE;
  set_initialize(&(RR->ZeroSet),cone->m);
  for (j = 0; j < cone->d; j++){
    dd_set(RR->Ray[j],p[j]);
  }
  for (i = 1; i <= cone->m; i++) {
    k=cone->OrderVector[i];
    AValue(&temp, cone->d, cone->A, p, k);
    if (localdebug) {
      printf("StoreRay1: AValue at row %ld =",k);
      dd_WriteNumber(stdout, temp);
      printf("\n");
    }
    if (dd_EqualToZero(temp)) {
      set_addelem(RR->ZeroSet, k);
      if (localdebug) {
        printf("recognized zero!\n");
      }
    }
    if (dd_Negative(temp)){
      if (localdebug) {
        printf("recognized negative!\n");
      }
      *feasible = FALSE;
      if (fii>cone->m) fii=i;  /* the first violating inequality index */
      if (localdebug) {
        printf("this ray is not feasible, neg comp = %ld\n", fii);
        dd_WriteNumber(stdout, temp);  printf("\n");
      }
    }
  }
  RR->FirstInfeasIndex=fii;
  RR->feasible = *feasible;
  dd_clear(temp);
}

void StoreRay2(dd_ConePtr cone, mytype *p, 
    boolean *feasible, boolean *weaklyfeasible)
   /* Ray storing routine when RelaxedEnumeration is TRUE.
       weaklyfeasible is true iff it is feasible with
       the strict_inequality conditions deleted. */
{
  dd_RayPtr RR;
  dd_rowrange i,k,fii=cone->m+1;
  dd_colrange j;
  mytype temp;
  boolean localdebug=debug;

  dd_init(temp);
  RR=cone->LastRay;
  if (debug) localdebug=TRUE;
  *feasible = TRUE;
  *weaklyfeasible = TRUE;
  set_initialize(&(RR->ZeroSet),cone->m);
  for (j = 0; j < cone->d; j++){
    dd_set(RR->Ray[j],p[j]);
  }
  for (i = 1; i <= cone->m; i++) {
    k=cone->OrderVector[i];
    AValue(&temp, cone->d, cone->A, p, k);
    if (dd_EqualToZero(temp)){
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
  dd_clear(temp);
}


void AddRay(dd_ConePtr cone, mytype *p)
{  
  boolean feasible, weaklyfeasible;
  dd_colrange j;

  if (cone->FirstRay == NULL) {
    cone->FirstRay = (dd_RayPtr) malloc(sizeof(dd_RayType));
    cone->FirstRay->Ray = (mytype *) calloc(cone->d, sizeof(mytype));
    for (j=0; j<cone->d; j++) dd_init(cone->FirstRay->Ray[j]);
    dd_init(cone->FirstRay->ARay);
    if (debug)
      printf("Create the first ray pointer\n");
    cone->LastRay = cone->FirstRay;
    cone->ArtificialRay->Next = cone->FirstRay;
  } else {
    cone->LastRay->Next = (dd_RayPtr) malloc(sizeof(dd_RayType));
    cone->LastRay->Next->Ray = (mytype *) calloc(cone->d, sizeof(mytype));
    for (j=0; j<cone->d; j++) dd_init(cone->LastRay->Next->Ray[j]);
    dd_init(cone->LastRay->Next->ARay);
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
  dd_colrange j,d1;
  boolean feasible;

  if (cone->d<=0) d1=1; else d1=cone->d;
  dd_InitializeArow(d1, &zerovector);
  if (cone->ArtificialRay != NULL) {
    printf("Warning !!!  FirstRay in not nil.  Illegal Call\n");
    free(zerovector); /* 086 */
    return;
  }
  cone->ArtificialRay = (dd_RayPtr) malloc(sizeof(dd_RayType));
  cone->ArtificialRay->Ray = (mytype *) calloc(d1, sizeof(mytype));
  for (j=0; j<d1; j++) dd_init(cone->ArtificialRay->Ray[j]);
  dd_init(cone->ArtificialRay->ARay);

  if (debug)
    printf("Create the artificial ray pointer\n");
  for (j = 0; j < d1; j++){
    dd_init(zerovector[j]);
  }
  cone->LastRay=cone->ArtificialRay;
  StoreRay1(cone, zerovector, &feasible);  
    /* This stores a vector to the record pointed by cone->LastRay */
  cone->ArtificialRay->Next = NULL;
  for (j = 0; j < d1; j++){
    dd_clear(zerovector[j]);
  }
  free(zerovector); /* 086 */
}

void ConditionalAddEdge(dd_ConePtr cone, 
    dd_RayPtr Ray1, dd_RayPtr Ray2, dd_RayPtr ValidFirstRay)
{
  long it,it_row,fii1,fii2,fmin,fmax;
  boolean adjacent,lastchance;
  dd_RayPtr TempRay,Rmin,Rmax;
  dd_AdjacencyType *NewEdge;
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
  /*  printf("Warning: CreateInitialEdges called with NULL pointer(s)\n"); */
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

void FreeDDMemory0(dd_ConePtr cone)
{
  dd_RayPtr Ptr, PrevPtr;
  long count;
  dd_rowrange i;
  boolean localdebug=FALSE;
  
  /* THIS SHOULD BE REWRITTEN carefully */
  PrevPtr=cone->ArtificialRay;
  if (PrevPtr!=NULL){
    count=0;
    for (Ptr=cone->ArtificialRay->Next; Ptr!=NULL; Ptr=Ptr->Next){
      free(PrevPtr->Ray);
      free(PrevPtr->ZeroSet);
      free(PrevPtr);
      count++;
      PrevPtr=Ptr;
    };
    cone->FirstRay=NULL;
/* must add (by Sato) */
    free(cone->LastRay->Ray);
    cone->LastRay->Ray = NULL;
    set_free(cone->LastRay->ZeroSet);
    cone->LastRay->ZeroSet = NULL;
    free(cone->LastRay);
/*    */
    cone->LastRay = NULL;
    cone->ArtificialRay=NULL;
    if (localdebug) printf("%ld ray storage spaces freed\n",count);
  }
/* must add (by Sato) */
  free(cone->Edges);
  
  set_free(cone->GroundSet); 
  set_free(cone->EqualitySet); 
  set_free(cone->NonequalitySet); 
  set_free(cone->AddedHalfspaces); 
  set_free(cone->WeaklyAddedHalfspaces); 
  set_free(cone->InitialHalfspaces);
  free(cone->InitialRayIndex);
  free(cone->OrderVector);
  free(cone->newcol);
  dd_FreeBmatrix(cone->d,cone->B);
  dd_FreeBmatrix(cone->d,cone->Bsave);

/*must replace (by Sato) */
  for (i=0; i<cone->m_alloc; i++){
    free(cone->A[i]);
    cone->A[i] = NULL;
  }

  free(cone->A);
  cone->A = NULL;

  free(cone);
}

void dd_FreeDDMemory(dd_PolyhedraPtr poly)
{
  FreeDDMemory0(poly->child);
  poly->child=NULL;
}

void dd_FreePolyhedra(dd_PolyhedraPtr poly)
{
  dd_bigrange i;

  if ((poly)->child != NULL) dd_FreeDDMemory(poly);
  dd_FreeAmatrix((poly)->m_alloc,poly->d_alloc, poly->A);
  dd_FreeArow((poly)->d_alloc,(poly)->c);
  free((poly)->EqualityIndex);
  if (poly->AincGenerated){
    for (i=1; i<=poly->m1; i++){
      set_free(poly->Ainc[i-1]);
    }
    free(poly->Ainc);
    set_free(poly->Ared);
    set_free(poly->Adom);
    poly->Ainc=NULL;
  }

  free(poly);
}

void Normalize(dd_colrange d_size, mytype *V)
{
  long j,jmin=0;
  mytype temp,min;
  boolean nonzerofound=FALSE;

  if (d_size>0){
    dd_init(min);  dd_init(temp);
    dd_abs(min,V[0]);  jmin=0; /* set the minmizer to 0 */
    if (dd_Positive(min)) nonzerofound=TRUE;
    for (j = 1; j < d_size; j++) {
      dd_abs(temp,V[j]);
      if (dd_Positive(temp)){
        if (!nonzerofound || dd_Smaller(temp,min)){
          dd_set(min, temp);  jmin=j;
        }
      }
    }
    if (dd_Positive(min)){
      for (j = 0; j < d_size; j++) dd_div(V[j], V[j], min);
    }
    dd_clear(min); dd_clear(temp);
  }
}


void ZeroIndexSet(dd_rowrange m_size, dd_colrange d_size, dd_Amatrix A, mytype *x, dd_rowset ZS)
{
  dd_rowrange i;
  mytype temp;

  set_emptyset(ZS);
  for (i = 1; i <= m_size; i++) {
    AValue(&temp, d_size, A, x, i);
/*    if (fabs(temp) < zero)*/
    if (dd_EqualToZero(temp))
      set_addelem(ZS, i);
  }
}

void CopyBmatrix(dd_colrange d_size, dd_Bmatrix T, dd_Bmatrix TCOPY)
{
  dd_rowrange i;
  dd_colrange j;

  for (i=0; i < d_size; i++) {
    for (j=0; j < d_size; j++) {
      dd_set(TCOPY[i][j],T[i][j]);
    }
  }
}


void CopyArow(mytype *acopy, mytype *a, dd_colrange d)
{
  dd_colrange j;

  for (j = 0; j < d; j++) {
    dd_set(acopy[j],a[j]);
  }
}

void CopyAmatrix(mytype **Acopy, mytype **A, dd_rowrange m, dd_colrange d)
{
  dd_rowrange i;

  for (i = 0; i< m; i++) {
    CopyArow(Acopy[i],A[i],d);
  }
}

void dd_InitializeArow(dd_colrange d,dd_Arow *a)
{
  dd_colrange j;

  if (d>0) *a=(mytype*) calloc(d,sizeof(mytype));
  for (j = 0; j < d; j++) {
      dd_init((*a)[j]);
  }
}

void dd_InitializeAmatrix(dd_rowrange m,dd_colrange d,dd_Amatrix *A)
{
  dd_rowrange i;
  dd_colrange j;

  if (m>0) (*A)=(mytype**) calloc(m,sizeof(mytype*));
  for (i = 0; i < m; i++) {
    dd_InitializeArow(d,&((*A)[i]));
  }
  for (i = 0; i < m; i++) {
    for (j = 0; j < d; j++) {
      dd_init((*A)[i][j]);
    }
  }
}

void dd_FreeAmatrix(dd_rowrange m,dd_colrange d,dd_Amatrix A)
{
  dd_rowrange i;
  dd_colrange j;

  for (i = 0; i < m; i++) {
    for (j = 0; j < d; j++) {
      dd_clear(A[i][j]);
    }
  }
  if (A!=NULL) {
    for (i = 0; i < m; i++) {
      free(A[i]);
    }
    free(A);
  }
}

void dd_FreeArow(dd_colrange d, dd_Arow a)
{
  dd_colrange j;

  for (j = 0; j < d; j++) {
    dd_clear(a[j]);
  }
  free(a);
}


void dd_InitializeBmatrix(dd_colrange d,dd_Bmatrix *B)
{
  dd_colrange i,j;

  (*B)=(mytype**) calloc(d,sizeof(mytype*));
  for (j = 0; j < d; j++) {
    (*B)[j]=(mytype*) calloc(d,sizeof(mytype));
  }
  for (i = 0; i < d; i++) {
    for (j = 0; j < d; j++) {
      dd_init((*B)[i][j]);
    }
  }
}

void dd_FreeBmatrix(dd_colrange d,dd_Bmatrix B)
{
  dd_colrange i,j;

  for (i = 0; i < d; i++) {
    for (j = 0; j < d; j++) {
      dd_clear(B[i][j]);
    }
  }
  if (B!=NULL) {
    for (j = 0; j < d; j++) {
      free(B[j]);
    }
    free(B);
  }
}

dd_SetFamilyPtr dd_CreateSetFamily(dd_bigrange fsize, dd_bigrange ssize)
{
  dd_SetFamilyPtr F;
  dd_bigrange i,f0,f1,s0,s1;

  if (fsize<=0) {
    f0=0; f1=1;  
    /* if fsize<=0, the fsize is set to zero and the created size is one */
  } else {
    f0=fsize; f1=fsize;
  }
  if (ssize<=0) {
    s0=0; s1=1;  
    /* if ssize<=0, the ssize is set to zero and the created size is one */
  } else {
    s0=ssize; s1=ssize;
  }

  F=(dd_SetFamilyPtr) malloc (sizeof(dd_SetFamilyType));
  F->set=(set_type*) calloc(f1,sizeof(set_type));
  for (i=0; i<f1; i++) {
    set_initialize(&(F->set[i]), s1);
  }
  F->famsize=f0;
  F->setsize=s0;
  return F;
}


void dd_FreeSetFamily(dd_SetFamilyPtr F)
{
  dd_bigrange i,f1;

  if (F->famsize<=0) f1=1; else f1=F->famsize; 
    /* the smallest created size is one */
  if (F!=NULL){
    for (i=0; i<f1; i++) {
      set_free(F->set[i]);
    }
    free(F->set);
    free(F);
  }
}

dd_MatrixPtr dd_CreateMatrix(dd_rowrange m_size,dd_colrange d_size)
{
  dd_MatrixPtr M;
  dd_rowrange i,m0,m1;
  dd_colrange j,d0,d1;

  if (m_size<=0){ 
    m0=0; m1=1;  
    /* if m_size <=0, the number of rows is set to zero, the actual size is 1 */
  } else {
    m0=m_size; m1=m_size;
  }
  if (d_size<=0){ 
    d0=0; d1=1;  
    /* if d_size <=0, the number of cols is set to zero, the actual size is 1 */
  } else {
    d0=d_size; d1=d_size;
  }
  M=(dd_MatrixPtr) malloc (sizeof(dd_MatrixType));
  dd_InitializeAmatrix(m1,d1,&(M->matrix));
  dd_InitializeArow(d1,&(M->rowvec));
  M->rowsize=m0;
  set_initialize(&(M->linset), m1);
  M->colsize=d0;
  for (i = 1; i <= m1; i++) {
    for (j = 1; j <= d1; j++) {
      dd_init(M->matrix[i - 1][j - 1]);
    }
  }
  for (j = 1; j <= d1; j++) {
    dd_init(M->rowvec[j - 1]);
  }
  return M;
}

void dd_FreeMatrix(dd_MatrixPtr M)
{
  dd_rowrange m1;
  dd_colrange d1;

  if (M->rowsize<=0) m1=1; else m1=M->rowsize;
  if (M->colsize<=0) d1=1; else d1=M->colsize;
  if (M!=NULL) {
    dd_FreeAmatrix(m1,d1,M->matrix);
    dd_FreeArow(d1,M->rowvec);
    set_free(M->linset);
    free(M);
  }
}

void dd_SetToIdentity(dd_colrange d_size, dd_Bmatrix T)
{
  dd_colrange j1, j2;

  for (j1 = 1; j1 <= d_size; j1++) {
    for (j2 = 1; j2 <= d_size; j2++) {
      if (j1 == j2)
        dd_set(T[j1 - 1][j2 - 1],dd_one);
      else
        dd_set(T[j1 - 1][j2 - 1],dd_purezero);
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
        for (i=1; i<=cone->m; i++) dd_set(cone->A[i-1][j1-1],cone->A[i-1][j-1]);
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
        for (i=1; i<=cone->d; i++) dd_WriteNumber(stdout, cone->B[i-1][j-1]);
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
  boolean stop, chosen, localdebug=debug;
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
  do {   /* Find a set of rows for a basis */
      SelectPivot2(cone->m, cone->d,cone->A,cone->B,cone->HalfspaceOrder,cone->OrderVector,
       cone->EqualitySet,cone->m, NopivotRow, ColSelected, &r, &s, &chosen);
      if (debug && chosen) 
        printf("Procedure FindBasis: pivot on (r,s) =(%ld, %ld).\n", r, s);
      if (chosen) {
        set_addelem(cone->InitialHalfspaces, r);
        set_addelem(NopivotRow, r);
        set_addelem(ColSelected, s);
        cone->InitialRayIndex[s]=r;    /* cone->InitialRayIndex[s] stores the corr. row index */
        (*rank)++;
        GaussianColumnPivot(cone->m, cone->d, cone->A, cone->B, r, s);
        if (localdebug) dd_WriteBmatrix(stdout,cone->d,cone->B);
      } else {
        stop=TRUE;
      }
      if (*rank==cone->d) stop = TRUE;
  } while (!stop);
  set_free(ColSelected);
  set_free(NopivotRow);
}


void FindInitialRays(dd_ConePtr cone, boolean *found)
{
  dd_rowset CandidateRows;
  dd_rowrange i;
  long rank;
  dd_RowOrderType roworder_save=LexMin;

  *found = FALSE;
  set_initialize(&CandidateRows, cone->m);
  if (cone->parent->InitBasisAtBottom==TRUE) {
    roworder_save=cone->HalfspaceOrder;
    cone->HalfspaceOrder=MaxIndex;
    cone->PreOrderedRun=FALSE;
  }
  else cone->PreOrderedRun=TRUE;
  if (debug) dd_WriteBmatrix(stdout, cone->d, cone->B);
  for (i = 1; i <= cone->m; i++)
    if (!set_member(i,cone->NonequalitySet)) set_addelem(CandidateRows, i);
    /*all rows not in NonequalitySet are candidates for initial cone*/
  FindBasis(cone, &rank);
  if (debug) dd_WriteBmatrix(stdout, cone->d, cone->B);
  if (debug) printf("FindInitialRays: rank of Amatrix = %ld\n", rank);
  cone->LinearityDim=cone->d - rank;
  if (debug) printf("Linearity Dimension = %ld\n", cone->LinearityDim);
  if (cone->LinearityDim > 0) {
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
    if (!dd_Equal((*RP1)->Ray[j - 1],(*RP2)->Ray[j - 1]))
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
  mytype a1, a2, v1, v2;
  static dd_Arow NewRay;
  static dd_colrange last_d=0;
  boolean localdebug=debug;

  dd_init(a1); dd_init(a2); dd_init(v1); dd_init(v2);
  if (last_d!=cone->d){
    if (last_d>0) {
      for (j=0; j<last_d; j++) dd_clear(NewRay[j]);
      free(NewRay);
    }
    NewRay=(mytype*)calloc(cone->d,sizeof(mytype));
    for (j=0; j<cone->d; j++) dd_init(NewRay[j]);
    last_d=cone->d;
  }

  AValue(&a1, cone->d, cone->A, Ptr1->Ray, ii);
  AValue(&a2, cone->d, cone->A, Ptr2->Ray, ii);
  if (localdebug) {
    printf("CreatNewRay: Ray1 ="); WriteArow(stdout, cone->d, Ptr1->Ray);
    printf("CreatNewRay: Ray2 ="); WriteArow(stdout, cone->d, Ptr2->Ray);
  }
  dd_abs(v1,a1);
  dd_abs(v2,a2);
  if (localdebug){
    printf("AValue1 and ABS");  dd_WriteNumber(stdout,a1); dd_WriteNumber(stdout,v1); printf("\n");
    printf("AValue2 and ABS");  dd_WriteNumber(stdout,a2); dd_WriteNumber(stdout,v2); printf("\n");
  }
  for (j = 0; j < cone->d; j++){
    dd_lincomb(NewRay[j], Ptr1->Ray[j],v2,Ptr2->Ray[j],v1);
  }
  if (localdebug) {
    printf("CreatNewRay: New ray ="); WriteArow(stdout, cone->d, NewRay);
  }
  Normalize(cone->d, NewRay);
  if (localdebug) {
    printf("CreatNewRay: Normalized ray ="); WriteArow(stdout, cone->d, NewRay);
  }
  AddRay(cone, NewRay);
  dd_clear(a1); dd_clear(a2); dd_clear(v1); dd_clear(v2);
}

void EvaluateARay1(dd_rowrange i, dd_ConePtr cone)
/* Evaluate the ith component of the vector  A x RD.Ray 
    and rearrange the linked list so that
    the infeasible rays with respect to  i  will be
    placed consecutively from First 
 */
{
  dd_colrange j;
  mytype temp,tnext;
  dd_RayPtr Ptr, PrevPtr, TempPtr;

  dd_init(temp); dd_init(tnext);
  Ptr = cone->FirstRay;
  PrevPtr = cone->ArtificialRay;
  if (PrevPtr->Next != Ptr) {
    printf("Error.  Artificial Ray does not point to FirstRay!!!\n");
  }
  while (Ptr != NULL) {
    dd_set(temp,dd_purezero);
    for (j = 0; j < cone->d; j++){
      dd_mul(tnext,cone->A[i - 1][j],Ptr->Ray[j]);
      dd_add(temp,temp,tnext);
    }
    dd_set(Ptr->ARay,temp);
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
  dd_clear(temp); dd_clear(tnext);
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
  mytype temp,tnext;
  dd_RayPtr Ptr, NextPtr;
  boolean zerofound=FALSE,negfound=FALSE,posfound=FALSE;

  dd_init(temp); dd_init(tnext);
  cone->PosHead=NULL;cone->ZeroHead=NULL;cone->NegHead=NULL;
  cone->PosLast=NULL;cone->ZeroLast=NULL;cone->NegLast=NULL;
  Ptr = cone->FirstRay;
  while (Ptr != NULL) {
    NextPtr=Ptr->Next;  /* remember the Next record */
    Ptr->Next=NULL;     /* then clear the Next pointer */
    dd_set(temp,dd_purezero);
    for (j = 0; j < cone->d; j++){
      dd_mul(tnext,cone->A[i - 1][j],Ptr->Ray[j]);
      dd_add(temp,temp,tnext);
    }
    dd_set(Ptr->ARay,temp);
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
  dd_clear(temp); dd_clear(tnext);
}

void DeleteNegativeRays(dd_ConePtr cone)
/* Eliminate the infeasible rays with respect to  i  which
   are supposed to be consecutive from the head of the dd_Ray list,
   and sort the zero list assumed to be consecutive at the
   end of the list.
 */
{
  dd_rowrange fii,fiitest;
  mytype temp;
  dd_RayPtr Ptr, PrevPtr,NextPtr,ZeroPtr1,ZeroPtr0;
  boolean found, completed, zerofound=FALSE,negfound=FALSE,posfound=FALSE;
  boolean localdebug=FALSE;
  
  dd_init(temp);
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
    dd_set(temp,Ptr->ARay);
    if (localdebug) {printf("Ptr->ARay :"); dd_WriteNumber(stdout, temp);}
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
  dd_clear(temp);
}

void FeasibilityIndices(long *fnum, long *infnum, dd_rowrange i, dd_ConePtr cone)
{
  /*Evaluate the number of feasible rays and infeasible rays*/
  /*  w.r.t the hyperplane  i*/
  dd_colrange j;
  mytype temp, tnext;
  dd_RayPtr Ptr;

  dd_init(temp); dd_init(tnext);
  *fnum = 0;
  *infnum = 0;
  Ptr = cone->FirstRay;
  while (Ptr != NULL) {
    dd_set(temp,dd_purezero);
    for (j = 0; j < cone->d; j++){
      dd_mul(tnext, cone->A[i - 1][j],Ptr->Ray[j]);
      dd_add(temp, temp, tnext);
    }
    if (temp >= 0)
      (*fnum)++;
    else
      (*infnum)++;
    Ptr = Ptr->Next;
  }
  dd_clear(temp); dd_clear(tnext);
}

boolean LexSmaller(mytype *v1, mytype *v2, long dmax)
{ /* dmax is the size of vectors v1,v2 */
  boolean determined, smaller;
  dd_colrange j;

  smaller = FALSE;
  determined = FALSE;
  j = 1;
  do {
    if (!dd_Equal(v1[j - 1],v2[j - 1])) {  /* 086 */
      if (dd_Smaller(v1[j - 1],v2[j - 1])) {  /*086 */
	    smaller = TRUE;
	  }
      determined = TRUE;
    } else
      j++;
  } while (!(determined) && (j <= dmax));
  return smaller;
}


boolean LexLarger(mytype *v1, mytype *v2, long dmax)
{
  return LexSmaller(v2, v1, dmax);
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
  double prevprogress, progress;
  mytype value1, value2;
  boolean adj, equal, completed;

  dd_init(value1); dd_init(value2);
  EvaluateARay1(hnew, cone);  
   /*Check feasibility of rays w.r.t. hnew 
     and put all infeasible ones consecutively */

  RayPtr0 = cone->ArtificialRay;   /*Pointer pointing RayPrt1*/
  RayPtr1 = cone->FirstRay;        /*1st hnew-infeasible ray to scan and compare with feasible rays*/
  dd_set(value1,cone->FirstRay->ARay);
  if (dd_Nonnegative(value1)) {
    if (cone->RayCount==cone->WeaklyFeasibleRayCount) cone->CompStatus=AllFound;
    goto _L99;        /* Sicne there is no hnew-infeasible ray and nothing to do */
  }
  else {
    RayPtr2s = RayPtr1->Next;/* RayPtr2s must point the first feasible ray */
    pos2=1;
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
    dd_set(value1,RayPtr1->ARay);
    dd_set(value2,RayPtr2->ARay);
    CheckEquality(cone->d, &RayPtr1, &RayPtr2, &equal);
    if ((dd_Positive(value1) && dd_Negative(value2)) || (dd_Negative(value1) && dd_Positive(value2))){
      CheckAdjacency(cone, &RayPtr1, &RayPtr2, &adj);
      if (adj) CreateNewRay(cone, RayPtr1, RayPtr2, hnew);
    }
    if (RayPtr2 != RayPtr3) {
      RayPtr2 = RayPtr2->Next;
      continue;
    }
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
  dd_clear(value1); dd_clear(value2);
}

void AddNewHalfspace2(dd_ConePtr cone, dd_rowrange hnew)
/* This procedure must be used under PreOrderedRun mode */
{
  dd_RayPtr RayPtr0,RayPtr1,RayPtr2;
  dd_AdjacencyType *EdgePtr, *EdgePtr0;
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
      dd_WriteNumber(stdout,RayPtr0->ARay); 
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
  mytype *v1, *v2;

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
  mytype *v1, *v2;

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
  mytype *x;
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
  long i,itemp;

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
    RandomPermutation(cone->OrderVector, cone->m, cone->rseed);
    QuickSort(cone->OrderVector, 1, cone->m, cone->A, cone->d);
    break;

  case LexMax:
    for(i=1; i<=cone->m; i++) cone->OrderVector[i]=i;
    RandomPermutation(cone->OrderVector, cone->m, cone->rseed);
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

boolean dd_Nonnegative(mytype val)
{
/*  if (val>=-dd_zero) return TRUE;  */
  if (dd_cmp(val,dd_minuszero)>=0) return TRUE;
  else return FALSE;
}

boolean dd_Nonpositive(mytype val)
{
/*  if (val<=dd_zero) return TRUE;  */
  if (dd_cmp(val,dd_zero)<=0) return TRUE;
  else return FALSE;
}

boolean dd_Positive(mytype val)
{
  return !dd_Nonpositive(val);
}

boolean dd_Negative(mytype val)
{
  return !dd_Nonnegative(val);
}

boolean dd_EqualToZero(mytype val)
{
  return (dd_Nonnegative(val) && dd_Nonpositive(val));
}

boolean dd_Nonzero(mytype val)
{
  return (dd_Positive(val) || dd_Negative(val));
}

boolean dd_Equal(mytype val1,mytype val2)
{
  return (!dd_Larger(val1,val2) && !dd_Smaller(val1,val2));
}

boolean dd_Larger(mytype val1,mytype val2)
{
  mytype temp;

  dd_init(temp);
  dd_sub(temp,val1, val2);
  return dd_Positive(temp);
  dd_clear(temp);
}

boolean dd_Smaller(mytype val1,mytype val2)
{
  return dd_Larger(val2,val1);
}

void dd_abs(mytype absval, mytype val)
{
  if (dd_Negative(val)) dd_neg(absval,val);
  else dd_set(absval,val); 
}

void dd_lincomb(mytype lc, mytype v1, mytype c1, mytype v2, mytype c2)
/*  lc := v1 * c1 + v2 * c2   */
{
  mytype temp;

  dd_init(temp);
  dd_mul(lc,v1,c1);
  dd_mul(temp,v2,c2); 
  dd_add(lc,lc,temp);
  dd_clear(temp);
}

/* end of cddcore.c */


