/* cddio.c:  Basic Input and Output Procedures for cddlib.c
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

static void fread_rational_value (FILE *, double *);
void SetLinearity(dd_PolyhedraPtr, char *);

dd_NumberType GetNumberType(char *line)
{
  dd_NumberType nt;

  if (strncmp(line, "integer", 7)==0) {
    nt = Integer;
  }
  else if (strncmp(line, "rational", 8)==0) {
    nt = Rational;
  }
  else if (strncmp(line, "real", 4)==0) {
    nt = Real;
  }
  else { 
    nt=Unknown;
  }
  return nt;
}

void ProcessCommandLine(dd_PolyhedraPtr poly, char *line)
{
  dd_colrange j;
  dd_rowrange eqsize,var;
  char *newline;


  if (strncmp(line, "hull", 4)==0) {
    poly->Representation = Generator;
    return;
  }
  if (strncmp(line, "debug", 5)==0) {
    debug = TRUE;
    return;
  }
  if ((strncmp(line, "partial_enum", 12)==0 || 
       strncmp(line, "equality", 8)==0  ||
       strncmp(line, "linearity", 9)==0 )
    && poly->RestrictedEnumeration==FALSE) {
    fgets(newline,dd_linelenmax,stdin);
    SetLinearity(poly,newline);
    return;
  }

}

void AddInequalities(dd_PolyhedraPtr poly, dd_MatrixPtr M)
{
  dd_Amatrix Anew;
  dd_rowrange i, m_new;
  dd_colrange j, d_new;

  poly->child->CompStatus=InProgress;
  
  m_new=poly->m + M->rowsize;
  d_new=poly->d_alloc;

  if (poly->m_alloc < m_new){
    dd_InitializeAmatrix(m_new,d_new,&(Anew));
    CopyAmatrix(Anew, poly->A, poly->m, poly->d);
    dd_FreeAmatrix(poly->m_alloc,&(poly->A));
    poly->A=Anew;
    poly->m_alloc=m_new;
  }
  for (i=0; i<M->rowsize; i++){
    for (j=0; j<poly->d; j++) poly->A[(poly->m)+i][j]=M->matrix[i][j];
  }
  poly->m=m_new;  
}


dd_PolyhedraPtr CreatePolyhedraData(dd_rowrange m, dd_colrange d)
{
  dd_rowrange i;
  dd_PolyhedraPtr poly=NULL;

  poly=(dd_PolyhedraPtr) malloc (sizeof(dd_Polyhedra));
  poly->child       =NULL; /* this links the homogenized cone data */
  poly->m           =m;
  poly->d           =d;  
  poly->m_alloc     =m+2; /* the allocated row size of matrix A */
  poly->d_alloc     =d;   /* the allocated col size of matrix A */
  poly->Number=Real;
  dd_InitializeAmatrix(poly->m_alloc,poly->d_alloc,&(poly->A));
  dd_InitializeArow(d,&(poly->c));           /* cost vector */
  poly->Representation       =Inequality;
  poly->Homogeneous =FALSE;

  poly->EqualityIndex=(int *)calloc(m+1, sizeof(int));  
    /* ith component is 1 if it is equality, -1 if it is strict inequality, 0 otherwise. */
  for (i = 0; i <= m; i++) poly->EqualityIndex[i]=0;

  poly->NondegAssumed           = FALSE;
  poly->InitBasisAtBottom       = FALSE;
  poly->RestrictedEnumeration   = FALSE;
  poly->RelaxedEnumeration      = FALSE;

  return poly;
}

boolean InitializeConeData(dd_rowrange m, dd_colrange d, dd_ConePtr *cone)
{
  boolean success=TRUE;
  dd_colrange j;

  (*cone)=(dd_ConePtr)calloc(1, sizeof(dd_Cone));

/* INPUT: A given representation of a cone: inequality */
  (*cone)->m=m;
  (*cone)->d=d;
  (*cone)->m_alloc=m+2; /* allocated row size of matrix A */
  (*cone)->d_alloc=d;   /* allocated col size of matrix A */
  (*cone)->Number=Real;
  (*cone)->parent=NULL;

/* CONTROL: variables to control computation */
  (*cone)->Iteration=0;

  (*cone)->HalfspaceOrder=LexMin;

  (*cone)->ArtificialRay=NULL;
  (*cone)->FirstRay=NULL;
  (*cone)->LastRay=NULL; /* The second description: Generator */
  (*cone)->PosHead=NULL;
  (*cone)->ZeroHead=NULL;
  (*cone)->NegHead=NULL;
  (*cone)->PosLast=NULL;
  (*cone)->ZeroLast=NULL;
  (*cone)->NegLast=NULL;
  (*cone)->RecomputeRowOrder  = TRUE;
  (*cone)->PreOrderedRun      = FALSE;
  set_initialize(&((*cone)->GroundSet),(*cone)->m_alloc);
  set_initialize(&((*cone)->EqualitySet),(*cone)->m_alloc);
  set_initialize(&((*cone)->NonequalitySet),(*cone)->m_alloc);
  set_initialize(&((*cone)->AddedHalfspaces),(*cone)->m_alloc);
  set_initialize(&((*cone)->WeaklyAddedHalfspaces),(*cone)->m_alloc);
  set_initialize(&((*cone)->InitialHalfspaces),(*cone)->m_alloc);
  (*cone)->RayCount=0;
  (*cone)->FeasibleRayCount=0;
  (*cone)->WeaklyFeasibleRayCount=0;
  (*cone)->TotalRayCount=0;
  (*cone)->ZeroRayCount=0;
  (*cone)->EdgeCount=0;
  (*cone)->TotalEdgeCount=0;
  (*cone)->count_int=0;
  (*cone)->count_int_good=0;
  (*cone)->count_int_bad=0;
  (*cone)->rseed=1;  /* random seed for random row permutation */
 
  dd_InitializeBmatrix((*cone)->d, &((*cone)->B));
  dd_InitializeBmatrix((*cone)->d, &((*cone)->Bsave));
  dd_InitializeAmatrix((*cone)->m_alloc,(*cone)->d_alloc,&((*cone)->A));

  (*cone)->Edges
     =(dd_Adjacency**) calloc((*cone)->m_alloc,sizeof(dd_Adjacency*));
  (*cone)->InitialRayIndex=(long*)calloc(d+1,sizeof(long));
  (*cone)->OrderVector=(long*)calloc((*cone)->m_alloc+1,sizeof(long));


  (*cone)->newcol=(long*)calloc(((*cone)->d)+1,sizeof(long));
  for (j=0; j<=(*cone)->d; j++) (*cone)->newcol[j]=j;  /* identity map, initially */
  (*cone)->LinearityDim = -2; /* -2 if it is not computed */
  (*cone)->ColReduced   = FALSE;
  (*cone)->d_orig = d;

/* STATES: variables to represent current state. */
/*(*cone)->Error;
  (*cone)->CompStatus;
  (*cone)->starttime;
  (*cone)->endtime;
*/
    
  return success;
}

dd_ConePtr ConeDataLoad(dd_PolyhedraPtr poly)
{
  dd_ConePtr cone=NULL;
  dd_colrange d,j;
  dd_rowrange m,i;

  m=poly->m;
  d=poly->d;
  if (!(poly->Homogeneous) && poly->Representation==Inequality) m=poly->m+1;

  InitializeConeData(m, d, &cone);
  cone->Representation=poly->Representation;

/* Points to the original polyhedra data, and reversely */
  cone->parent=poly;
  poly->child=cone;

  for (i=1; i<=poly->m; i++)
    for (j=1; j<=cone->d; j++)
      cone->A[i-1][j-1]=poly->A[i-1][j-1];  
  
  if (poly->Representation==Inequality && !(poly->Homogeneous)){
    cone->A[m-1][0]=1.0;
    for (j=2; j<=d; j++) cone->A[m-1][j-1]=0.0;
  }

  return cone;
}

void SetLinearity(dd_PolyhedraPtr poly, char *line)
{
  int i,j;
  dd_rowrange eqsize,var;
  char *rest;

  j=0; rest=line;
  while (line[j]==' '){
    j++;
    rest=&line[j];
  }
  eqsize=atol(rest); 
  if (debug) printf("CARDINALITY:%ld\n",eqsize);
  i=0;
  while (i < eqsize && strchr(rest,' ')!=NULL) {
     i++;
     rest=strchr(rest,' ')+1;
     var=atol(rest);
     if (debug) printf(" LINEARITY:%ld\n",var);
     poly->EqualityIndex[var]=1;
  }
  if (i==eqsize) {
    poly->RestrictedEnumeration=TRUE;
  } else {
    printf("* Warning: Linearity setting failed.\n");
  }
  return;
}

boolean dd_PolyhedraInput(dd_PolyhedraPtr *poly, dd_ErrorType *Error)
{
  dd_rowrange m_input,i;
  dd_colrange d_input,j;
  dd_RepresentationType rep=Inequality;
  double value;
  boolean found=FALSE, newformat=FALSE, successful=FALSE, linearity=FALSE;
  char command[dd_linelenmax], comsave[dd_linelenmax], numbtype[dd_wordlenmax], line[dd_linelenmax];
  dd_NumberType NT;
  dd_rowrange eqsize,var;


  (*Error)=None;
  while (!found)
  {
    if (fscanf(stdin,"%s",command)==EOF) {
      (*Error)=ImproperInputFormat;
      goto _L99;
    }
    else {
      if (strncmp(command, "V-representation", 16)==0) {
        rep=Generator; newformat=TRUE;
      }
      if (strncmp(command, "H-representation", 16)==0){
        rep=Inequality; newformat=TRUE;
      }
      if (strncmp(command, "partial_enum", 12)==0 || 
          strncmp(command, "equality", 8)==0  ||
          strncmp(command, "linearity", 9)==0 ) {
        linearity=TRUE;
        fgets(comsave,dd_linelenmax,stdin);
      }
      if (strncmp(command, "begin", 5)==0) found=TRUE;
    }
  }
  fscanf(stdin, "%ld %ld %s", &m_input, &d_input, numbtype);
  printf("size = %ld x %ld\nNumber Type = %s\n", m_input, d_input, numbtype);
  NT=GetNumberType(numbtype);
  if (NT==Unknown) {
      (*Error)=ImproperInputFormat;
      goto _L99;
    } 
  (*poly)=CreatePolyhedraData(m_input, d_input);
  (*poly)->Representation=rep;
  (*poly)->Homogeneous=TRUE;

  for (i = 1; i <= m_input; i++) {
    for (j = 1; j <= d_input; j++) {
      if (NT==Real) {
        fscanf(stdin, "%lf", &value);
      } else {
        fread_rational_value (stdin, &value);
      }
      (*poly)->A[i-1][j - 1] = value;
      if (j==1 && dd_Nonzero(value)) (*poly)->Homogeneous = FALSE;
      if (debug) printf("a(%3ld,%5ld) = %10.4f\n",i,j,value);
    }  /*of j*/
  }  /*of i*/
  if (fscanf(stdin,"%s",command)==EOF) {
   	 (*Error)=ImproperInputFormat;
  	 goto _L99;
  }
  else if (strncmp(command, "end", 3)!=0) {
     if (debug) printf("'end' missing or illegal extra data: %s\n",command);
     (*Error)=ImproperInputFormat;
     goto _L99;
  }
  
  successful=TRUE;
  if (linearity) {
    SetLinearity(*poly,comsave);
  }
  while (!feof(stdin)) {
    fscanf(stdin,"%s", command);
    ProcessCommandLine(*poly, command);
  } 

_L99: ;
  if (stdin!=NULL) fclose(stdin);
  return successful;
}

void dd_PolyhedraLoadMatrix(dd_PolyhedraPtr *poly,
dd_RepresentationType rep, dd_MatrixPtr M)
{
  dd_rowrange i;
  dd_colrange j;

  (*poly)=CreatePolyhedraData(M->rowsize, M->colsize);
  (*poly)->Representation=rep;
  (*poly)->Homogeneous=TRUE;

  for (i = 1; i <= M->rowsize; i++) {
    if (set_member(i, M->linset)) {
      (*poly)->EqualityIndex[i]=1;
    }
    for (j = 1; j <= M->colsize; j++) {
       (*poly)->A[i-1][j-1] = M->matrix[i-1][j-1];
      if (j==1 && dd_Nonzero(M->matrix[i-1][j-1])) (*poly)->Homogeneous = FALSE;
    }  /*of j*/
  }  /*of i*/
  
}

void MatrixIntegerFilter(dd_MatrixPtr M)
{   /* setting an almost integer to the integer. */
  dd_rowrange i;
  dd_colrange j;

  for (i=0; i< M->rowsize; i++)
    for (j=0; j< M->colsize; j++)
       M->matrix[i][j]=SnapToInteger(M->matrix[i][j]);    
}

void CopyRay(double *a, dd_colrange d_origsize, dd_RayPtr RR, 
  dd_RepresentationType rep, dd_colindex reducedcol)
{
  long j,j1;
  double b;

  for (j = 1; j <= d_origsize; j++){
    j1=reducedcol[j];
    if (j1>0){
      a[j-1]=RR->Ray[j1-1]; 
        /* the original column j is mapped to j1, and thus
           copy the corresponding component */
    } else {
      a[j-1]=0.;  
        /* original column is redundant and removed for computation */
    }
  }

  b=a[0];
  if (rep==Generator && dd_Nonzero(b)){
    a[0]=1.;
    for (j = 2; j <= d_origsize; j++)
       a[j-1]=a[j-1]/b;    /* normalization for generators */
  }

}

void WriteRay(FILE *f, dd_colrange d_origsize, dd_RayPtr RR, dd_RepresentationType rep, dd_colindex reducedcol)
{
  dd_colrange j;
  static dd_colrange d_last=0;
  static dd_Arow a;

  if (d_last< d_origsize){
    if (d_last>0) free(a);
    dd_InitializeArow(d_origsize+1, &a);
    d_last=d_origsize+1;
  }

  CopyRay(a, d_origsize, RR, rep, reducedcol);
  for (j = 0; j < d_origsize; j++) dd_WriteReal(f, a[j]);
  fprintf(f, "\n");
}


void dd_WriteAmatrix(FILE *f, dd_Amatrix A, long rowmax, long colmax)
{
  long i,j;

  if (A==NULL){
    fprintf(f, "WriteAmatrix: The requested matrix is empty\n");
    goto _L99;
  }
  fprintf(f, "begin\n");
  fprintf(f, " %5ld  %5ld    real\n",rowmax, colmax);
  for (i=1; i <= rowmax; i++) {
    for (j=1; j <= colmax; j++) {
      dd_WriteReal(f, A[i-1][j-1]);
    }
    fprintf(f,"\n");
  }
  fprintf(f, "end\n");
_L99:;
}

void dd_WriteBmatrix(FILE *f, dd_colrange d_size, dd_Bmatrix B)
{
  dd_colrange j1, j2;

  if (B==NULL){
    fprintf(f, "WriteBmatrix: The requested matrix is empty\n");
    goto _L99;
  }
  for (j1 = 0; j1 < d_size; j1++) {
    for (j2 = 0; j2 < d_size; j2++) {
      fprintf(f, "%15.7f ", B[j1][j2]);
    }  /*of j2*/
    putc('\n', f);
  }  /*of j1*/
  putc('\n', f);
_L99:;
}

void dd_WriteSetFamily(FILE *f, dd_SetFamilyPtr F)
{
  dd_bigrange i;

  if (F==NULL){
    fprintf(f, "WriteSetFamily: The requested family is empty\n");
    goto _L99;
  }
  fprintf(f,"begin\n");
  fprintf(f,"  %ld    %ld\n", F->famsize, F->setsize);
  for (i=0; i<F->famsize; i++) {
    fprintf(f, "  %ld : ", set_card(F->set[i]));
    set_fwrite(f, F->set[i]);
  }
  fprintf(f,"end\n");
_L99:;
}

void dd_WriteSetFamilyWithNumbers(FILE *f, dd_SetFamilyPtr F)
{
  dd_bigrange i;

  if (F==NULL){
    fprintf(f, "WriteSetFamily: The requested family is empty\n");
    goto _L99;
  }
  fprintf(f,"begin\n");
  fprintf(f,"  %ld    %ld\n", F->famsize, F->setsize);
  for (i=0; i<F->famsize; i++) {
    fprintf(f, "  %ld  %ld : ", i+1, set_card(F->set[i]));
    set_fwrite(f, F->set[i]);
  }
  fprintf(f,"end\n");
_L99:;
}

void dd_WriteMatrix(FILE *f, dd_MatrixPtr M)
{
  dd_rowrange i, linsize;

  if (M==NULL){
    fprintf(f, "WriteMmatrix: The requested matrix is empty\n");
    goto _L99;
  }
  linsize=set_card(M->linset);
  if (linsize>0) {
    fprintf(f, "linearity %ld ", linsize);
    for (i=1; i<=M->rowsize; i++) 
      if (set_member(i, M->linset)) fprintf(f, " %ld", i);
    fprintf(f, "\n");
  }
  dd_WriteAmatrix(f, M->matrix, M->rowsize, M->colsize);
_L99:;
}

double SnapToInteger(double x)
{
  long ix1,ix2,ix;

  ix1= fabs(x)* (double)dd_magfac + 0.5;
  ix2= (fabs(x) + 0.5);
  ix2= ix2*dd_magfac;
  if ( ix1 == ix2) {
    if (x>0) {
      ix = x + 0.5;
    } else {
      ix = -x + 0.5;
      ix = -ix;
    }
    return (double)ix;
  } else
    return x;
}


void dd_WriteReal(FILE *f, double x)
{
  long ix1,ix2,ix;

  ix1= fabs(x) * 10000. + 0.5;
  ix2= (fabs(x) + 0.5);
  ix2= ix2*10000;
  if ( ix1 == ix2) {
    if (x>0) {
      ix = x + 0.5;
    } else {
      ix = -x + 0.5;
      ix = -ix;
    }
    fprintf(f, " %2ld", ix);
  } else
    fprintf(f, " % .9E", x);
}


void WriteIncidence(FILE *f, dd_ConePtr cone, dd_RayPtr RR)
{
  dd_rowset cset;
  long zcar;

  set_initialize(&cset,cone->m);
  zcar = set_card(RR->ZeroSet);
  if (cone->m - zcar >= zcar) {
    fprintf(f, " %1ld : ", zcar);
    set_fwrite(f, RR->ZeroSet);
  } else {
    set_diff(cset, cone->GroundSet, RR->ZeroSet);
    fprintf(f, " %1ld : ", zcar - cone->m);
    set_fwrite(f, cset);
  }
  set_free(cset);
}


void WriteProgramDescription(FILE *f)
{
  fprintf(f, "* cdd: Double Description Method C-Code:%s\n", DDVERSION);
  fprintf(f,"* %s\n",COPYRIGHT);
}

void dd_WriteRunningMode(FILE *f, dd_PolyhedraPtr poly)
{

  switch (poly->child->HalfspaceOrder) {

    case MinIndex:
      fprintf(f, "minindex\n");
      break;

    case MaxIndex:
      fprintf(f, "maxindex\n");
      break;

    case MinCutoff:
      fprintf(f, "mincutoff\n");
      break;

    case MaxCutoff:
      fprintf(f, "maxcutoff\n");
    break;

    case MixCutoff:
      fprintf(f, "mixcutoff\n");
      break;

    case LexMin:
      fprintf(f, "lexmin\n");
      break;

    case LexMax:
      fprintf(f, "lexmax\n");
      break;

    case RandomRow:
      fprintf(f, "random  %d\n",poly->child->rseed);
      break;

    case LineShelling:
      fprintf(f, "lineshelling\n");
      break;

    default: break;
  }
}


void WriteCompletionStatus(FILE *f, dd_ConePtr cone)
{
  if (cone->Iteration<cone->m && cone->CompStatus==AllFound) {
    fprintf(f,"*Computation completed at Iteration %4ld.\n", cone->Iteration);
  } 
  if (cone->CompStatus == RegionEmpty) {
    fprintf(f,"*Computation completed at Iteration %4ld because the region found empty.\n",cone->Iteration);
  }   
}

void dd_WritePolyhedraFile(FILE *f, dd_PolyhedraPtr poly)
{
  dd_WriteAmatrix(f,poly->A,poly->m,poly->d);
}


void dd_WriteErrorMessages(FILE *f, dd_ErrorType Error)
{
  switch (Error) {

  case LowColumnRank:
      fprintf(f,"*Input Error: Input matrix (b, -A) is not column full rank => no vertices,\n");
      fprintf(f,"*Input Error: Or The polytope (convex hull) is not full dimensional.\n");
      break;
 
  case DimensionTooLarge:
    fprintf(f, "*Input Error: Input matrix is too large:\n");
    fprintf(f, "*Please increase MMAX and/or NMAX in the source code and recompile.\n");
    break;

  case FileNotFound:
    fprintf(f, "*Input Error: Specified input file does not exist.\n");
    break;

  case ImproperInputFormat:
    fprintf(f,"*Input Error: Input format is not correct.\n");
    fprintf(f,"*Format:\n");
    fprintf(f," begin\n");
    fprintf(f,"   m   n  NumberType(real, rational or integer)\n");
    fprintf(f,"   b  -A\n");
    fprintf(f," end\n");
    break;

  case EmptyVrepresentation:
    fprintf(f, "*Input Error: V-representation is empty:\n");
    fprintf(f, "*cddlib does not accept this trivial case for which output can be any inconsistent system.\n");
    break;

  case None:
    fprintf(f,"*No Error found.\n");
    break;
  }
}

dd_SetFamilyPtr dd_CopyIncidence(dd_PolyhedraPtr poly)
{
  dd_RayPtr RayPtr;
  dd_SetFamilyPtr F=NULL;
  dd_bigrange i=0;

  if (poly->child==NULL || poly->child->CompStatus!=AllFound) goto _L99;
  F=dd_CreateSetFamily(poly->child->FeasibleRayCount, poly->child->m);
  RayPtr = poly->child->FirstRay;
  while (RayPtr != NULL) {
    if (RayPtr->feasible) {
      set_copy(F->set[i], RayPtr->ZeroSet);
    }
    RayPtr = RayPtr->Next; i++;
  }
_L99:;
  return F;
}

dd_SetFamilyPtr dd_CopyAdjacency(dd_PolyhedraPtr poly)
{
  dd_RayPtr RayPtr1,RayPtr2;
  dd_SetFamilyPtr F=NULL;
  long pos1, pos2;
  boolean adj;

  if (poly->child==NULL || poly->child->CompStatus!=AllFound) goto _L99;
  F=dd_CreateSetFamily(poly->child->FeasibleRayCount, 
      poly->child->FeasibleRayCount);
  poly->child->LastRay->Next=NULL;
  for (RayPtr1=poly->child->FirstRay, pos1=1;RayPtr1 != NULL; 
				RayPtr1 = RayPtr1->Next, pos1++){
    for (RayPtr2=poly->child->FirstRay, pos2=1; RayPtr2 != NULL; 
					RayPtr2 = RayPtr2->Next, pos2++){
      if (RayPtr1!=RayPtr2){
        CheckAdjacency(poly->child, &RayPtr1, &RayPtr2, &adj);
        if (adj){
          set_addelem(F->set[pos1-1], pos2);
        }
      }
    }
  }
_L99:;
  return F;
}

dd_MatrixPtr dd_CopyGenerators(dd_PolyhedraPtr poly)
{
  dd_RayPtr RayPtr;
  dd_MatrixPtr M=NULL;
  dd_rowrange i=0,total;
  dd_colrange j,j1;
  double b;

  total=poly->child->LinearityDim + poly->child->FeasibleRayCount;
  if (poly->child->newcol[1]==0) total=total-1;
  if (poly->child==NULL || poly->child->CompStatus!=AllFound) goto _L99;
  if (poly->Representation==Inequality){
    M=dd_CreateMatrix(total, poly->d);

    RayPtr = poly->child->FirstRay;
    while (RayPtr != NULL) {
      if (RayPtr->feasible) {
        CopyRay(M->matrix[i], poly->d, RayPtr, Generator, poly->child->newcol);
      }
      RayPtr = RayPtr->Next; i++;
    }
    for (j=2; j<=poly->d; j++){
      if (poly->child->newcol[j]==0){
        b=poly->child->Bsave[0][j-1];
        if (dd_Nonzero(b)){
          M->matrix[i][0]=1.0;  /* Normalize */
          for (j1=1; j1<poly->d; j1++) 
            M->matrix[i][j1]=(poly->child->Bsave[j1][j-1])/b;
        } else {
          for (j1=0; j1<poly->d; j1++)
            M->matrix[i][j1]=poly->child->Bsave[j1][j-1];
        }
        set_addelem(M->linset, i+1);
        i++;
      }     
    }
  } else {
    M=dd_CreateMatrix(poly->m, poly->d);
    CopyAmatrix(M->matrix, poly->A, poly->m, poly->d);
    for (i=1; i<=poly->m; i++) 
      if (poly->EqualityIndex[i]==1) set_addelem(M->linset,i);
  }
  MatrixIntegerFilter(M);
_L99:;
  return M;
}

dd_MatrixPtr dd_CopyInequalities(dd_PolyhedraPtr poly)
{
  dd_RayPtr RayPtr;
  dd_MatrixPtr M=NULL;
  dd_rowrange i=0,total;
  dd_colrange j,j1;

  total=poly->child->LinearityDim + poly->child->FeasibleRayCount;
  if (poly->child->newcol[1]==0) total=total-1;
  if (poly->child==NULL || poly->child->CompStatus!=AllFound) goto _L99;
  if (poly->Representation==Generator){
    M=dd_CreateMatrix(total, poly->d);

    RayPtr = poly->child->FirstRay;
    while (RayPtr != NULL) {
      if (RayPtr->feasible) {
        CopyRay(M->matrix[i], poly->d, RayPtr, Inequality, poly->child->newcol);
      }
      RayPtr = RayPtr->Next; i++;
    }
    for (j=2; j<=poly->d; j++){
      if (poly->child->newcol[j]==0){
        for (j1=0; j1<poly->d; j1++) M->matrix[i][j1]=poly->child->Bsave[j1][j-1];
        set_addelem(M->linset, i+1);
        i++;
      }     
    }
  } else {
    M=dd_CreateMatrix(poly->m, poly->d);
    CopyAmatrix(M->matrix, poly->A, poly->m, poly->d);
    for (i=1; i<=poly->m; i++) 
      if (poly->EqualityIndex[i]==1) set_addelem(M->linset,i);
  }
  MatrixIntegerFilter(M);
_L99:;
  return M;
}

/****************************************************************************************/
/*  rational number (a/b) read is taken from Vinci by Benno Bueeler and Andreas Enge    */
/****************************************************************************************/
typedef double rational;

void sread_rational_value (char *s, rational *value)
   /* reads a rational value from the specified string "s" and assigns it to "value"    */
   
{
   char     *numerator_s, *denominator_s, *position, token;
   int      sign = 1, i;
   rational numerator, denominator;
   
   /* determine the sign of the number */
   numerator_s = s;
   if (s [0] == '-')
   {  sign = -1;
      numerator_s++;
   }
   else if (s [0] == '+')
      numerator_s++;
      
   /* look for a sign '/' and eventually split the number in numerator and denominator */
   position = strchr (numerator_s, '/');
   if (position != NULL)
   {  *position = '\0'; /* terminates the numerator */
      denominator_s = position + 1;
   };

   /* determine the floating point values of numerator and denominator */
   numerator = 0;
   for (i = 0; i < strlen (numerator_s); i++)
   {  token = numerator_s [i];
      if (strchr ("0123456789", token)) /* token is a cypher */
         numerator = 10 * numerator + atoi (&token);
   }
   
   if (position != NULL)
   {  denominator = 0;
      for (i = 0; i < strlen (denominator_s); i++)
      {  token = denominator_s [i];
         if (strchr ("0123456789", token)) /* token is a cypher */
            denominator = 10 * denominator + atoi (&token);
      }    
   }
   else denominator = 1;
   
   *value = sign * numerator / denominator;
}
   

static void fread_rational_value (FILE *f, double *value)
   /* reads a rational value from the specified file "f" and assigns it to "value"      */
   
{
   char     number_s [255];
   rational rational_value;
   
   fscanf (f, "%s ", number_s);
   sread_rational_value (number_s, &rational_value);
   *value = rational_value;
   
}
   
/****************************************************************************************/


/* end of cddio.c */

