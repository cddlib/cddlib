/* dplex.c:  dual simplex method c-code
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.80, March 13, 1999
*/

/* dplex.c : C-Implementation of the dual simplex method for
   solving an LP: max/min  A_(m-1).x subject to  x in P, where
   P= {x :  A_i.x >= 0, i=0,...,m-2, and  x_0=1}, and
   A_i is the i-th row of an m x n matrix A.
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#include "setoper.h"  /* set operation library header (Ver. March 16,1995 or later) */
#include "dplex.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define COPYRIGHT   "Copyright (C) 1997, Komei Fukuda, fukuda@ifor.math.ethz.ch"
#define DPLEXVERSION   "Version 0.62dev (June 1998)"

#define dp_FALSE 0
#define dp_TRUE 1

typedef set_type rowset;  /* set_type defined in setoper.h */
typedef set_type colset;

typedef struct lpdata {
  dp_FilenameType filename;
  dp_LPConversionType conv;
  dp_LPSolverType solver; 
  dp_rowrange m;
  dp_colrange d;
  dp_Amatrix A;
  dp_Bmatrix B;
  dp_rowrange objrow;
  dp_colrange rhscol;
  dp_NumberType number;

  dp_LPStatusType LPS;  /* the current solution status */
  dp_rowrange m_alloc; /* the allocated row size of matrix A */
  dp_colrange d_alloc; /* the allocated col size of matrix A */
  double optvalue;  /* optimal value */
  dp_Arow sol;   /* primal solution */
  dp_Arow dsol;  /* dual solution */
  dp_colindex nbindex;  /* current basis represented by nonbasic indices */
  dp_rowrange re;  /* row index as a certificate in the case of inconsistency */
  dp_colrange se;  /* col index as a certificate in the case of dual inconsistency */
  long anticycle_iter;
  long phase1_iter;
  long phase2_iter;
  long total_iter;
  int use_given_basis;  /* switch to indicate the use of the given basis */
  dp_colindex given_nbindex;  /* given basis represented by nonbasic indices */
} dp_LPDataType;


typedef enum {
  dp_MaxIndex,dp_MinIndex,dp_LexMin,dp_LexMax,dp_RandomRow,dp_LineShelling
} dp_HyperplaneOrderType;

void dp_InitializeArow(dp_colrange,dp_Arow *);
void dp_InitializeAmatrix(dp_rowrange,dp_colrange,dp_Amatrix*);
void dp_InitializeBmatrix(dp_colrange,dp_Bmatrix*);
void dp_FreeAmatrix(dp_rowrange,dp_Amatrix *);
void dp_FreeBmatrix(dp_colrange,dp_Bmatrix *);
void CrissCrossMinimize(dp_rowrange,dp_colrange,dp_Amatrix,dp_Bmatrix T,
  dp_rowrange,dp_colrange,int,
  dp_LPStatusType *,double *optvalue,dp_Arow,dp_Arow,dp_colindex,
  dp_rowrange *,dp_colrange *,long *,dp_ErrorType *);
void CrissCrossMaximize(dp_rowrange,dp_colrange,dp_Amatrix,dp_Bmatrix T,
  dp_rowrange,dp_colrange,int,
  dp_LPStatusType *,double *optvalue,dp_Arow,dp_Arow,dp_colindex,
  dp_rowrange *,dp_colrange *,long *,dp_ErrorType *);
void DualSimplexMinimize(dp_rowrange,dp_colrange, dp_Amatrix,dp_Bmatrix,
  dp_rowrange,dp_colrange,int,
  dp_LPStatusType *,double *,dp_Arow,dp_Arow,dp_colindex,
  dp_rowrange *,dp_colrange *,long *,dp_ErrorType *);
void DualSimplexMaximize(dp_rowrange,dp_colrange,dp_Amatrix,dp_Bmatrix,
  dp_rowrange,dp_colrange,int,
  dp_LPStatusType *,double *,dp_Arow,dp_Arow,dp_colindex,
  dp_rowrange *,dp_colrange *,long *,dp_ErrorType *);
void FindLPBasis(dp_rowrange,dp_colrange,dp_Amatrix,dp_Bmatrix,dp_rowindex,
    dp_colindex,dp_rowindex,dp_rowrange,dp_colrange,int,
    dp_colrange *,int *,dp_LPStatusType *,long *);
void FindDualFeasibleBasis(dp_rowrange,dp_colrange,dp_Amatrix,dp_Bmatrix,dp_rowindex,
    dp_colindex,long *,dp_rowrange,dp_colrange,
    dp_colrange *,int *,dp_LPStatusType *,long *);

void dp_WriteBmatrix(FILE *f,dp_colrange d_size,dp_Bmatrix T);
void dp_SetNumberType(char *line,dp_NumberType *number,dp_ErrorType *Error);
void dp_ComputeRowOrderVector(dp_rowrange m_size,dp_colrange d_size,dp_Amatrix A,
    dp_rowindex OV,dp_HyperplaneOrderType ho,unsigned int rseed);
void dp_SelectPreorderedNext(dp_rowrange m_size,dp_colrange d_size,
    rowset excluded,dp_rowindex OV,dp_rowrange *hnext);
void dp_WriteReal(FILE *f,double x);
void dp_SetInputFile(FILE **f,dp_FilenameType inputfile, dp_ErrorType *);
void SetSolutions(dp_rowrange,dp_colrange,
   dp_Amatrix,dp_Bmatrix,dp_rowrange,dp_colrange,dp_LPStatusType,
   double *,dp_Arow,dp_Arow,dp_colindex,dp_rowrange,dp_colrange);


dp_LPSolutionPtr dp_LPSolutionLoad(dp_LPPtr lp)
{
  dp_LPSolutionPtr lps;
  dp_colrange j;
  long i;

  lps=(dp_LPSolutionPtr) calloc(1,sizeof(dp_LPSolutionType));
  for (i=1; i<=dp_filenamelen; i++) lps->filename[i-1]=lp->filename[i-1];
  lps->conv=lp->conv;
  lps->solver=lp->solver; 
  lps->m=lp->m;
  lps->d=lp->d;
  lps->number=lp->number;

  lps->LPS=lp->LPS;  /* the current solution status */
  lps->optvalue=lp->optvalue;  /* optimal value */
  lps->sol=(double*) calloc(lp->d,sizeof(double));   /* primal solution */
  lps->dsol=(double*) calloc(lp->d,sizeof(double));  /* dual solution */
  lps->nbindex=(long*) calloc((lp->d)+1,sizeof(long));  /* dual solution */
  for (j=1; j<=lp->d; j++){
    lps->sol[j-1]=lp->sol[j-1];
    lps->dsol[j-1]=lp->dsol[j-1];
    lps->nbindex[j]=lp->nbindex[j];
  }
  lps->anticycle_iter=lp->anticycle_iter;
  lps->phase1_iter=lp->phase1_iter;
  lps->phase2_iter=lp->phase2_iter;
  lps->total_iter=lp->total_iter;

  return lps;
}

void dp_SetInputFile(FILE **f,dp_FilenameType inputfile,dp_ErrorType *Error)
{
  int opened=0,stop,quit=0;
  int i,dotpos=0,trial=0;
  char ch;
  char *tempname;
  
  
  *Error=dp_None;
  while (!opened && !quit) {
    printf("\n>> Input file (*.ine) : ");
    scanf("%s",inputfile);
    ch=getchar();
    stop=dp_FALSE;
    for (i=0; i<dp_filenamelen && !stop; i++){
      ch=inputfile[i];
      switch (ch) {
        case '.': 
          dotpos=i+1;
          break;
        case ';':  case ' ':  case '\0':  case '\n':  case '\t':     
          stop=dp_TRUE;
          tempname=(char*)calloc(dp_filenamelen,sizeof(ch));
          strncpy(tempname,inputfile,i);
          strcpy(inputfile,tempname);
          break;
      }
    }
    if ( ( *f = fopen(inputfile,"r") )!= NULL) {
      printf("input file %s is open\n",inputfile);
      opened=1;
      *Error=dp_None;
    }
    else{
      printf("The file %s not found\n",inputfile);
      trial++;
      if (trial>5) {
        *Error=dp_FileNotFound;
        quit=1;
      }
    }
  }
}


void LPLoadInit(dp_LPPtr lp, dp_LPConversionType CV,
   dp_NumberType NUMB,dp_rowrange m,dp_colrange d,
   dp_rowrange OBJ,dp_colrange RHS)
{
  lp->solver=dp_DualSimplex;  /* set the default lp solver */
  lp->d=d;
  lp->m=m;
  lp->number=NUMB;
  lp->objrow=OBJ;
  lp->rhscol=RHS;
  lp->LPS=dp_LPSundecided;

  lp->sol=(double*) calloc(d,sizeof(double));
  lp->dsol=(double*) calloc(d,sizeof(double));
  lp->nbindex=(long*) calloc(d+1,sizeof(long));
  lp->given_nbindex=(long*) calloc(d+1,sizeof(long));
}

dp_LPPtr dp_LPLoad(dp_LPConversionType CV,
   dp_NumberType NUMB,dp_rowrange m,dp_colrange d,dp_Amatrix A,
   dp_rowrange OBJ,dp_colrange RHS,dp_ErrorType *err)
{
  dp_LPPtr lp=NULL;
  dp_rowrange i;
  dp_colrange j;

  lp=(dp_LPPtr) calloc(1,sizeof(dp_LPDataType));
  LPLoadInit(lp,CV,NUMB,m,d,OBJ,RHS);
  lp->m_alloc=lp->m+2;
  lp->d_alloc=lp->d+2;
  switch (CV) {
    case dp_LPmax:  case dp_LPmin:
      lp->conv=CV;break;
    default: lp->conv=dp_LPmax;
  }
  dp_InitializeBmatrix(lp->d_alloc,&(lp->B));

  dp_InitializeAmatrix(lp->m_alloc,lp->d_alloc,&(lp->A));
  for (i=0; i<m; i++) {
    for (j=0; j<d; j++) lp->A[i][j]=A[i][j];   /* copy the A matrix */
  }
  return lp;
}

dp_LPPtr dp_LPDirectLoad(dp_LPConversionType CV,
   dp_NumberType NUMB,dp_rowrange m,dp_colrange d,dp_Amatrix *A,
   dp_rowrange OBJ,dp_colrange RHS,dp_ErrorType *err)
{
/* This loads an LP with given data and ERASES the given matrix A at
   the same time.  This saves one duplicate space of maxtrix A but
   must be used with causion.  Use dp_alloc  when space is not
   a problem.
*/
  dp_LPPtr lp=NULL;
  dp_rowrange i;

  lp=(dp_LPPtr) calloc(1,sizeof(dp_LPDataType));
  LPLoadInit(lp,CV,NUMB,m,d,OBJ,RHS);
  lp->m_alloc=lp->m+2;
  lp->d_alloc=lp->d+2;
  switch (CV) {
    case dp_LPmax:  case dp_LPmin:
      lp->conv=CV;break;
    default: lp->conv=dp_LPmax;
  }

  dp_InitializeBmatrix(lp->d_alloc,&(lp->B));

  dp_InitializeAmatrix(lp->m_alloc,0,&(lp->A));
  for (i=0; i<m; i++) {
      lp->A[i]=(*A)[i];   /* direct link to the i-th row of A */
  } 
  for (i=m; i<lp->m_alloc; i++) {
     dp_InitializeArow(d,&(lp->A[i]));  /* create a new row space of the rest */
  }
  free(*A);  /* remove the space created for *A */
  return lp;
}


void dp_FreeLPData(dp_LPPtr *lp)
{
  if ((*lp)!=NULL){
    dp_FreeAmatrix((*lp)->m_alloc,&((*lp)->A));
    dp_FreeBmatrix((*lp)->d_alloc,&((*lp)->B));
    free((*lp)->sol);
    free((*lp)->dsol);
    free((*lp)->nbindex);
    free((*lp)->given_nbindex);
    free(*lp);
  }
}

void dp_FreeLPSolution(dp_LPSolutionPtr *lps)
{
  if ((*lps)!=NULL){
    free((*lps)->sol);
    free((*lps)->dsol);
    free(*lps);
  }
}

int dp_LPReverseRow(dp_LPPtr lp, dp_rowrange i)
{
  dp_colrange j;
  int success=0;

  if (i>=1 && i<=lp->m){
    lp->LPS=dp_LPSundecided;
    for (j=1; j<=lp->d; j++) {
      lp->A[i-1][j-1]=-lp->A[i-1][j-1];
      /* negating the i-th constraint of A */
    }
    success=1;
  }
  return success;
}

int dp_LPReplaceRow(dp_LPPtr lp, dp_rowrange i, dp_Arow a)
{
  dp_colrange j;
  int success=0;

  if (i>=1 && i<=lp->m){
    lp->LPS=dp_LPSundecided;
    for (j=1; j<=lp->d; j++) {
      lp->A[i-1][j-1]=a[j-1];
      /* replacing the i-th constraint by a */
    }
    success=1;
  }
  return success;
}

dp_Arow dp_LPCopyRow(dp_LPPtr lp, dp_rowrange i)
{
  dp_colrange j;
  dp_Arow a;

  if (i>=1 && i<=lp->m){
    dp_InitializeArow(lp->d, &a);
    for (j=1; j<=lp->d; j++) {
      a[j-1]=lp->A[i-1][j-1];
      /* copying the i-th row to a */
    }
  }
  return a;
}

dp_LPPtr dp_LPCopy(dp_LPPtr lp)
{
  dp_LPPtr lpcopy=NULL;
  dp_rowrange i;
  dp_colrange j;

  lpcopy=(dp_LPPtr) calloc(1,sizeof(dp_LPDataType));
  LPLoadInit(lpcopy,lp->conv,lp->number,lp->m,lp->d,lp->objrow,lp->rhscol);
  lpcopy->m_alloc=lpcopy->m+2;
  lpcopy->d_alloc=lpcopy->d+2;
  dp_InitializeBmatrix(lpcopy->d_alloc,&(lpcopy->B));

  dp_InitializeAmatrix(lpcopy->m_alloc,lpcopy->d_alloc,&(lpcopy->A));
  for (i=0; i<lp->m; i++) {
    for (j=0; j<lp->d; j++) lpcopy->A[i][j]=lp->A[i][j];   /* copy the A matrix */
  }
  return lpcopy;
}


dp_LPPtr dp_LPInput(FILE **f,dp_ErrorType *err)
{
  dp_LPPtr lp=NULL;
  long i,j;
  dp_rowrange m_input,m;
  dp_colrange d_input,d;
  dp_NumberType NUMB;
  dp_LPConversionType CONV;
  dp_Amatrix A;
  dp_rowrange OBJ; 
  dp_colrange RHS;
  dp_FilenameType filename;
  double value,cost;
  int found=0,localdebug=0;
  char command[dp_wordlenmax],numbtype[dp_wordlenmax],line[dp_linelenmax];


  *err=dp_None;

  dp_SetInputFile(f,filename,err);

  if (*err!=dp_None){
    goto _L99;
  }

  while (!found)
  {
    if (fscanf(*f,"%s",command)==EOF) {
      *err=dp_ImproperInputFormat;
      goto _L99;
    }
    else if (strncmp(command,"begin",5)==0) {
      found=dp_TRUE;
    }
  }
  fscanf(*f,"%ld %ld %s",&m_input,&d_input,numbtype);
  printf("size = %ld x %ld\nNumber Type = %s\n",m_input,d_input,numbtype);
  m=m_input+1;  /* an additional row for the objective function */
  d=d_input;
  dp_SetNumberType(numbtype,&(NUMB),err);
  if (NUMB==dp_Unknown || NUMB== dp_Rational) {
      goto _L99;
  }
/*  The size restriction is eliminated from cdd-062
  if (N+1 > dp_NMAX || M > dp_MMAX) {
    *err = dp_DimensionTooLarge;
    goto _L99;
  }
*/

  dp_InitializeAmatrix(m,d,&A);
  for (i = 1; i <= m-1; i++) {
    for (j = 1; j <= d; j++) {
      fscanf(*f,"%lf",&value);
      A[i-1][j-1] = value;
      if (localdebug) printf("a(%3ld,%5ld) = %10.4f\n",i,j,value);
    }  /*of j*/
    fgets(line,dp_linelenmax,*f);
    if (localdebug) printf("comments to be skipped: %s\n",line);
    if (localdebug) putchar('\n');
  }  /*of i*/
  if (fscanf(*f,"%s",command)==EOF) {
   	 *err=dp_ImproperInputFormat;
  	 goto _L99;
  }
  else if (strncmp(command,"end",3)!=0) {
     if (localdebug) printf("'end' missing or illegal extra data: %s\n",command);
   	 *err=dp_ImproperInputFormat;
  	 goto _L99;
  }
  
  OBJ=m;
  RHS=1L;
  found=0;
  CONV=dp_LPmax;

  while (!found)
  {
    if (fscanf(*f,"%s",command)==EOF) {
      *err=dp_ImproperInputFormat;
      goto _L99;
    }
    if (strncmp(command,"maximize",8)==0) {
      CONV=dp_LPmax;
      found=dp_TRUE;
    }
    if (strncmp(command,"minimize",8)==0) {
      CONV=dp_LPmin;
      found=dp_TRUE;
    }
  }
  for (j=1;j<=d;j++) {
    fscanf(*f,"%lf",&cost);
    A[m-1][j-1]=cost;
    if (localdebug) printf(" cost[%ld] = %.9E\n",j,A[m-1][j]);
  }
/*  Direct but risky way to load an lp.  The next removes A as well. */
  lp=dp_LPDirectLoad(CONV,NUMB,m,d,&A,OBJ,RHS,err); 

/*  Safe way to load an lp */
/*  lp=dp_LPLoad(CONV,NUMB,m,d,A,OBJ,RHS,err); */
/*  dp_FreeAmatrix(m,&A); */
  
_L99: ;
  if (*f!=NULL) fclose(*f);
  return lp;
}

int dp_Nonnegative(double val)
{
  if (val>=-dp_zero) return dp_TRUE;
  else return dp_FALSE;
}

int dp_Nonpositive(double val)
{
  if (val<=dp_zero) return dp_TRUE;
  else return dp_FALSE;
}

int dp_Positive(double val)
{
  return !dp_Nonpositive(val);
}

int dp_Negative(double val)
{
  return !dp_Nonnegative(val);
}

int dp_Zero(double val)
{
  return (dp_Nonnegative(val) && dp_Nonpositive(val));
}

int dp_Nonzero(double val)
{
  return (dp_Positive(val) || dp_Negative(val));
}


void dp_SetNumberType(char *line,dp_NumberType *number,dp_ErrorType *Error)
{
  if (strncmp(line,"integer",7)==0) {
    *number = dp_Integer;
    return;
  }
  else if (strncmp(line,"rational",8)==0) {
    *number = dp_Rational;
    *Error=dp_ImproperInputFormat;  /* Rational Input not supported */
    return;
  }
  else if (strncmp(line,"real",4)==0) {
    *number = dp_Real;
    return;
  }
  else { 
    *number=dp_Unknown;
    *Error=dp_ImproperInputFormat;
  }
}

void dp_WriteErrorMessages(FILE *f,dp_ErrorType Error)
{
  switch (Error) {
    case dp_DimensionTooLarge:
      fprintf(f,"dp_Error: LP size is too large.  Modify dp_NMAX and/or dp_MMAX.\n");
      break;
      
    case dp_LowColumnRank:
      fprintf(f,"dp_Error: The matrix A is not column full rank.\n");
      break;
 
    case dp_ImproperInputFormat:
      fprintf(f,"dp_Error: Input file format is not correct.\n");
      break;

    case dp_FileNotFound:
      fprintf(f,"dp_Error: The input file does not exist.\n");
      break;

    case dp_LPLoadAmatrixError:
      fprintf(f,"dp_Error: LPLoad detects A-matrix has a wrong size.\n");
      break;
    
    case dp_None:
      fprintf(f,"dp_Error: No error occured.\n");
      break;

    default:
      fprintf(f,"dp_Error: Unknown error found.\n");
      break;
  }
}


double dp_TableauEntry(dp_rowrange m_size,dp_colrange d_size,dp_Amatrix A,dp_Bmatrix T,
				dp_rowrange r,dp_colrange s)
/* Compute the (r,s) entry of A.T   */
{
  dp_colrange j;
  double val;
  
  val=0;
  for (j=0; j< d_size; j++) {
    val = val + A[r-1][j] * T[j][s-1];
  }
  return val;
}

void dp_WriteTableau(FILE *f,dp_rowrange m_size,dp_colrange d_size,dp_Amatrix A,dp_Bmatrix T,
  dp_colindex nbindex,dp_rowindex bflag)
/* Write the tableau  A.T   */
{
  dp_colrange j;
  dp_rowrange i;
  
  fprintf(f,"  %ld   %ld    real\n",m_size,d_size);
  fprintf(f,"          |");
  for (j=1; j<= d_size; j++) {
    fprintf(f," %12ld",nbindex[j]);
  } fprintf(f,"\n");
  for (j=1; j<= d_size+1; j++) {
    fprintf(f," ------------");
  } fprintf(f,"\n");
  for (i=1; i<= m_size; i++) {
    fprintf(f," %3ld(%3ld) |",i,bflag[i]);  
    for (j=1; j<= d_size; j++) {
      fprintf(f," %12.3f",dp_TableauEntry(m_size,d_size,A,T,i,j));
    }
    fprintf(f,"\n");
  }
  fprintf(f,"end\n");
}


void SelectDualSimplexPivot(dp_rowrange m_size,dp_colrange d_size,
    int Phase1,dp_Amatrix A,dp_Bmatrix T,dp_rowindex OV,
    dp_colindex nbindex,dp_rowindex bflag,
    dp_rowrange objrow,dp_colrange rhscol,
    dp_rowrange *r,dp_colrange *s,int *selected,dp_LPStatusType *lps)
{ 
  /* selects a dual simplex pivot (*r,*s) if the current
     basis is dual feasible and not optimal. If not dual feasible,
     the procedure returns *selected=dp_FALSE and *lps=LPSundecided.
     If Phase1=dp_TRUE, the RHS column will be considered as the negative
     of the column of the largest variable (==m_size).  For this case, it is assumed
     that the caller used the auxiliary row (with variable m_size) to make the current
     dictionary dual feasible before calling this routine so that the nonbasic
     column for m_size corresponds to the auxiliary variable.
  */
  int colselected=dp_FALSE,rowselected=dp_FALSE,
    dualfeasible=dp_TRUE,localdebug=dp_FALSE;
  dp_rowrange i;
  dp_colrange j;
  double val=0,minval=0,rat=0,minrat=0;
  static dp_Arow rcost;
  static dp_colrange lastnn=0;

  if (lastnn<d_size) {
    if (rcost!=NULL) free(rcost);
    rcost=(double*) calloc(d_size,sizeof(double));
  }
  lastnn=d_size;

  *r=0; *s=0;
  *selected=dp_FALSE;
  *lps=dp_LPSundecided;
  for (j=1; j<=d_size; j++){
    if (j!=rhscol){
      rcost[j-1]=dp_TableauEntry(m_size,d_size,A,T,objrow,j);
      if (dp_Positive(rcost[j-1])) { 
        dualfeasible=dp_FALSE;
      }
    }
  }
  if (dualfeasible){
    while ((*lps==dp_LPSundecided) && (!rowselected) && (!colselected)) {
      for (i=1; i<=m_size; i++) {
        if (i!=objrow && bflag[i]==-1) {  /* i is a basic variable */
          if (Phase1){
            val=-dp_TableauEntry(m_size,d_size,A,T,i,bflag[m_size]); 
            /* for dual Phase I */
          } 
          else {val=dp_TableauEntry(m_size,d_size,A,T,i,rhscol);}
          if (val < minval) {
            *r=i;
            minval=val;
          }
        }
      }
      if (dp_Nonnegative(minval)) {
        *lps=dp_Optimal;
      }
      else {
        rowselected=dp_TRUE;
        for (j=1; j<=d_size; j++){
          val=dp_TableauEntry(m_size,d_size,A,T,*r,j);
          if (j!=rhscol && dp_Positive(val)) {
            rat=-rcost[j-1]/val;
            if (*s==0 || rat < minrat){
              minrat=rat;
              *s=j;
            }
          }
        }
        if (*s>0) {colselected=dp_TRUE; *selected=dp_TRUE;}
        else *lps=dp_Inconsistent;
      }
    } /* end of while */
  }
  if (localdebug) {
     if (Phase1) printf("Phase 1 : select %ld,%ld\n",*r,*s);
     else printf("Phase 2 : select %ld,%ld\n",*r,*s);
  }
}

void dp_SelectPivot2(dp_rowrange m_size,dp_colrange d_size,dp_Amatrix A,dp_Bmatrix T,
            dp_HyperplaneOrderType roworder,dp_rowindex ordervec,
            dp_rowrange rowmax,rowset NopivotRow,
            colset NopivotCol,dp_rowrange *r,dp_colrange *s,
            int *selected)
/* Select a position (*r,*s) in the matrix A.T such that (A.T)[*r][*s] is nonzero
   The choice is feasible, i.e., not on NopivotRow and NopivotCol, and
   best with respect to the specified roworder 
 */
{
  int stop;
  dp_rowrange r_val;
  rowset rowexcluded;
  double Xtemp;

  stop = dp_FALSE;
  set_initialize(&rowexcluded,m_size);
  set_copy(rowexcluded,NopivotRow);
  for (r_val=rowmax+1;r_val<=m_size;r_val++) {
    set_addelem(rowexcluded,r_val);   /* cannot pivot on any row > rmax */
  }
  *selected = dp_FALSE;
  do {
    r_val=0;
    dp_SelectPreorderedNext(m_size,d_size,rowexcluded,ordervec,&r_val);
    if (r_val>=1) {
      *r=r_val;
      *s=1;
      while (*s <= d_size && !*selected) {
        Xtemp=dp_TableauEntry(m_size,d_size,A,T,*r,*s);
        if (!set_member(*s,NopivotCol) && dp_Nonzero(Xtemp)) {
          *selected = dp_TRUE;
          stop = dp_TRUE;
        } else {
          (*s)++;
        }
      }
      if (!*selected) {
        set_addelem(rowexcluded,r_val);
      }
    }
    else {
      *r = 0;
      *s = 0;
      stop = dp_TRUE;
    }
  } while (!stop);
  set_free(rowexcluded);
}

void dp_GaussianColumnPivot(dp_rowrange m_size,dp_colrange d_size,
    dp_Amatrix A,dp_Bmatrix T,dp_colindex nbindex,dp_rowindex bflag,dp_rowrange r,dp_colrange s)
/* Update the Transformation matrix T with the pivot operation on (r,s) 
   This procedure performs a implicit pivot operation on the matrix A by
   updating the dual basis inverse  T.
 */
{
  int localdebug=dp_FALSE;
  long j,j1,entering;
  static dp_Arow Rtemp;
  double Xtemp0,Xtemp;
  static dp_colrange lastnn=0;

  if (lastnn<d_size) {
    if (Rtemp!=NULL) free(Rtemp);
    Rtemp=(double*) calloc(d_size,sizeof(double));
  }
  lastnn=d_size;

  for (j=1; j<=d_size; j++) Rtemp[j-1]=dp_TableauEntry(m_size,d_size,A,T,r,j);
  Xtemp0 = Rtemp[s-1];
  for (j = 1; j <= d_size; j++) {
    if (j != s) {
      Xtemp = Rtemp[j-1];
      for (j1 = 1; j1 <= d_size; j1++)
        T[j1-1][j-1] -= T[j1-1][s - 1] * Xtemp / Xtemp0;
    }
  }
  for (j = 1; j <= d_size; j++) T[j-1][s - 1] /= Xtemp0;

  entering=nbindex[s];
  bflag[r]=s;     /* the nonbasic variable r corresponds to column s */
  nbindex[s]=r;   /* the nonbasic variable on s column is r */
  if (localdebug) {
    fprintf(stdout,"Column pivot: (leaving, entering) = (%ld, %ld)\n", r,entering);
    dp_WriteBmatrix(stdout,d_size,T);
    dp_WriteTableau(stdout,m_size,d_size,A,T,nbindex,bflag);
  }

  if (entering>0) bflag[entering]=-1;
     /* original variables have negative index and should not affect the row index */

}

void dp_InitializeArow(dp_colrange d,dp_Arow *a)
{
  if (d>0) *a=(double*) calloc(d,sizeof(double));
}

void dp_InitializeAmatrix(dp_rowrange m,dp_colrange d,dp_Amatrix *A)
{
  dp_rowrange i;

  (*A)=(double**) calloc(m,sizeof(double*));
  for (i = 0; i < m; i++) {
    dp_InitializeArow(d,&((*A)[i]));
  }
}

void dp_FreeAmatrix(dp_rowrange m,dp_Amatrix *A)
{
  dp_rowrange i;

  for (i = 0; i < m; i++) {
    free((*A)[i]);
  }
  free(*A);
}

void dp_InitializeBmatrix(dp_colrange d,dp_Bmatrix *B)
{
  dp_colrange j;

  (*B)=(double**) calloc(d,sizeof(double*));
  for (j = 0; j < d; j++) {
    (*B)[j]=(double*) calloc(d,sizeof(double));
  }
}

void dp_FreeBmatrix(dp_colrange d,dp_Bmatrix *B)
{
  dp_colrange j;

  for (j = 0; j < d; j++) {
    free((*B)[j]);
  }
  free(*B);
}

void dp_SetToIdentity(dp_colrange d,dp_Bmatrix T)
{
  dp_colrange j1,j2;

  for (j1 = 1; j1 <= d; j1++) {
    for (j2 = 1; j2 <= d; j2++) {
      if (j1 == j2)
        T[j1 - 1][j2 - 1] = 1.0;
      else
        T[j1 - 1][j2 - 1] = 0.0;
    }
  }
}

void ResetTableau(dp_rowrange m_size,dp_colrange d_size,dp_Bmatrix T,
    dp_colindex nbindex,dp_rowindex bflag,dp_rowrange objrow,dp_colrange rhscol,
   int UsePrevBasis)
{
  dp_rowrange i;
  dp_colrange j;
  
  if (!UsePrevBasis) {  /* Initialize T and nbindex */
    for (j=1; j<=d_size; j++) nbindex[j]=-j;
    nbindex[rhscol]=0; 
      /* RHS is already in nonbasis and is considered to be associated
         with the zero-th row of input. */
     dp_SetToIdentity(d_size,T);
  }
  
  /* Set the bflag according to nbindex */
  for (i=1; i<=m_size; i++) bflag[i]=-1;  
    /* all basic variables have index -1 */
  bflag[objrow]= 0; 
    /* bflag of the objective variable is 0,
       different from other basic variables which have -1 */
  for (j=1; j<=d_size; j++) if (nbindex[j]>0) bflag[nbindex[j]]=j;
    /* bflag of a nonbasic variable is its column number */

}

void dp_SelectCrissCrossPivot(dp_rowrange m_size,dp_colrange d_size,dp_Amatrix A,dp_Bmatrix T,
    dp_rowindex bflag,dp_rowrange objrow,dp_colrange rhscol,
    dp_rowrange *r,dp_colrange *s,
    int *selected,dp_LPStatusType *lps)
{
  int colselected=dp_FALSE,rowselected=dp_FALSE;
  dp_rowrange i;
  double val;
  
  *selected=dp_FALSE;
  *lps=dp_LPSundecided;
  while ((*lps==dp_LPSundecided) && (!rowselected) && (!colselected)) {
    for (i=1; i<=m_size; i++) {
      if (i!=objrow && bflag[i]==-1) {  /* i is a basic variable */
        val=dp_TableauEntry(m_size,d_size,A,T,i,rhscol);
        if (dp_Negative(val)) {
          rowselected=dp_TRUE;
          *r=i;
          break;
        }
      }
      else if (bflag[i] >0) { /* i is nonbasic variable */
        val=dp_TableauEntry(m_size,d_size,A,T,objrow,bflag[i]);
        if (dp_Positive(val)) {
          colselected=dp_TRUE;
          *s=bflag[i];
          break;
        }
      }
    }
    if  ((!rowselected) && (!colselected)) {
      *lps=dp_Optimal;
      return;
    }
    else if (rowselected) {
     for (i=1; i<=m_size; i++) {
       if (bflag[i] >0) { /* i is nonbasic variable */
          val=dp_TableauEntry(m_size,d_size,A,T,*r,bflag[i]);
          if (dp_Positive(val)) {
            colselected=dp_TRUE;
            *s=bflag[i];
            *selected=dp_TRUE;
            break;
          }
        }
      }
    }
    else if (colselected) {
      for (i=1; i<=m_size; i++) {
        if (i!=objrow && bflag[i]==-1) {  /* i is a basic variable */
          val=dp_TableauEntry(m_size,d_size,A,T,i,*s);
          if (dp_Negative(val)) {
            rowselected=dp_TRUE;
            *r=i;
            *selected=dp_TRUE;
            break;
          }
        }
      }
    }
    if (!rowselected) {
      *lps=dp_DualInconsistent;
    }
    else if (!colselected) {
      *lps=dp_Inconsistent;
    }
  }
}

void CrissCrossMinimize(dp_rowrange m_size,dp_colrange d_size,
    dp_Amatrix A,dp_Bmatrix T,
    dp_rowrange objrow,dp_colrange rhscol,int UsePrevBasis,dp_LPStatusType *LPS,
    double *optvalue,dp_Arow sol,dp_Arow dsol,dp_colindex nbindex,
    dp_rowrange *re,dp_colrange *se,long *iter,dp_ErrorType *err)
{
   dp_colrange j;

   *err=dp_None;
   for (j=1; j<=d_size; j++)
     A[objrow-1][j-1]=-A[objrow-1][j-1];
   CrissCrossMaximize(m_size,d_size,A,T,objrow,rhscol,
     UsePrevBasis,LPS,optvalue,sol,dsol,nbindex,re,se,iter,err);
   *optvalue=-*optvalue;
   for (j=1; j<=d_size; j++){
     dsol[j-1]=-dsol[j-1];
     A[objrow-1][j-1]=-A[objrow-1][j-1];
   }
}

void CrissCrossMaximize(dp_rowrange m_size,dp_colrange d_size,
    dp_Amatrix A,dp_Bmatrix T,
    dp_rowrange objrow,dp_colrange rhscol,int UsePrevBasis,dp_LPStatusType *LPS,
    double *optvalue,dp_Arow sol,dp_Arow dsol,dp_colindex nbindex,
    dp_rowrange *re,dp_colrange *se,long *iter,dp_ErrorType *err)
/* 
When LP is inconsistent then *re returns the evidence row.
When LP is dual-inconsistent then *se returns the evidence column.
*/
{
  int stop,chosen,found;
  long pivots_p0;
  dp_rowrange r;
  dp_colrange s;
  static dp_rowindex bflag;
  static long mlast=0;
  static dp_rowindex OrderVector;  /* the permutation vector to store a preordered row indeces */
  unsigned int rseed=1;

  *err=dp_None;
  if (bflag==NULL || mlast!=m_size){
     if (mlast!=m_size) {
       free(bflag);   /* called previously with different m_size */
       free(OrderVector);
     }
     bflag=(long *) calloc(m_size+1,sizeof(long*));
     OrderVector=(long *)calloc(m_size+1,sizeof(long*)); 
     /* initialize only for the first time or when a larger space is needed */
     mlast=m_size;
  }
  /* Initializing control variables. */
  dp_ComputeRowOrderVector(m_size,d_size,A,OrderVector,dp_MinIndex,rseed);

  *re=0; *se=0; *iter=0;

  ResetTableau(m_size,d_size,T,nbindex,bflag,
    objrow,rhscol,UsePrevBasis);

  FindLPBasis(m_size,d_size,A,T,
      OrderVector,nbindex,bflag,
      objrow,rhscol,UsePrevBasis,&s,&found,LPS,&pivots_p0);
  *iter+=pivots_p0;

  if (!found){
     *se=s;
     goto _L99;
     /* No LP basis is found, and thus Inconsistent.  
     Output the evidence column. */
  }

  stop=dp_FALSE;
  do {   /* Criss-Cross Method */
    dp_SelectCrissCrossPivot(m_size,d_size,A,T,bflag,
       objrow,rhscol,&r,&s,&chosen,LPS);
    if (chosen) {
      dp_GaussianColumnPivot(m_size,d_size,A,T,nbindex,bflag,r,s);
      (*iter)++;
    } else {
      switch (*LPS){
        case dp_Inconsistent: *re=r;
        case dp_DualInconsistent: *se=s;
        default: break;
      }
      stop=dp_TRUE;
    }
  } while(!stop);
  
_L99:

  SetSolutions(m_size,d_size,A,T,
   objrow,rhscol,*LPS,optvalue,sol,dsol,nbindex,*re,*se);

}


int dp_LexSmaller(double *v1,double *v2,long dmax)
{ /* dmax is the size of vectors v1,v2 */
  int determined,smaller;
  dp_colrange j;

  smaller = dp_FALSE;
  determined = dp_FALSE;
  j = 1;
  do {
    if (dp_Nonzero(v1[j - 1]-v2[j - 1])) {
      if (v1[j - 1] < v2[j - 1]) {
	    smaller = dp_TRUE;
	  }
      determined = dp_TRUE;
    } else
      j++;
  } while (!(determined) && (j <= dmax));
  return smaller;
}

int dp_LexLarger(double *v1,double *v2,long dmax)
{ /* dmax is the size of vectors v1,v2 */
  int determined,larger;
  dp_colrange j;

  larger = dp_FALSE;
  determined = dp_FALSE;
  j = 1;
  do {
    if (dp_Nonzero(v1[j - 1]-v2[j - 1])) {
      if (v1[j - 1] > v2[j - 1]) {
	    larger = dp_TRUE;
	  }
      determined = dp_TRUE;
    } else
      j++;
  } while (!(determined) && (j <= dmax));
  return larger;
}

void FindLPBasis(dp_rowrange m_size,dp_colrange d_size,
    dp_Amatrix A, dp_Bmatrix T,dp_rowindex OV,dp_colindex nbindex,
    dp_rowindex bflag,dp_rowrange objrow,dp_colrange rhscol,int useprevbasis,
    dp_colrange *cs,int *found,dp_LPStatusType *lps,long *pivot_no)
{ 
  /* Find a LP basis using Gaussian pivots.
     If the problem has an LP basis,
     the procedure returns *found=dp_TRUE,*lps=LPSundecided and an LP basis.
     If the constraint matrix A (excluding the rhs and objective) is not
     column indepent, there are two cases.  If the dependency gives a dual
     inconsistency, this returns *found=dp_FALSE, *lps=dp_StrucDualInconsistent and 
     the evidence column *s.  Otherwise, this returns *found=dp_TRUE, 
     *lps=LPSundecided and an LP basis of size less than d_size.  Columns j
     that do not belong to the basis (i.e. cannot be chosen as pivot because
     they are all zero) will be indicated in nbindex vector: nbindex[j] will
     be negative and set to -j.
  */
  int localdebug=dp_FALSE,chosen,stop;
  long pivots_p0=0,rank;
  colset ColSelected;
  rowset RowSelected;
  double val;

  dp_rowrange i,r;
  dp_colrange j,s;

  *found=dp_FALSE; *cs=0; rank=0;
  *lps=dp_LPSundecided;

  set_initialize(&RowSelected,m_size);
  set_initialize(&ColSelected,d_size);
  set_addelem(RowSelected,objrow);
  set_addelem(ColSelected,rhscol);
 
  stop=dp_FALSE;
  do {   /* Find a LP basis */
    dp_SelectPivot2(m_size,d_size,A,T,dp_MinIndex,OV,
      m_size,RowSelected,ColSelected,&r,&s,&chosen);
    if (chosen) {
      set_addelem(RowSelected,r);
      set_addelem(ColSelected,s);
      dp_GaussianColumnPivot(m_size,d_size,A,T,nbindex,bflag,r,s);
      pivots_p0++;
      rank++;
    } else {
      for (j=1;j<=d_size  && *lps==dp_LPSundecided; j++) {
        if (j!=rhscol && nbindex[j]<0){
          val=dp_TableauEntry(m_size,d_size,A,T,objrow,j);
          if (dp_Nonzero(val)){  /* dual inconsistent */
            *lps=dp_StrucDualInconsistent;
            *cs=j;
            /* dual inconsistent because the nonzero reduced cost */
          }
        }
      }
      if (*lps==dp_LPSundecided) *found=dp_TRUE;  
         /* dependent columns but not dual inconsistent. */
      stop=dp_TRUE;
    }
    if (rank==d_size-1) {
      stop = dp_TRUE;
      *found=dp_TRUE;
    }
  } while (!stop);

/* Check whether the objrow is in the basis in case of UsePrevBasis. */
  if (useprevbasis && (s=bflag[objrow])>0){ /* objrow in the nonbasis. */
    for (j=0;j<=d_size;j++){
      if (j!=s) set_addelem(ColSelected,j);
      if ((i=nbindex[j])>0) set_addelem(RowSelected,i);
    }
    if (localdebug) printf("UsePrevBasis but the current basis does not contain objrow.\n");
    dp_SelectPivot2(m_size,d_size,A,T,dp_MinIndex,OV,
      m_size,RowSelected,ColSelected,&r,&s,&chosen);
    if (localdebug && chosen) printf("Procedure FindBasis: pivot on (r,s) =(%ld, %ld).\n",r,s);
    if (chosen) {
      set_addelem(RowSelected,r);
      set_addelem(ColSelected,s);
      dp_GaussianColumnPivot(m_size,d_size,A,T,nbindex,bflag,r,s);
      pivots_p0++;
    } else {
      /* objrow is nonbasic and the corresponding column is zero. */
      *found=dp_FALSE;
      *lps=dp_StrucDualInconsistent;
      *cs=s;
    }
  }

  *pivot_no=pivots_p0;
  set_free(RowSelected);
  set_free(ColSelected);
}

void FindDualFeasibleBasis(dp_rowrange m_size,dp_colrange d_size,
    dp_Amatrix A,dp_Bmatrix T,dp_rowindex OV,
    dp_colindex nbindex,dp_rowindex bflag,dp_rowrange objrow,dp_colrange rhscol,
    dp_colrange *s,int *found,dp_LPStatusType *lps,long *pivot_no)
{ 
  /* Find a dual feasible basis using Phase I of Dual Simplex method.
     If the problem is dual feasible,
     the procedure returns *found=dp_TRUE, *lps=LPSundecided and a dual feasible
     basis.   If the problem is dual infeasible, this returns
     *found=dp_FALSE, *lps=DualInconsistent and the evidence column *s.
     Caution: matrix A must have at least one extra row:  the row space A[m_size] must
     have been allocated.
  */
  int phase1,dualfeasible=dp_TRUE,localdebug=dp_FALSE,chosen,stop;
  dp_LPStatusType LPSphase1;
  long pivots_p1=0;
  dp_rowrange i,r_val;
  dp_colrange j,l,ms=0,s_val,local_m_size;
  double val=0,purezero=0,maxcost=-1;
  static dp_colrange lastnn=0;
  static dp_Arow rcost;

  if (lastnn<d_size) {
    if (rcost!=NULL) free(rcost);
    rcost=(double*) calloc(d_size,sizeof(double));
  }
  lastnn=d_size;

  *found=dp_TRUE; *lps=dp_LPSundecided; *s=0;
  local_m_size=m_size+1;  /* increase m_size by 1 */

  ms=0;  /* ms will be the index of column which has the largest reduced cost */
  for (j=1; j<=d_size; j++){
    if (j!=rhscol){
      if (localdebug) printf("checking the column %ld var %ld\n",j,nbindex[j]); 
      rcost[j-1]=dp_TableauEntry(local_m_size,d_size,A,T,objrow,j);
      if (localdebug) printf("reduced cost =  %f\n",rcost[j-1]); 
      if (rcost[j-1] > maxcost) {maxcost=rcost[j-1]; ms = j;}
    }
  }
  if (dp_Positive(maxcost)) dualfeasible=dp_FALSE;

  if (!dualfeasible){
    for (j=1; j<=d_size; j++){
      A[local_m_size-1][j-1]=purezero;
      for (l=1; l<=d_size; l++){
        if (nbindex[l]>0) {
          A[local_m_size-1][j-1]-=A[nbindex[l]-1][j-1]; 
          /* To make the auxiliary row (0,-1,-1,...,-1).  */
        }
      }
    }
    if (localdebug){
      printf("FindDualFeasibleBasis: curruent basis is not dual feasible.\n");
      printf("because of the column %ld assoc. with var %ld   dual cost =%f\n",
       ms,nbindex[ms],maxcost);
    }

    /* Pivot on (local_m_size,ms) so that the dual basic solution becomes feasible */
    dp_GaussianColumnPivot(local_m_size,d_size,A,T,nbindex,bflag,local_m_size,ms);
    pivots_p1=pivots_p1+1;

    phase1=dp_TRUE; stop=dp_FALSE;
    do {   /* Dual Simplex Phase I */
      chosen=dp_FALSE; LPSphase1=dp_LPSundecided;
      SelectDualSimplexPivot(local_m_size,d_size,phase1,A,T,OV,nbindex,bflag,
        objrow,rhscol,&r_val,&s_val,&chosen,&LPSphase1);
      if (!chosen) {
        /* The current dictionary is terminal.  There are two cases:
           dp_TableauEntry(local_m_size,d_size,A,T,objrow,ms) is negative or zero.
           The first case implies dual infeasible,
           and the latter implies dual feasible but local_m_size is still in nonbasis.
           We must pivot in the auxiliary variable local_m_size. 
        */

        double minval=0;
        r_val=0;
        for (i=1; i<=local_m_size; i++){
          if (bflag[i]<0) { 
             /* i is basic and not the objective variable */
            val=dp_TableauEntry(local_m_size,d_size,A,T,i,ms);  /* auxiliary column*/
            if (val < minval) {
              r_val=i;
              minval=val;
              if (localdebug) printf("update minval with = %f  r_val = %ld\n",minval,r_val);
            }
          }
        }

        dp_GaussianColumnPivot(local_m_size,d_size,A,T,nbindex,bflag,r_val,ms);
        pivots_p1=pivots_p1+1;

        if (dp_Negative(dp_TableauEntry(local_m_size,d_size,A,T,objrow,ms))){
          if (localdebug){
            printf("Dual infeasible.\n");
            printf("obj-ms: %f  dp_zero = %f\n",
              dp_TableauEntry(local_m_size,d_size,A,T,objrow,ms),dp_zero);
          }
          *found=dp_FALSE; *lps=dp_DualInconsistent;  *s=ms;
        }
        stop=dp_TRUE;
      } else {
        dp_GaussianColumnPivot(local_m_size,d_size,A,T,nbindex,bflag,r_val,s_val);
        pivots_p1=pivots_p1+1;
        if (bflag[local_m_size]<0) {
          stop=dp_TRUE; 
          if (localdebug) 
            printf("Dual Phase I: the auxiliary variable entered the basis, go to phase II\n");
        }
      }
    } while(!stop);
  }
  *pivot_no=pivots_p1;
}

void DualSimplexMinimize(dp_rowrange m_size,dp_colrange d_size,
   dp_Amatrix A,dp_Bmatrix T,
   dp_rowrange objrow,dp_colrange rhscol,int UsePrevBasis,dp_LPStatusType *LPS,
   double *optvalue,dp_Arow sol,dp_Arow dsol,dp_colindex nbindex,
   dp_rowrange *re,dp_colrange *se,long *iter,dp_ErrorType *err)
{
   dp_colrange j;

   *err=dp_None;
   for (j=1; j<=d_size; j++)
     A[objrow-1][j-1]=-A[objrow-1][j-1];
   DualSimplexMaximize(m_size,d_size,A,T,objrow,rhscol,UsePrevBasis,
     LPS,optvalue,sol,dsol,nbindex,re,se,iter,err);
   *optvalue=-*optvalue;
   for (j=1; j<=d_size; j++){
     dsol[j-1]=-dsol[j-1];
     A[objrow-1][j-1]=-A[objrow-1][j-1];
   }
}

void DualSimplexMaximize(dp_rowrange m_size,dp_colrange d_size,
   dp_Amatrix A,dp_Bmatrix T,
   dp_rowrange objrow,dp_colrange rhscol,int UsePrevBasis,dp_LPStatusType *LPS,
   double *optvalue,dp_Arow sol,dp_Arow dsol,dp_colindex nbindex,
   dp_rowrange *re,dp_colrange *se,long *iter,dp_ErrorType *err)
/* 
When LP is inconsistent then *re returns the evidence row.
When LP is dual-inconsistent then *se returns the evidence column.
*/
{
  int stop,chosen,phase1,found;
  long pivots_ds=0,pivots_p0=0,pivots_p1=0,pivots_pc=0,maxpivots,maxpivfactor=70;
  dp_rowrange r;
  dp_colrange s;
  static dp_rowindex bflag;
  static long mlast=0,nlast=0;
  int localdebug=dp_FALSE;
  static dp_rowindex OrderVector;  /* the permutation vector to store a preordered row indeces */
  unsigned int rseed=1;
  
  *err=dp_None;
  maxpivots=maxpivfactor*d_size;  /* maximum pivots to be performed before cc pivot is applied. */
  if (mlast!=m_size || nlast!=d_size){
     if (mlast>0) { /* called previously with different m_size */
       free(OrderVector);
       free(bflag);
     }
     OrderVector=(long *)calloc(m_size+1,sizeof(*OrderVector));
     bflag=(long *) calloc(m_size+2,sizeof(*bflag));  /* one more element for an auxiliary variable  */
     mlast=m_size;nlast=d_size;
  }
  /* Initializing control variables. */
  dp_ComputeRowOrderVector(m_size,d_size,A,OrderVector,dp_MinIndex,rseed);

  *re=0; *se=0; *iter=0;
  
  ResetTableau(m_size,d_size,T,nbindex,bflag,objrow,rhscol,UsePrevBasis);
   
  FindLPBasis(m_size,d_size,A,T,OrderVector,nbindex,bflag,
      objrow,rhscol,UsePrevBasis,&s,&found,LPS,&pivots_p0);
  *iter+=pivots_p0;

  if (!found){
     *se=s;
     goto _L99;
     /* No LP basis is found, and thus Inconsistent.  
     Output the evidence column. */
  }

  FindDualFeasibleBasis(m_size,d_size,A,T,OrderVector,nbindex,bflag,
      objrow,rhscol,&s,&found,LPS,&pivots_p1);
  *iter+=pivots_p1;

  if (!found){
     *se=s;
     goto _L99;
     /* No dual feasible basis is found, and thus DualInconsistent.  
     Output the evidence column. */
  }
  
  /* Dual Simplex Method */
  stop=dp_FALSE;
  do {
    chosen=dp_FALSE; *LPS=dp_LPSundecided; phase1=dp_FALSE;
    if (pivots_ds<maxpivots) {
      SelectDualSimplexPivot(m_size,d_size,phase1,A,T,OrderVector,nbindex,bflag,
        objrow,rhscol,&r,&s,&chosen,LPS);
    }
    if (chosen) pivots_ds=pivots_ds+1;
    if (!chosen && *LPS==dp_LPSundecided) {  
      /* In principle this should not be executed because we already have dual feasibility
         attained and dual simplex pivot should have been chosen.  This might occur
         under floating point computation, or the case of cycling.
      */
      dp_SelectCrissCrossPivot(m_size,d_size,A,T,bflag,
        objrow,rhscol,&r,&s,&chosen,LPS);
      if (chosen) pivots_pc=pivots_pc+1;
    }
    if (chosen) {
      dp_GaussianColumnPivot(m_size,d_size,A,T,nbindex,bflag,r,s);
      (*iter)++;
    } else {
      switch (*LPS){
        case dp_Inconsistent: *re=r;
        case dp_DualInconsistent: *se=s;
        default: break;
      }
      stop=dp_TRUE;
    }
  } while(!stop);

_L99:

  if (localdebug){
     printf("LP solved with %ld pivots. (ds pivt#= %ld,  p1 piv#= %ld",*iter,pivots_ds,pivots_p1);
     if (pivots_pc > 0) printf(", cc piv#= %ld",pivots_pc);
     printf(")\n");
  }
  
  SetSolutions(m_size,d_size,A,T,
   objrow,rhscol,*LPS,optvalue,sol,dsol,nbindex,*re,*se);

}

void SetSolutions(dp_rowrange m_size,dp_colrange d_size,
   dp_Amatrix A,dp_Bmatrix T,
   dp_rowrange objrow,dp_colrange rhscol,dp_LPStatusType LPS,
   double *optvalue,dp_Arow sol,dp_Arow dsol,dp_colindex nbindex,
   dp_rowrange re,dp_colrange se)
/* 
Assign the solution vectors to sol,dsol,*optvalue after solving
the LP.
*/
{
  dp_colrange j;
  double sw;
  int localdebug=dp_FALSE;
  
  switch (LPS){
  case dp_Optimal:
    for (j=1;j<=d_size; j++) {
      sol[j-1]=T[j-1][rhscol-1];
      dsol[j-1]=-dp_TableauEntry(m_size,d_size,A,T,objrow,j);
      *optvalue=dp_TableauEntry(m_size,d_size,A,T,objrow,rhscol);
      if (localdebug) printf("dsol[%ld]= %f\n",nbindex[j],dsol[j-1]);
    }
    break;
  case dp_Inconsistent:
    if (localdebug) printf("DualSimplexSolve: LP is inconsistent.\n");
    for (j=1;j<=d_size; j++) {
      sol[j-1]=T[j-1][rhscol-1];
      dsol[j-1]=-dp_TableauEntry(m_size,d_size,A,T,re,j);
      if (localdebug)  printf("dsol[%ld]= %f\n",nbindex[j],dsol[j-1]);
    }
    break;
  case dp_DualInconsistent:
    for (j=1;j<=d_size; j++) {
      sol[j-1]=T[j-1][se-1];
      dsol[j-1]=-dp_TableauEntry(m_size,d_size,A,T,objrow,j);
      if (localdebug)  printf("dsol[%ld]= %f\n",nbindex[j],dsol[j-1]);
    }
    if (localdebug) printf( "DualSimplexSolve: LP is dual inconsistent.\n");
    break;

  case dp_StrucDualInconsistent:
    if (dp_Positive(dp_TableauEntry(m_size,d_size,A,T,objrow,se))) sw=1;
    else sw=-1;
    for (j=1;j<=d_size; j++) {
      sol[j-1]=sw*T[j-1][se-1];
      dsol[j-1]=-dp_TableauEntry(m_size,d_size,A,T,objrow,j);
      if (localdebug)  printf("dsol[%ld]= %f\n",nbindex[j],dsol[j-1]);
    }
    if (localdebug) printf( "DualSimplexSolve: LP is dual inconsistent.\n");
    break;

  default:break;
  }
}


long dp_Partition(dp_rowindex OV,long p,long r,dp_Amatrix A,long dmax)
{
  double *x;
  long i,j,ovi;
  
  x=A[OV[p]-1];
  i=p-1;
  j=r+1;
  while (dp_TRUE){
    do{
      j--;
    } while (dp_LexLarger(A[OV[j]-1],x,dmax));
    do{
      i++;
    } while (dp_LexSmaller(A[OV[i]-1],x,dmax));
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

void dp_QuickSort(dp_rowindex OV,long p,long r,dp_Amatrix A,long dmax)
{
  long q;
  
  if (p < r){
    q = dp_Partition(OV,p,r,A,dmax);
    dp_QuickSort(OV,p,q,A,dmax);
    dp_QuickSort(OV,q+1,r,A,dmax);
  }
}

void dp_LineShellingOrder(dp_rowrange m_size,dp_colrange d_size,dp_Amatrix A,dp_rowindex OV,double *z,double *d)
/* find the shelling ordering induced by a point 
   z (interior point, i.e. A z > 0) and a direction vector  d */
{
  long i,j;
  double temp1,temp2,infinity=10.0e+20;
  static double *beta;
  static long mlast=0;
  int localdebug=dp_FALSE;
  
  if ( mlast<m_size ){
    if (beta!=NULL) free(beta);
    beta=(double *)calloc(m_size,sizeof(*beta));
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
    if (dp_Nonzero(temp1)) A[i-1][0]=temp2/temp1;  
    else if (temp1*temp2 > 0) A[i-1][0]= infinity;
    else A[i-1][0]= -infinity;
     /* use the first column of A tentatively */
  }
  if (localdebug) 
    for (i=1; i<= m_size; i++){
      printf("set A[%ld] = %g\n",i,A[i-1][0]);
    }
  dp_QuickSort(OV,1,m_size,A,1);
  for (i=1; i<= m_size; i++) {
    A[i-1][0]=beta[i-1]; 
     /* restore the first column of A */ 
    if (localdebug) printf("restore A[%ld] with %g\n",i,A[i-1][0]);
  }
}


#ifndef RAND_MAX 
#define RAND_MAX 32767 
#endif

void dp_RandomPermutation(dp_rowindex OV,long t,unsigned int seed)
{
  long k,j,ovj;
  double u,xk,r,rand_max=(double) RAND_MAX;
  int localdebug=dp_FALSE;

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

void dp_ComputeRowOrderVector(dp_rowrange m_size,dp_colrange d_size,dp_Amatrix A,
    dp_rowindex OV,dp_HyperplaneOrderType ho,unsigned int rseed)
{
  long i,itemp,j;
  static dp_Arow zvec,dvec;
  static dp_colrange lastnn=0;

  if (lastnn<d_size) {
    if (zvec!=NULL) {free(zvec); free(dvec);}
    zvec=(double*) calloc(d_size,sizeof(double));
    dvec=(double*) calloc(d_size,sizeof(double));
  }
  lastnn=d_size;
  
  OV[0]=0;
  switch (ho){
  case dp_MaxIndex:
    for(i=1; i<=m_size; i++) OV[i]=m_size-i+1;
    break;

  case dp_MinIndex: 
    for(i=1; i<=m_size; i++) OV[i]=i;
    break;

  case dp_LexMin:
    for(i=1; i<=m_size; i++) OV[i]=i;
    dp_QuickSort(OV,1,m_size,A,d_size);
    break;

  case dp_LexMax:
    for(i=1; i<=m_size; i++) OV[i]=i;
    dp_QuickSort(OV,1,m_size,A,d_size);
    for(i=1; i<=m_size/2;i++){   /* just reverse the order */
      itemp=OV[i];
      OV[i]=OV[m_size-i+1];
      OV[m_size-i+1]=itemp;
    }
    break;

  case dp_RandomRow:
    for(i=1; i<=m_size; i++) OV[i]=i;
    if (rseed<=0) rseed=1;
    dp_RandomPermutation(OV,m_size,rseed);
    break;

  case dp_LineShelling:
    for(i=1; i<=m_size; i++) OV[i]=i;
    zvec[0]=1;
    dvec[0]=0;
    if (rseed<=0) rseed=1;
    srand(rseed);
    for(j=2; j<=d_size; j++){
      zvec[j-1]=0;
      dvec[j-1]=d_size-j+1;
      /* dvec[j-1]=rand(); */
    }
    dp_LineShellingOrder(m_size,d_size,A,OV,zvec,dvec);
    break;
  }
}


void dp_SelectPreorderedNext(dp_rowrange m_size,dp_colrange d_size,
    rowset excluded,dp_rowindex OV,dp_rowrange *hnext)
{
  dp_rowrange i,k;
  
  *hnext=0;
  for (i=1; i<=m_size && *hnext==0; i++){
    k=OV[i];
    if (!set_member(k,excluded)) *hnext=k ;
  }
}



void dp_LPSolve(dp_LPPtr lp,dp_ErrorType *err)
/* 
When LP is inconsistent then *re returns the evidence row.
When LP is dual-inconsistent then *se returns the evidence column.
*/
{
  *err=dp_None;
  switch (lp->conv) {
    case dp_LPmax:
      if (lp->solver==dp_CrissCross)
         CrissCrossMaximize(lp->m,lp->d,lp->A,lp->B,lp->objrow,lp->rhscol,lp->use_given_basis,
           &(lp->LPS),&(lp->optvalue),lp->sol,lp->dsol,lp->nbindex,&(lp->re),&(lp->se),&(lp->total_iter),err);
      else
         DualSimplexMaximize(lp->m,lp->d,lp->A,lp->B,lp->objrow,lp->rhscol,lp->use_given_basis,
           &(lp->LPS),&(lp->optvalue),lp->sol,lp->dsol,lp->nbindex,&(lp->re),&(lp->se),&(lp->total_iter),err);
      break;
      
    case dp_LPmin:
      if (lp->solver==dp_CrissCross)
         CrissCrossMinimize(lp->m,lp->d,lp->A,lp->B,lp->objrow,lp->rhscol,lp->use_given_basis,
           &(lp->LPS),&(lp->optvalue),lp->sol,lp->dsol,lp->nbindex,&(lp->re),&(lp->se),&(lp->total_iter),err);
      else
         DualSimplexMinimize(lp->m,lp->d,lp->A,lp->B,lp->objrow,lp->rhscol,lp->use_given_basis,
           &(lp->LPS),&(lp->optvalue),lp->sol,lp->dsol,lp->nbindex,&(lp->re),&(lp->se),&(lp->total_iter),err);
      break;
  }
}

dp_LPPtr dp_MakeLPforInteriorFinding(dp_LPPtr lp)
/* Delete the objective row,
   add an extra column with -1's to the matrix A,
   add an extra row with (bceil, 0,...,0,-1),
   add an objective row with (0,...,0,1), and 
   rows & columns, and change m_size and d_size accordingly, to output new_A.
  This sets up the LP:
  maximize      x_{d+1}
  s.t.    A x + x_{d+1}  <=  b
                x_{d+1}  <=  bm * bmax,
  where bm is set to 2 by default, and bmax=max{1, b[1],...,b[m_size]}.
*/
{
  dp_rowrange m;
  dp_colrange d;
  dp_NumberType NUMB;
  dp_LPConversionType CONV;
  dp_Amatrix A;
  dp_rowrange OBJ; 
  dp_colrange RHS;
  dp_LPPtr lpnew;
  dp_rowrange i; 
  dp_colrange j;
  double bm=2.0,bmax=1,bceil;
  dp_ErrorType err;
  int localdebug=dp_FALSE;

  NUMB=lp->number;
  m=lp->m+1;
  d=lp->d+1;
  OBJ=lp->objrow;
  RHS=lp->rhscol;
  CONV=dp_LPmax;

  dp_InitializeAmatrix(m,d,&A);

  for (i=1; i<=lp->m; i++) {
    if (lp->A[i-1][lp->rhscol-1]>bmax) bmax = lp->A[i-1][lp->rhscol-1];
  }
  bceil=bm*bmax;
  if (localdebug) printf("bceil is set to %g\n",bceil);
  
  for (i=1; i <= m-1; i++) {
    for (j=1; j <= d-1; j++) {
      A[i-1][j-1]=lp->A[i-1][j-1];
      if (localdebug) dp_WriteReal(stdout,A[i-1][j-1]);
    }
    if (localdebug) fprintf(stdout,"\n");
  }
  for (i=1;i<=m; i++) {
    if (i!=OBJ) A[i-1][d-1]=-1.0;  /* new column with all minus one's */
  }
  for (j=1;j<=d-1;j++) {
    A[m-1][j-1]=0.0;    /* new row (bceil, 0,...,0,-1) */
  }
  A[m-1][RHS-1]=bceil;  /* new row (bceil, 0,...,0,-1) */
  for (j=1;j<= d-1;j++) {
    A[OBJ-1][j-1]=0.0;  /* new obj row with (0,...,0,1) */
  }
  A[OBJ-1][d-1]=1.0;    /* new obj row with (0,...,0,1) */

/* A direct but risky way to load an LP.  The next removes (frees) A as well. */
  lpnew=dp_LPDirectLoad(CONV,NUMB,m,d,&A,OBJ,RHS,&err);

/* A safe way to load an LP */
/*  lpnew=dp_LPLoad(CONV,NUMB,m,d,A,OBJ,RHS,&err); */
/*  dp_FreeAmatrix(m,&A); */
 
  return lpnew;
}

void dp_WriteLPResult(FILE *f,dp_LPPtr lp,dp_ErrorType err)
{
  long j;

  fprintf(f,"\n*dplex LP result\n");
  
  if (err!=dp_None) {
    dp_WriteErrorMessages(f,err);
    goto _L99;
  }

  fprintf(f,"* #constraints = %ld\n",lp->m-1);
  fprintf(f,"* #variables   = %ld\n",lp->d-1);

  switch (lp->solver) {
    case dp_DualSimplex:
      fprintf(f,"*Algorithm: dual simplex algorithm\n");break; 
    case dp_CrissCross:
      fprintf(f,"*Algorithm: criss-cross method\n");break;
  }

  switch (lp->conv) {
    case dp_LPmax:
      fprintf(f,"*maximization is chosen\n");break; 
    case dp_LPmin:
      fprintf(f,"*minimization is chosen\n");break;
  }
  
  if (lp->conv==dp_LPmax||lp->conv==dp_LPmin){
    fprintf(f,"*Objective function is\n");  
    for (j=0; j<lp->d; j++){
      if (j>0 && lp->A[lp->objrow-1][j]>=0 ) fprintf(f," +");
      if (j>0 && (j % 5) == 0) fprintf(f,"\n");
      dp_WriteReal(f,lp->A[lp->objrow-1][j]);
      if (j>0) fprintf(f," X[%3ld]",j);
    }
    fprintf(f,"\n");
  }

  switch (lp->LPS){
  case dp_Optimal:
    fprintf(f,"*LP status: a dual pair (x,y) of optimal solutions found.\n");
    fprintf(f,"begin\n");
    fprintf(f,"  primal_solution\n");
    for (j=1; j<lp->d; j++) {
      fprintf(f,"  %3ld : ",j);
      dp_WriteReal(f,lp->sol[j]);
      fprintf(f,"\n");
    }
    fprintf(f,"  dual_solution\n");
    for (j=1; j<lp->d; j++){
      if (lp->nbindex[j+1]>0) {
        fprintf(f,"  %3ld : ",lp->nbindex[j+1]);
        dp_WriteReal(f,lp->dsol[j]); fprintf(f,"\n");
      }
    }
    fprintf(f,"  optimal_value : % .9E\n",lp->optvalue);
    fprintf(f,"end\n");
    break;

  case dp_Inconsistent:
    fprintf(f,"*LP status: LP is inconsistent.\n");
    fprintf(f,"*The positive combination of original inequalities with\n");
    fprintf(f,"*the following coefficients will prove the inconsistency.\n");
    fprintf(f,"begin\n");
    fprintf(f,"  dual_direction\n");
    fprintf(f,"  %3ld : ",lp->re);
    dp_WriteReal(f,1.0);  fprintf(f,"\n");
    for (j=1; j<lp->d; j++){
      if (lp->nbindex[j+1]>0) {
        fprintf(f,"  %3ld : ",lp->nbindex[j+1]);
        dp_WriteReal(f,lp->dsol[j]); fprintf(f,"\n");
      }
    }
    fprintf(f,"end\n");
    break;

  case dp_DualInconsistent: case dp_StrucDualInconsistent:
    fprintf(f,"*LP status: LP is dual inconsistent.\n");
    fprintf(f,"*The linear combination of columns with\n");
    fprintf(f,"*the following coefficients will prove the dual inconsistency.\n");
    fprintf(f,"*(It is also an unbounded direction for the primal LP.)\n");
    fprintf(f,"begin\n");
    fprintf(f,"  primal_direction\n");
    for (j=1; j<lp->d; j++) {
      fprintf(f,"  %3ld : ",j);
      dp_WriteReal(f,lp->sol[j]);
      fprintf(f,"\n");
    }
    fprintf(f,"end\n");
    break;

  default:
    break;
  }
  fprintf(f,"*number of pivot operations = %ld (ph1 = %ld, ph2 = %ld)\n",lp->total_iter,lp->phase1_iter,lp->phase2_iter);
_L99:;
}

void dp_WriteBmatrix(FILE *f,dp_colrange d_size,dp_Bmatrix T)
{
  dp_colrange j1,j2;

  for (j1 = 0; j1 < d_size; j1++) {
    for (j2 = 0; j2 < d_size; j2++) {
      fprintf(f,"%15.7f ",T[j1][j2]);
    }  /*of j2*/
    putc('\n',f);
  }  /*of j1*/
  putc('\n',f);
}

void dp_WriteReal(FILE *f,double x)
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
    fprintf(f," %2ld",ix);
  } else
    fprintf(f," % .9E",x);
}


/* end of dplex.c */


