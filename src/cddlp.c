/* dplex.c:  dual simplex method c-code
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.90, May 28, 2000
*/

/* dplex.c : C-Implementation of the dual simplex method for
   solving an LP: max/min  A_(m-1).x subject to  x in P, where
   P= {x :  A_i.x >= 0, i=0,...,m-2, and  x_0=1}, and
   A_i is the i-th row of an m x n matrix A.
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

#define CDDLPVERSION   "Version 0.90dev (May 12, 2000)"

#define FALSE 0
#define TRUE 1

typedef set_type rowset;  /* set_type defined in setoper.h */
typedef set_type colset;

void CrissCrossMinimize(dd_rowrange,dd_colrange,dd_Amatrix,dd_Bmatrix T,dd_rowset,
  dd_rowrange,dd_colrange,
  dd_LPStatusType *,mytype *optvalue,dd_Arow,dd_Arow,dd_colindex,
  dd_rowrange *,dd_colrange *,long *,dd_ErrorType *);
void CrissCrossMaximize(dd_rowrange,dd_colrange,dd_Amatrix,dd_Bmatrix T,dd_rowset,
  dd_rowrange,dd_colrange,
  dd_LPStatusType *,mytype *optvalue,dd_Arow,dd_Arow,dd_colindex,
  dd_rowrange *,dd_colrange *,long *,dd_ErrorType *);
void DualSimplexMinimize(dd_rowrange,dd_colrange, dd_Amatrix,dd_Bmatrix,dd_rowset,
  dd_rowrange,dd_colrange,
  dd_LPStatusType *,mytype *,dd_Arow,dd_Arow,dd_colindex,
  dd_rowrange *,dd_colrange *,long *,dd_ErrorType *);
void DualSimplexMaximize(dd_rowrange,dd_colrange,dd_Amatrix,dd_Bmatrix,dd_rowset,
  dd_rowrange,dd_colrange,
  dd_LPStatusType *,mytype *,dd_Arow,dd_Arow,dd_colindex,
  dd_rowrange *,dd_colrange *,long *,dd_ErrorType *);
void FindLPBasis(dd_rowrange,dd_colrange,dd_Amatrix,dd_Bmatrix,dd_rowindex,dd_rowset,
    dd_colindex,dd_rowindex,dd_rowrange,dd_colrange,
    dd_colrange *,int *,dd_LPStatusType *,long *);
void FindDualFeasibleBasis(dd_rowrange,dd_colrange,dd_Amatrix,dd_Bmatrix,dd_rowindex,
    dd_colindex,long *,dd_rowrange,dd_colrange,
    dd_colrange *,int *,dd_LPStatusType *,long *);

void dd_WriteBmatrix(FILE *f,dd_colrange d_size,dd_Bmatrix T);
void dd_SetNumberType(char *line,dd_NumberType *number,dd_ErrorType *Error);
void ComputeRowOrderVector2(dd_rowrange m_size,dd_colrange d_size,dd_Amatrix A,
    dd_rowindex OV,dd_RowOrderType ho,unsigned int rseed);
void SelectPreorderedNext2(dd_rowrange m_size,dd_colrange d_size,
    rowset excluded,dd_rowindex OV,dd_rowrange *hnext);
void SetSolutions(dd_rowrange,dd_colrange,
   dd_Amatrix,dd_Bmatrix,dd_rowrange,dd_colrange,dd_LPStatusType,
   mytype *,dd_Arow,dd_Arow,dd_colindex,dd_rowrange,dd_colrange);


dd_LPSolutionPtr dd_LPSolutionLoad(dd_LPPtr lp)
{
  dd_LPSolutionPtr lps;
  dd_colrange j;
  long i;

  lps=(dd_LPSolutionPtr) calloc(1,sizeof(dd_LPSolutionType));
  for (i=1; i<=dd_filenamelen; i++) lps->filename[i-1]=lp->filename[i-1];
  lps->objective=lp->objective;
  lps->solver=lp->solver; 
  lps->m=lp->m;
  lps->d=lp->d;
  lps->numbtype=lp->numbtype;

  lps->LPS=lp->LPS;  /* the current solution status */
  dd_init(lps->optvalue);
  dd_set(lps->optvalue,lp->optvalue);  /* optimal value */
  lps->sol=(mytype*) calloc(lp->d+1,sizeof(mytype));   /* primal solution */
  lps->dsol=(mytype*) calloc(lp->d+1,sizeof(mytype));  /* dual solution */
  lps->nbindex=(long*) calloc((lp->d)+1,sizeof(long));  /* dual solution */
  for (j=0; j<=lp->d; j++){
    dd_init(lps->sol[j]);
    dd_init(lps->dsol[j]); 
    dd_set(lps->sol[j],lp->sol[j]);
    dd_set(lps->dsol[j],lp->dsol[j]);
    lps->nbindex[j]=lp->nbindex[j];
  }
  lps->pivots[0]=lp->pivots[0];
  lps->pivots[1]=lp->pivots[1];
  lps->pivots[2]=lp->pivots[2];
  lps->pivots[3]=lp->pivots[3];
  lps->total_pivots=lp->total_pivots;

  return lps;
}


dd_LPPtr CreateLPData(dd_LPObjectiveType obj,
   dd_NumberType nt,dd_rowrange m,dd_colrange d)
{
  dd_LPType *lp;

  lp=(dd_LPPtr) calloc(1,sizeof(dd_LPType));
  lp->solver=DualSimplex;  /* set the default lp solver */
  lp->d=d;
  lp->m=m;
  lp->numbtype=nt;
  lp->objrow=m;
  lp->rhscol=1L;
  lp->objective=LPnone;
  lp->LPS=LPSundecided;
  lp->eqnumber=0;  /* the number of equalities */

  lp->nbindex=(long*) calloc(d+1,sizeof(long));
  lp->given_nbindex=(long*) calloc(d+1,sizeof(long));
  set_initialize(&(lp->equalityset),m);  
    /* ith component is 1 if it is equality, -1 if it is strict inequality, 0 otherwise. */

  lp->m_alloc=lp->m+2;
  lp->d_alloc=lp->d+2;
  lp->objective=obj;
  dd_InitializeBmatrix(lp->d_alloc,&(lp->B));
  dd_InitializeAmatrix(lp->m_alloc,lp->d_alloc,&(lp->A));
  dd_InitializeArow(lp->d_alloc,&(lp->sol));
  dd_InitializeArow(lp->d_alloc,&(lp->dsol));
  dd_init(lp->optvalue);
  return lp;
}


dd_LPPtr dd_Matrix2LP(dd_MatrixPtr M, dd_ErrorType *err)
{
  dd_rowrange m, i, irev, linc;
  dd_colrange d, j;
  dd_LPType *lp;
  boolean localdebug=FALSE;

  *err=NoError;
  linc=set_card(M->linset);
  m=M->rowsize+1+linc; 
     /* We represent each equation by two inequalities.
        This is not the best way but makes the code simple. */
  d=M->colsize;
  if (localdebug) printf("number of equalities = %ld\n", linc);
  
  lp=CreateLPData(M->objective, M->numbtype, m, d);
  lp->Homogeneous = TRUE;
  lp->eqnumber=linc;  /* this records the number of equations */

  irev=M->rowsize; /* the first row of the linc reversed inequalities. */
  for (i = 1; i <= M->rowsize; i++) {
    if (set_member(i, M->linset)) {
      irev=irev+1;
      set_addelem(lp->equalityset,i);    /* it is equality. */
                                         /* the reversed row irev is not in the equality set. */
      for (j = 1; j <= M->colsize; j++) {
        dd_neg(lp->A[irev-1][j-1],M->matrix[i-1][j-1]);
      }  /*of j*/
      if (localdebug) printf("equality row %ld generates the reverse row %ld.\n",i,irev);
    }
    for (j = 1; j <= M->colsize; j++) {
      dd_set(lp->A[i-1][j-1],M->matrix[i-1][j-1]);
      if (j==1 && i<M->rowsize && dd_Nonzero(M->matrix[i-1][j-1])) lp->Homogeneous = FALSE;
    }  /*of j*/
  }  /*of i*/
  for (j = 1; j <= M->colsize; j++) {
    dd_set(lp->A[m-1][j-1],M->rowvec[j-1]);  /* objective row */
  }  /*of j*/

  return lp;
}


void dd_FreeLPData(dd_LPPtr lp)
{
  if ((lp)!=NULL){
    dd_FreeAmatrix(lp->m_alloc,lp->d_alloc,lp->A);
    dd_FreeBmatrix(lp->d_alloc,lp->B);
    dd_FreeArow(lp->d_alloc,lp->sol);
    dd_FreeArow(lp->d_alloc,lp->dsol);
    set_free(lp->equalityset);
    dd_clear(lp->optvalue);
    free(lp->nbindex);
    free(lp->given_nbindex);
    free(lp);
  }
}

void dd_FreeLPSolution(dd_LPSolutionPtr lps)
{
  if (lps!=NULL){
    free(lps->sol);
    free(lps->dsol);
    free(lps);
  }
}

int dd_LPReverseRow(dd_LPPtr lp, dd_rowrange i)
{
  dd_colrange j;
  int success=0;

  if (i>=1 && i<=lp->m){
    lp->LPS=LPSundecided;
    for (j=1; j<=lp->d; j++) {
      dd_neg(lp->A[i-1][j-1],lp->A[i-1][j-1]);
      /* negating the i-th constraint of A */
    }
    success=1;
  }
  return success;
}

int dd_LPReplaceRow(dd_LPPtr lp, dd_rowrange i, dd_Arow a)
{
  dd_colrange j;
  int success=0;

  if (i>=1 && i<=lp->m){
    lp->LPS=LPSundecided;
    for (j=1; j<=lp->d; j++) {
      dd_set(lp->A[i-1][j-1],a[j-1]);
      /* replacing the i-th constraint by a */
    }
    success=1;
  }
  return success;
}

dd_Arow dd_LPCopyRow(dd_LPPtr lp, dd_rowrange i)
{
  dd_colrange j;
  dd_Arow a;

  if (i>=1 && i<=lp->m){
    dd_InitializeArow(lp->d, &a);
    for (j=1; j<=lp->d; j++) {
      dd_set(a[j-1],lp->A[i-1][j-1]);
      /* copying the i-th row to a */
    }
  }
  return a;
}


void dd_SetNumberType(char *line,dd_NumberType *number,dd_ErrorType *Error)
{
  if (strncmp(line,"integer",7)==0) {
    *number = Integer;
    return;
  }
  else if (strncmp(line,"rational",8)==0) {
    *number = Rational;
    return;
  }
  else if (strncmp(line,"real",4)==0) {
    *number = Real;
    return;
  }
  else { 
    *number=Unknown;
    *Error=ImproperInputFormat;
  }
}


void dd_WriteTableau(FILE *f,dd_rowrange m_size,dd_colrange d_size,dd_Amatrix A,dd_Bmatrix T,
  dd_colindex nbindex,dd_rowindex bflag)
/* Write the tableau  A.T   */
{
  dd_colrange j;
  dd_rowrange i;
  mytype x;
  
  dd_init(x);
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
      TableauEntry(&x,m_size,d_size,A,T,i,j);
      dd_WriteNumber(f,x);
    }
    fprintf(f,"\n");
  }
  fprintf(f,"end\n");
  dd_clear(x);
}


void SelectDualSimplexPivot(dd_rowrange m_size,dd_colrange d_size,
    int Phase1,dd_Amatrix A,dd_Bmatrix T,dd_rowindex OV,
    dd_colindex nbindex,dd_rowindex bflag,
    dd_rowrange objrow,dd_colrange rhscol,
    dd_rowrange *r,dd_colrange *s,int *selected,dd_LPStatusType *lps)
{ 
  /* selects a dual simplex pivot (*r,*s) if the current
     basis is dual feasible and not optimal. If not dual feasible,
     the procedure returns *selected=FALSE and *lps=LPSundecided.
     If Phase1=TRUE, the RHS column will be considered as the negative
     of the column of the largest variable (==m_size).  For this case, it is assumed
     that the caller used the auxiliary row (with variable m_size) to make the current
     dictionary dual feasible before calling this routine so that the nonbasic
     column for m_size corresponds to the auxiliary variable.
  */
  int colselected=FALSE,rowselected=FALSE,
    dualfeasible=TRUE,localdebug=FALSE;
  dd_rowrange i;
  dd_colrange j;
  mytype val,minval,rat,minrat;
  static dd_Arow rcost;
  static dd_colrange d_last=0;

  dd_init(val); dd_init(minval); dd_init(rat); dd_init(minrat);
  if (d_last<d_size) {
    if (d_last>0) {
      for (j=1; j<=d_last; j++){ dd_clear(rcost[j-1]);}
      free(rcost);
    }
    rcost=(mytype*) calloc(d_size,sizeof(mytype));
    for (j=1; j<=d_size; j++){ dd_init(rcost[j-1]);}
  }
  d_last=d_size;

  *r=0; *s=0;
  *selected=FALSE;
  *lps=LPSundecided;
  for (j=1; j<=d_size; j++){
    if (j!=rhscol){
      TableauEntry(&(rcost[j-1]),m_size,d_size,A,T,objrow,j);
      if (dd_Positive(rcost[j-1])) { 
        dualfeasible=FALSE;
      }
    }
  }
  if (dualfeasible){
    while ((*lps==LPSundecided) && (!rowselected) && (!colselected)) {
      for (i=1; i<=m_size; i++) {
        if (i!=objrow && bflag[i]==-1) {  /* i is a basic variable */
          if (Phase1){
            TableauEntry(&val, m_size,d_size,A,T,i,bflag[m_size]);
            dd_neg(val,val);
            /* for dual Phase I */
          } 
          else {TableauEntry(&val,m_size,d_size,A,T,i,rhscol);}
          if (dd_Smaller(val,minval)) {
            *r=i;
            dd_set(minval,val);
          }
        }
      }
      if (dd_Nonnegative(minval)) {
        *lps=Optimal;
      }
      else {
        rowselected=TRUE;
        for (j=1; j<=d_size; j++){
          TableauEntry(&val,m_size,d_size,A,T,*r,j);
          if (j!=rhscol && dd_Positive(val)) {
            dd_div(rat,rcost[j-1],val);
            dd_neg(rat,rat);
            if (*s==0 || dd_Smaller(rat,minrat)){
              dd_set(minrat,rat);
              *s=j;
            }
          }
        }
        if (*s>0) {colselected=TRUE; *selected=TRUE;}
        else *lps=Inconsistent;
      }
    } /* end of while */
  }
  if (localdebug) {
     if (Phase1) printf("Phase 1 : select %ld,%ld\n",*r,*s);
     else printf("Phase 2 : select %ld,%ld\n",*r,*s);
  }
  dd_clear(val); dd_clear(minval); dd_clear(rat); dd_clear(minrat);
}

void TableauEntry(mytype *x,dd_rowrange m_size, dd_colrange d_size, dd_Amatrix X, dd_Bmatrix T,
				dd_rowrange r, dd_colrange s)
/* Compute the (r,s) entry of X.T   */
{
  dd_colrange j;
  mytype temp;

  dd_init(temp);
  dd_set(*x,dd_purezero);
  for (j=0; j< d_size; j++) {
    dd_mul(temp,X[r-1][j], T[j][s-1]);
    dd_add(*x, *x, temp);
  }
  dd_clear(temp);
}

void SelectPivot2(dd_rowrange m_size,dd_colrange d_size,dd_Amatrix A,dd_Bmatrix T,
            dd_RowOrderType roworder,dd_rowindex ordervec, rowset equalityset,
            dd_rowrange rowmax,rowset NopivotRow,
            colset NopivotCol,dd_rowrange *r,dd_colrange *s,
            boolean *selected)
/* Select a position (*r,*s) in the matrix A.T such that (A.T)[*r][*s] is nonzero
   The choice is feasible, i.e., not on NopivotRow and NopivotCol, and
   best with respect to the specified roworder 
 */
{
  int stop;
  dd_rowrange i,rtemp;
  rowset rowexcluded;
  mytype Xtemp;
  boolean localdebug=FALSE;

  stop = FALSE;
  localdebug=debug;
  dd_init(Xtemp);
  set_initialize(&rowexcluded,m_size);
  set_copy(rowexcluded,NopivotRow);
  for (i=rowmax+1;i<=m_size;i++) {
    set_addelem(rowexcluded,i);   /* cannot pivot on any row > rmax */
  }
  *selected = FALSE;
  do {
    rtemp=0; i=1;
    while (i<=m_size && rtemp==0) {  /* equalityset vars have highest priorities */
      if (set_member(i,equalityset) && !set_member(i,rowexcluded)){
        if (localdebug) printf("marked set %ld chosen as a candidate\n",i);
        rtemp=i;
      }
      i++;
    }
    if (rtemp==0) SelectPreorderedNext2(m_size,d_size,rowexcluded,ordervec,&rtemp);;
    if (rtemp>=1) {
      *r=rtemp;
      *s=1;
      while (*s <= d_size && !*selected) {
        TableauEntry(&Xtemp,m_size,d_size,A,T,*r,*s);
        if (!set_member(*s,NopivotCol) && dd_Nonzero(Xtemp)) {
          *selected = TRUE;
          stop = TRUE;
        } else {
          (*s)++;
        }
      }
      if (!*selected) {
        set_addelem(rowexcluded,rtemp);
      }
    }
    else {
      *r = 0;
      *s = 0;
      stop = TRUE;
    }
  } while (!stop);
  set_free(rowexcluded); dd_clear(Xtemp);
}

void GaussianColumnPivot(dd_rowrange m_size, dd_colrange d_size, 
    dd_Amatrix X, dd_Bmatrix T, dd_rowrange r, dd_colrange s)
/* Update the Transformation matrix T with the pivot operation on (r,s) 
   This procedure performs a implicit pivot operation on the matrix X by
   updating the dual basis inverse  T.
 */
{
  dd_colrange j, j1;
  mytype Xtemp0, Xtemp1, Xtemp;
  static dd_Arow Rtemp;
  static dd_colrange last_d=0;
  boolean localdebug=debug;

  dd_init(Xtemp0); dd_init(Xtemp1); dd_init(Xtemp);
  if (last_d!=d_size){
    if (last_d>0) {
      for (j=1; j<=last_d; j++) dd_clear(Rtemp[j-1]);
      free(Rtemp);
    }
    Rtemp=(mytype*)calloc(d_size,sizeof(mytype));
    for (j=1; j<=d_size; j++) dd_init(Rtemp[j-1]);
    last_d=d_size;
  }

  for (j=1; j<=d_size; j++) {
    TableauEntry(&(Rtemp[j-1]), m_size, d_size, X, T, r,j);
  }
  dd_set(Xtemp0,Rtemp[s-1]);
  if (localdebug) {
    printf("Gaussian Pivot: pivot entry = "); dd_WriteNumber(stdout,Xtemp0);
    printf("\n");
  }
  for (j = 1; j <= d_size; j++) {
    if (j != s) {
      dd_div(Xtemp,Rtemp[j-1],Xtemp0);
      dd_set(Xtemp1,dd_purezero);
      for (j1 = 1; j1 <= d_size; j1++){
        dd_mul(Xtemp1,Xtemp,T[j1-1][s - 1]);
        dd_sub(T[j1-1][j-1],T[j1-1][j-1],Xtemp1);
 /*     T[j1-1][j-1] -= T[j1-1][s - 1] * Xtemp / Xtemp0;  */
        if (localdebug){
          dd_WriteNumber(stdout, T[j1-1][j-1]);
        }
      }
      if (localdebug) printf("\n");
    }
    if (localdebug) printf("\n");
  }
  for (j = 1; j <= d_size; j++)
    dd_div(T[j-1][s - 1],T[j-1][s - 1],Xtemp0);

  dd_clear(Xtemp0); dd_clear(Xtemp1); dd_clear(Xtemp);
}

void GaussianColumnPivot2(dd_rowrange m_size,dd_colrange d_size,
    dd_Amatrix A,dd_Bmatrix T,dd_colindex nbindex,dd_rowindex bflag,dd_rowrange r,dd_colrange s)
/* Update the Transformation matrix T with the pivot operation on (r,s) 
   This procedure performs a implicit pivot operation on the matrix A by
   updating the dual basis inverse  T.
 */
{
  int localdebug=debug;
  long entering;

  GaussianColumnPivot(m_size,d_size,A,T,r,s);
  entering=nbindex[s];
  bflag[r]=s;     /* the nonbasic variable r corresponds to column s */
  nbindex[s]=r;   /* the nonbasic variable on s column is r */
  if (localdebug) {
    fprintf(stdout,"Column pivot: (leaving, entering) = (%ld, %ld)\n", r,entering);
    if (m_size <=10){
      dd_WriteBmatrix(stdout,d_size,T);
      dd_WriteTableau(stdout,m_size,d_size,A,T,nbindex,bflag);
    }
  }

  if (entering>0) bflag[entering]=-1;
     /* original variables have negative index and should not affect the row index */
}


void ResetTableau(dd_rowrange m_size,dd_colrange d_size,dd_Bmatrix T,
    dd_colindex nbindex,dd_rowindex bflag,dd_rowrange objrow,dd_colrange rhscol)
{
  dd_rowrange i;
  dd_colrange j;
  
  /* Initialize T and nbindex */
  for (j=1; j<=d_size; j++) nbindex[j]=-j;
  nbindex[rhscol]=0; 
    /* RHS is already in nonbasis and is considered to be associated
       with the zero-th row of input. */
  dd_SetToIdentity(d_size,T);
  
  /* Set the bflag according to nbindex */
  for (i=1; i<=m_size; i++) bflag[i]=-1;  
    /* all basic variables have index -1 */
  bflag[objrow]= 0; 
    /* bflag of the objective variable is 0,
       different from other basic variables which have -1 */
  for (j=1; j<=d_size; j++) if (nbindex[j]>0) bflag[nbindex[j]]=j;
    /* bflag of a nonbasic variable is its column number */

}

void dd_SelectCrissCrossPivot(dd_rowrange m_size,dd_colrange d_size,dd_Amatrix A,dd_Bmatrix T,
    dd_rowindex bflag,dd_rowrange objrow,dd_colrange rhscol,
    dd_rowrange *r,dd_colrange *s,
    int *selected,dd_LPStatusType *lps)
{
  int colselected=FALSE,rowselected=FALSE;
  dd_rowrange i;
  mytype val;
  
  dd_init(val);
  *selected=FALSE;
  *lps=LPSundecided;
  while ((*lps==LPSundecided) && (!rowselected) && (!colselected)) {
    for (i=1; i<=m_size; i++) {
      if (i!=objrow && bflag[i]==-1) {  /* i is a basic variable */
        TableauEntry(&val,m_size,d_size,A,T,i,rhscol);
        if (dd_Negative(val)) {
          rowselected=TRUE;
          *r=i;
          break;
        }
      }
      else if (bflag[i] >0) { /* i is nonbasic variable */
        TableauEntry(&val,m_size,d_size,A,T,objrow,bflag[i]);
        if (dd_Positive(val)) {
          colselected=TRUE;
          *s=bflag[i];
          break;
        }
      }
    }
    if  ((!rowselected) && (!colselected)) {
      *lps=Optimal;
      return;
    }
    else if (rowselected) {
     for (i=1; i<=m_size; i++) {
       if (bflag[i] >0) { /* i is nonbasic variable */
          TableauEntry(&val,m_size,d_size,A,T,*r,bflag[i]);
          if (dd_Positive(val)) {
            colselected=TRUE;
            *s=bflag[i];
            *selected=TRUE;
            break;
          }
        }
      }
    }
    else if (colselected) {
      for (i=1; i<=m_size; i++) {
        if (i!=objrow && bflag[i]==-1) {  /* i is a basic variable */
          TableauEntry(&val,m_size,d_size,A,T,i,*s);
          if (dd_Negative(val)) {
            rowselected=TRUE;
            *r=i;
            *selected=TRUE;
            break;
          }
        }
      }
    }
    if (!rowselected) {
      *lps=DualInconsistent;
    }
    else if (!colselected) {
      *lps=Inconsistent;
    }
  }
  dd_clear(val);
}

void CrissCrossMinimize(dd_rowrange m_size,dd_colrange d_size,
    dd_Amatrix A,dd_Bmatrix T,dd_rowset equalityset,
    dd_rowrange objrow,dd_colrange rhscol,dd_LPStatusType *LPS,
    mytype *optvalue,dd_Arow sol,dd_Arow dsol,dd_colindex nbindex,
    dd_rowrange *re,dd_colrange *se,long *pivots,dd_ErrorType *err)
{
   dd_colrange j;

   *err=NoError;
   for (j=1; j<=d_size; j++)
     dd_neg(A[objrow-1][j-1],A[objrow-1][j-1]);
   CrissCrossMaximize(m_size,d_size,A,T,equalityset, objrow,rhscol,
     LPS,optvalue,sol,dsol,nbindex,re,se,pivots,err);
   dd_neg(*optvalue,*optvalue);
   for (j=1; j<=d_size; j++){
     dd_neg(dsol[j-1],dsol[j-1]);
     dd_neg(A[objrow-1][j-1],A[objrow-1][j-1]);
   }
}

void CrissCrossMaximize(dd_rowrange m_size,dd_colrange d_size,
    dd_Amatrix A,dd_Bmatrix T,dd_rowset equalityset,
    dd_rowrange objrow,dd_colrange rhscol,dd_LPStatusType *LPS,
    mytype *optvalue,dd_Arow sol,dd_Arow dsol,dd_colindex nbindex,
    dd_rowrange *re,dd_colrange *se,long *pivots,dd_ErrorType *err)
/* 
When LP is inconsistent then *re returns the evidence row.
When LP is dual-inconsistent then *se returns the evidence column.
*/
{
  int stop,chosen,found;
  long pivots0,pivots1;
  dd_rowrange i,r;
  dd_colrange s;
  static dd_rowindex bflag;
  static long mlast=0;
  static dd_rowindex OrderVector;  /* the permutation vector to store a preordered row indeces */
  unsigned int rseed=1;

  *err=NoError;
  for (i=0; i<= 3; i++) pivots[i]=0;
  if (bflag==NULL || mlast!=m_size){
     if (mlast!=m_size && mlast>0) {
       free(bflag);   /* called previously with different m_size */
       free(OrderVector);
     }
     bflag=(long *) calloc(m_size+1,sizeof(long*));
     OrderVector=(long *)calloc(m_size+1,sizeof(long*)); 
     /* initialize only for the first time or when a larger space is needed */
     mlast=m_size;
  }
  /* Initializing control variables. */
  ComputeRowOrderVector2(m_size,d_size,A,OrderVector,MinIndex,rseed);

  *re=0; *se=0; pivots1=0;

  ResetTableau(m_size,d_size,T,nbindex,bflag,objrow,rhscol);

  FindLPBasis(m_size,d_size,A,T,OrderVector, equalityset,
      nbindex,bflag,
      objrow,rhscol,&s,&found,LPS,&pivots0);
  pivots[0]=pivots0;

  if (!found){
     *se=s;
     goto _L99;
     /* No LP basis is found, and thus Inconsistent.  
     Output the evidence column. */
  }

  stop=FALSE;
  do {   /* Criss-Cross Method */
    dd_SelectCrissCrossPivot(m_size,d_size,A,T,bflag,
       objrow,rhscol,&r,&s,&chosen,LPS);
    if (chosen) {
      GaussianColumnPivot2(m_size,d_size,A,T,nbindex,bflag,r,s);
      pivots1++;
    } else {
      switch (*LPS){
        case Inconsistent: *re=r;
        case DualInconsistent: *se=s;
        default: break;
      }
      stop=TRUE;
    }
  } while(!stop);
  pivots[1]=pivots1;
  
_L99:

  SetSolutions(m_size,d_size,A,T,
   objrow,rhscol,*LPS,optvalue,sol,dsol,nbindex,*re,*se);

}


int dd_LexSmaller(mytype *v1,mytype *v2,long dmax)
{ /* dmax is the size of vectors v1,v2 */
  int determined,smaller;
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

int dd_LexLarger(mytype *v1,mytype *v2,long dmax)
{
  return dd_LexSmaller(v2,v1,dmax);
}

void FindLPBasis(dd_rowrange m_size,dd_colrange d_size,
    dd_Amatrix A, dd_Bmatrix T,dd_rowindex OV,dd_rowset equalityset, dd_colindex nbindex,
    dd_rowindex bflag,dd_rowrange objrow,dd_colrange rhscol,
    dd_colrange *cs,int *found,dd_LPStatusType *lps,long *pivot_no)
{ 
  /* Find a LP basis using Gaussian pivots.
     If the problem has an LP basis,
     the procedure returns *found=TRUE,*lps=LPSundecided and an LP basis.
     If the constraint matrix A (excluding the rhs and objective) is not
     column indepent, there are two cases.  If the dependency gives a dual
     inconsistency, this returns *found=FALSE, *lps=dd_StrucDualInconsistent and 
     the evidence column *s.  Otherwise, this returns *found=TRUE, 
     *lps=LPSundecided and an LP basis of size less than d_size.  Columns j
     that do not belong to the basis (i.e. cannot be chosen as pivot because
     they are all zero) will be indicated in nbindex vector: nbindex[j] will
     be negative and set to -j.
  */
  int chosen,stop;
  long pivots_p0=0,rank;
  colset ColSelected;
  rowset RowSelected;
  mytype val;

  dd_rowrange r;
  dd_colrange j,s;

  dd_init(val);
  *found=FALSE; *cs=0; rank=0;
  *lps=LPSundecided;

  set_initialize(&RowSelected,m_size);
  set_initialize(&ColSelected,d_size);
  set_addelem(RowSelected,objrow);
  set_addelem(ColSelected,rhscol);

  stop=FALSE;
  do {   /* Find a LP basis */
    SelectPivot2(m_size,d_size,A,T,MinIndex,OV,equalityset,
      m_size,RowSelected,ColSelected,&r,&s,&chosen);
    if (chosen) {
      set_addelem(RowSelected,r);
      set_addelem(ColSelected,s);
      GaussianColumnPivot2(m_size,d_size,A,T,nbindex,bflag,r,s);
      pivots_p0++;
      rank++;
    } else {
      for (j=1;j<=d_size  && *lps==LPSundecided; j++) {
        if (j!=rhscol && nbindex[j]<0){
          TableauEntry(&val,m_size,d_size,A,T,objrow,j);
          if (dd_Nonzero(val)){  /* dual inconsistent */
            *lps=StrucDualInconsistent;
            *cs=j;
            /* dual inconsistent because the nonzero reduced cost */
          }
        }
      }
      if (*lps==LPSundecided) *found=TRUE;  
         /* dependent columns but not dual inconsistent. */
      stop=TRUE;
    }
    if (rank==d_size-1) {
      stop = TRUE;
      *found=TRUE;
    }
  } while (!stop);

  *pivot_no=pivots_p0;
  set_free(RowSelected);
  set_free(ColSelected);
  dd_clear(val);
}

void FindDualFeasibleBasis(dd_rowrange m_size,dd_colrange d_size,
    dd_Amatrix A,dd_Bmatrix T,dd_rowindex OV,
    dd_colindex nbindex,dd_rowindex bflag,dd_rowrange objrow,dd_colrange rhscol,
    dd_colrange *s,int *found,dd_LPStatusType *lps,long *pivot_no)
{ 
  /* Find a dual feasible basis using Phase I of Dual Simplex method.
     If the problem is dual feasible,
     the procedure returns *found=TRUE, *lps=LPSundecided and a dual feasible
     basis.   If the problem is dual infeasible, this returns
     *found=FALSE, *lps=DualInconsistent and the evidence column *s.
     Caution: matrix A must have at least one extra row:  the row space A[m_size] must
     have been allocated.
  */
  int phase1,dualfeasible=TRUE,localdebug=FALSE,chosen,stop;
  dd_LPStatusType LPSphase1;
  long pivots_p1=0;
  dd_rowrange i,r_val;
  dd_colrange j,l,ms=0,s_val,local_m_size;
  mytype x,val,maxcost;
  static dd_colrange d_last=0;
  static dd_Arow rcost;

  dd_init(x); dd_init(val);
  dd_init(maxcost);  dd_set(maxcost,dd_minuszero);
  if (d_last<d_size) {
    if (d_last>0) {
      for (j=1; j<=d_last; j++){ dd_clear(rcost[j-1]);}
      free(rcost);
    }
    rcost=(mytype*) calloc(d_size,sizeof(mytype));
    for (j=1; j<=d_size; j++){ dd_init(rcost[j-1]);}
  }
  d_last=d_size;

  *found=TRUE; *lps=LPSundecided; *s=0;
  local_m_size=m_size+1;  /* increase m_size by 1 */

  ms=0;  /* ms will be the index of column which has the largest reduced cost */
  for (j=1; j<=d_size; j++){
    if (j!=rhscol){
      if (localdebug) printf("checking the column %ld var %ld\n",j,nbindex[j]); 
      TableauEntry(&(rcost[j-1]),local_m_size,d_size,A,T,objrow,j);
      if (localdebug) {printf("reduced cost = "); dd_WriteNumber(stdout, rcost[j-1]); }
      if (dd_Larger(rcost[j-1],maxcost)) {dd_set(maxcost,rcost[j-1]); ms = j;}
    }
  }
  if (dd_Positive(maxcost)) dualfeasible=FALSE;

  if (!dualfeasible){
    for (j=1; j<=d_size; j++){
      dd_set(A[local_m_size-1][j-1], dd_purezero);
      for (l=1; l<=d_size; l++){
        if (nbindex[l]>0) {
          dd_sub(A[local_m_size-1][j-1],A[local_m_size-1][j-1],A[nbindex[l]-1][j-1]); 
          /* To make the auxiliary row (0,-1,-1,...,-1).  */
        }
      }
    }
    if (localdebug){
      printf("FindDualFeasibleBasis: curruent basis is not dual feasible.\n");
      printf("because of the column %ld assoc. with var %ld   dual cost =",
       ms,nbindex[ms]);
      dd_WriteNumber(stdout, maxcost);
    }

    /* Pivot on (local_m_size,ms) so that the dual basic solution becomes feasible */
    GaussianColumnPivot2(local_m_size,d_size,A,T,nbindex,bflag,local_m_size,ms);
    pivots_p1=pivots_p1+1;

    phase1=TRUE; stop=FALSE;
    do {   /* Dual Simplex Phase I */
      chosen=FALSE; LPSphase1=LPSundecided;
      SelectDualSimplexPivot(local_m_size,d_size,phase1,A,T,OV,nbindex,bflag,
        objrow,rhscol,&r_val,&s_val,&chosen,&LPSphase1);
      if (!chosen) {
        /* The current dictionary is terminal.  There are two cases:
           TableauEntry(local_m_size,d_size,A,T,objrow,ms) is negative or zero.
           The first case implies dual infeasible,
           and the latter implies dual feasible but local_m_size is still in nonbasis.
           We must pivot in the auxiliary variable local_m_size. 
        */

        mytype minval;
        dd_init(minval);
        r_val=0;
        for (i=1; i<=local_m_size; i++){
          if (bflag[i]<0) { 
             /* i is basic and not the objective variable */
            TableauEntry(&val,local_m_size,d_size,A,T,i,ms);  /* auxiliary column*/
            if (dd_Smaller(val, minval)) {
              r_val=i;
              dd_set(minval,val);
              if (localdebug) {
                printf("update minval with = ");
                dd_WriteNumber(stdout, minval);
                printf("  r_val = %ld\n",r_val);
              }
            }
          }
        }
        dd_clear(minval);

        GaussianColumnPivot2(local_m_size,d_size,A,T,nbindex,bflag,r_val,ms);
        pivots_p1=pivots_p1+1;

        TableauEntry(&x,local_m_size,d_size,A,T,objrow,ms);
        if (dd_Negative(x)){
          *found=FALSE; *lps=DualInconsistent;  *s=ms;
        }
        stop=TRUE;
      } else {
        GaussianColumnPivot2(local_m_size,d_size,A,T,nbindex,bflag,r_val,s_val);
        pivots_p1=pivots_p1+1;
        if (bflag[local_m_size]<0) {
          stop=TRUE; 
          if (localdebug) 
            printf("Dual Phase I: the auxiliary variable entered the basis, go to phase II\n");
        }
      }
    } while(!stop);
  }
  *pivot_no=pivots_p1;
  dd_clear(x); dd_clear(val); dd_clear(maxcost);
}

void DualSimplexMinimize(dd_rowrange m_size,dd_colrange d_size,
   dd_Amatrix A,dd_Bmatrix T,dd_rowset equalityset,
   dd_rowrange objrow,dd_colrange rhscol,dd_LPStatusType *LPS,
   mytype *optvalue,dd_Arow sol,dd_Arow dsol,dd_colindex nbindex,
   dd_rowrange *re,dd_colrange *se,long *pivots,dd_ErrorType *err)
{
   dd_colrange j;

   *err=NoError;
   for (j=1; j<=d_size; j++)
     dd_neg(A[objrow-1][j-1],A[objrow-1][j-1]);
   DualSimplexMaximize(m_size,d_size,A,T,equalityset, objrow,rhscol,
     LPS,optvalue,sol,dsol,nbindex,re,se,pivots,err);
   dd_neg(*optvalue,*optvalue);
   for (j=1; j<=d_size; j++){
     dd_neg(dsol[j-1],dsol[j-1]);
     dd_neg(A[objrow-1][j-1],A[objrow-1][j-1]);
   }
}

void DualSimplexMaximize(dd_rowrange m_size,dd_colrange d_size,
   dd_Amatrix A,dd_Bmatrix T,dd_rowset equalityset,
   dd_rowrange objrow,dd_colrange rhscol,dd_LPStatusType *LPS,
   mytype *optvalue,dd_Arow sol,dd_Arow dsol,dd_colindex nbindex,
   dd_rowrange *re,dd_colrange *se,long *pivots,dd_ErrorType *err)
/* 
When LP is inconsistent then *re returns the evidence row.
When LP is dual-inconsistent then *se returns the evidence column.
*/
{
  int stop,chosen,phase1,found;
  long pivots_ds=0,pivots_p0=0,pivots_p1=0,pivots_pc=0,maxpivots,maxpivfactor=70;
  dd_rowrange i,r;
  dd_colrange s;
  static dd_rowindex bflag;
  static long mlast=0,nlast=0;
  static dd_rowindex OrderVector;  /* the permutation vector to store a preordered row indeces */
  unsigned int rseed=1;
  
  *err=NoError;
  for (i=0; i<= 3; i++) pivots[i]=0;
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
  ComputeRowOrderVector2(m_size,d_size,A,OrderVector,MinIndex,rseed);

  *re=0; *se=0;
  
  ResetTableau(m_size,d_size,T,nbindex,bflag,objrow,rhscol);
   
  FindLPBasis(m_size,d_size,A,T,OrderVector,equalityset,nbindex,bflag,
      objrow,rhscol,&s,&found,LPS,&pivots_p0);
  pivots[0]=pivots_p0;

  if (!found){
     *se=s;
     goto _L99;
     /* No LP basis is found, and thus Inconsistent.  
     Output the evidence column. */
  }

  FindDualFeasibleBasis(m_size,d_size,A,T,OrderVector,nbindex,bflag,
      objrow,rhscol,&s,&found,LPS,&pivots_p1);
  pivots[1]=pivots_p1;

  if (!found){
     *se=s;
     goto _L99;
     /* No dual feasible basis is found, and thus DualInconsistent.  
     Output the evidence column. */
  }
  
  /* Dual Simplex Method */
  stop=FALSE;
  do {
    chosen=FALSE; *LPS=LPSundecided; phase1=FALSE;
    if (pivots_ds<maxpivots) {
      SelectDualSimplexPivot(m_size,d_size,phase1,A,T,OrderVector,nbindex,bflag,
        objrow,rhscol,&r,&s,&chosen,LPS);
    }
    if (chosen) pivots_ds=pivots_ds+1;
    if (!chosen && *LPS==LPSundecided) {  
      /* In principle this should not be executed because we already have dual feasibility
         attained and dual simplex pivot should have been chosen.  This might occur
         under floating point computation, or the case of cycling.
      */
      dd_SelectCrissCrossPivot(m_size,d_size,A,T,bflag,
        objrow,rhscol,&r,&s,&chosen,LPS);
      if (chosen) pivots_pc=pivots_pc+1;
    }
    if (chosen) {
      GaussianColumnPivot2(m_size,d_size,A,T,nbindex,bflag,r,s);
    } else {
      switch (*LPS){
        case Inconsistent: *re=r;
        case DualInconsistent: *se=s;
        default: break;
      }
      stop=TRUE;
    }
  } while(!stop);
  pivots[2]=pivots_ds;
  pivots[3]=pivots_pc;

_L99:
  
  SetSolutions(m_size,d_size,A,T,
   objrow,rhscol,*LPS,optvalue,sol,dsol,nbindex,*re,*se);

}

void SetSolutions(dd_rowrange m_size,dd_colrange d_size,
   dd_Amatrix A,dd_Bmatrix T,
   dd_rowrange objrow,dd_colrange rhscol,dd_LPStatusType LPS,
   mytype *optvalue,dd_Arow sol,dd_Arow dsol,dd_colindex nbindex,
   dd_rowrange re,dd_colrange se)
/* 
Assign the solution vectors to sol,dsol,*optvalue after solving
the LP.
*/
{
  dd_colrange j;
  mytype x,sw;
  int localdebug=FALSE;
  
  dd_init(x); dd_init(sw);
  switch (LPS){
  case Optimal:
    for (j=1;j<=d_size; j++) {
      dd_set(sol[j-1],T[j-1][rhscol-1]);
      TableauEntry(&x,m_size,d_size,A,T,objrow,j);
      dd_neg(dsol[j-1],x);
      TableauEntry(optvalue,m_size,d_size,A,T,objrow,rhscol);
      if (localdebug) {printf("dsol[%ld]= ",nbindex[j]); dd_WriteNumber(stdout, dsol[j-1]); }
    }
    break;
  case Inconsistent:
    if (localdebug) printf("DualSimplexSolve: LP is inconsistent.\n");
    for (j=1;j<=d_size; j++) {
      dd_set(sol[j-1],T[j-1][rhscol-1]);
      TableauEntry(&x,m_size,d_size,A,T,re,j);
      dd_neg(dsol[j-1],x);
      if (localdebug) {printf("dsol[%ld]= ",nbindex[j]); dd_WriteNumber(stdout,dsol[j-1]);}
    }
    break;
  case DualInconsistent:
    for (j=1;j<=d_size; j++) {
      dd_set(sol[j-1],T[j-1][se-1]);
      TableauEntry(&x,m_size,d_size,A,T,objrow,j);
      dd_neg(dsol[j-1],x);
      if (localdebug) {printf("dsol[%ld]= \n",nbindex[j]);dd_WriteNumber(stdout,dsol[j-1]);}
    }
    if (localdebug) printf( "DualSimplexSolve: LP is dual inconsistent.\n");
    break;

  case StrucDualInconsistent:
    TableauEntry(&x,m_size,d_size,A,T,objrow,se);
    if (dd_Positive(x)) dd_set(sw,dd_one);
    else dd_neg(sw,dd_one);
    for (j=1;j<=d_size; j++) {
      dd_mul(sol[j-1],sw,T[j-1][se-1]);
      TableauEntry(&x,m_size,d_size,A,T,objrow,j);
      dd_neg(dsol[j-1],x);
      if (localdebug) {printf("dsol[%ld]= ",nbindex[j]);dd_WriteNumber(stdout,dsol[j-1]);}
    }
    if (localdebug) printf( "DualSimplexSolve: LP is dual inconsistent.\n");
    break;

  default:break;
  }
  dd_clear(x); dd_clear(sw);
}


long dd_Partition(dd_rowindex OV,long p,long r,dd_Amatrix A,long dmax)
{
  mytype *x;
  long i,j,ovi;
  
  x=A[OV[p]-1];
  printf("x = "); for (j=1; j<=dmax; j++) dd_WriteNumber(stdout, x[j-1]);
  printf("\n");

  i=p-1;
  j=r+1;
  while (TRUE){
    do{
      j--;
    } while (dd_LexLarger(A[OV[j]-1],x,dmax));
    do{
      i++;
    } while (dd_LexSmaller(A[OV[i]-1],x,dmax));
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

void QuickSort2(dd_rowindex OV,long p,long r,dd_Amatrix A,long dmax)
{
  long q;
  
  if (p < r){
    q = dd_Partition(OV,p,r,A,dmax);
    QuickSort2(OV,p,q,A,dmax);
    QuickSort2(OV,q+1,r,A,dmax);
  }
}


#ifndef RAND_MAX 
#define RAND_MAX 32767 
#endif

void RandomPermutation2(dd_rowindex OV,long t,unsigned int seed)
{
  long k,j,ovj;
  double u,xk,r,rand_max=(double) RAND_MAX;
  int localdebug=FALSE;

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

void ComputeRowOrderVector2(dd_rowrange m_size,dd_colrange d_size,dd_Amatrix A,
    dd_rowindex OV,dd_RowOrderType ho,unsigned int rseed)
{
  long i,itemp;
  
  OV[0]=0;
  switch (ho){
  case MaxIndex:
    for(i=1; i<=m_size; i++) OV[i]=m_size-i+1;
    break;

  case LexMin:
    for(i=1; i<=m_size; i++) OV[i]=i;
    QuickSort2(OV,1,m_size,A,d_size);
   break;

  case LexMax:
    for(i=1; i<=m_size; i++) OV[i]=i;
    QuickSort2(OV,1,m_size,A,d_size);
    for(i=1; i<=m_size/2;i++){   /* just reverse the order */
      itemp=OV[i];
      OV[i]=OV[m_size-i+1];
      OV[m_size-i+1]=itemp;
    }
    break;

  case RandomRow:
    for(i=1; i<=m_size; i++) OV[i]=i;
    if (rseed<=0) rseed=1;
    RandomPermutation2(OV,m_size,rseed);
    break;

  case MinIndex: 
    for(i=1; i<=m_size; i++) OV[i]=i;
    break;

  default: 
    for(i=1; i<=m_size; i++) OV[i]=i;
    break;
 }
}

void SelectPreorderedNext2(dd_rowrange m_size,dd_colrange d_size,
    rowset excluded,dd_rowindex OV,dd_rowrange *hnext)
{
  dd_rowrange i,k;
  
  *hnext=0;
  for (i=1; i<=m_size && *hnext==0; i++){
    k=OV[i];
    if (!set_member(k,excluded)) *hnext=k ;
  }
}



boolean dd_LPSolve(dd_LPPtr lp,dd_LPSolverType solver,dd_ErrorType *err)
/* 
When LP is inconsistent then *re returns the evidence row.
When LP is dual-inconsistent then *se returns the evidence column.
*/
{
  int i;
  boolean found=FALSE;

  *err=NoError;
  lp->solver=solver;
  time(&lp->starttime);
  switch (lp->objective) {
    case LPmax:
      if (solver==CrissCross)
         CrissCrossMaximize(lp->m,lp->d,lp->A,lp->B,lp->equalityset,lp->objrow,lp->rhscol,
           &(lp->LPS),&(lp->optvalue),lp->sol,lp->dsol,lp->nbindex,&(lp->re),&(lp->se),lp->pivots,err);
      else
         DualSimplexMaximize(lp->m,lp->d,lp->A,lp->B,lp->equalityset,lp->objrow,lp->rhscol,
           &(lp->LPS),&(lp->optvalue),lp->sol,lp->dsol,lp->nbindex,&(lp->re),&(lp->se),lp->pivots,err);
      break;
      
    case LPmin:
      if (solver==CrissCross)
         CrissCrossMinimize(lp->m,lp->d,lp->A,lp->B,lp->equalityset,lp->objrow,lp->rhscol,
           &(lp->LPS),&(lp->optvalue),lp->sol,lp->dsol,lp->nbindex,&(lp->re),&(lp->se),lp->pivots,err);
      else
         DualSimplexMinimize(lp->m,lp->d,lp->A,lp->B,lp->equalityset,lp->objrow,lp->rhscol,
           &(lp->LPS),&(lp->optvalue),lp->sol,lp->dsol,lp->nbindex,&(lp->re),&(lp->se),lp->pivots,err);
      break;

    case LPnone: *err=NoLPObjective; break;
  }
  time(&lp->endtime);
  lp->total_pivots=0;
  for (i=0; i<=3; i++) lp->total_pivots+=lp->pivots[i];
  if (*err==NoError) found=TRUE;
  return found;
}

dd_LPPtr dd_MakeLPforInteriorFinding(dd_LPPtr lp)
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
  dd_rowrange m;
  dd_colrange d;
  dd_NumberType numbtype;
  dd_LPObjectiveType obj;
  dd_LPType *lpnew;
  dd_rowrange i; 
  dd_colrange j;
  mytype bm,bmax,bceil;
  int localdebug=FALSE;

  dd_init(bm); dd_init(bmax); dd_init(bceil);
  dd_add(bm,dd_one,dd_one); dd_set(bmax,dd_one);
  numbtype=lp->numbtype;
  m=lp->m+1;
  d=lp->d+1;
  obj=LPmax;

  lpnew=CreateLPData(obj, numbtype, m, d);

  for (i=1; i<=lp->m; i++) {
    if (dd_Larger(lp->A[i-1][lp->rhscol-1],bmax)) 
      dd_set(bmax,lp->A[i-1][lp->rhscol-1]);
  }
  dd_mul(bceil,bm,bmax);
  if (localdebug) {printf("bceil is set to "); dd_WriteNumber(stdout, bceil);}
  
  for (i=1; i <= lp->m; i++) {
    for (j=1; j <= lp->d; j++) {
      dd_set(lpnew->A[i-1][j-1],lp->A[i-1][j-1]);
    }
  }

  for (i=1;i<=lp->m; i++){
    dd_neg(lpnew->A[i-1][lp->d],dd_one);  /* new column with all minus one's */
  }

  for (j=1;j<=lp->d;j++){
    dd_set(lpnew->A[m-2][j-1],dd_purezero);   /* new row (bceil, 0,...,0,-1) */
  }
  dd_set(lpnew->A[m-2][0],bceil);  /* new row (bceil, 0,...,0,-1) */

  for (j=1;j<= d-1;j++) {
    dd_set(lpnew->A[m-1][j-1],dd_purezero);  /* new obj row with (0,...,0,1) */
  }
  dd_set(lpnew->A[m-1][d-1],dd_one);    /* new obj row with (0,...,0,1) */
 
  if (localdebug) dd_WriteAmatrix(stdout, lp->A, lp->m, lp->d);
  if (localdebug) dd_WriteAmatrix(stdout, lpnew->A, lpnew->m, lpnew->d);
  dd_clear(bm); dd_clear(bmax); dd_clear(bceil);

  return lpnew;
}

void dd_WriteLPResult(FILE *f,dd_LPPtr lp,dd_ErrorType err)
{
  long j;

  fprintf(f,"* cdd LP solver result\n");
  
  if (err!=NoError) {
    dd_WriteErrorMessages(f,err);
    goto _L99;
  }

  dd_WriteProgramDescription(f);

  fprintf(f,"* #constraints = %ld\n",lp->m-1);
  fprintf(f,"* #variables   = %ld\n",lp->d-1);

  switch (lp->solver) {
    case DualSimplex:
      fprintf(f,"* Algorithm: dual simplex algorithm\n");break; 
    case CrissCross:
      fprintf(f,"* Algorithm: criss-cross method\n");break;
  }

  switch (lp->objective) {
    case LPmax:
      fprintf(f,"* maximization is chosen\n");break; 
    case LPmin:
      fprintf(f,"* minimization is chosen\n");break;
    case LPnone:
      fprintf(f,"* no objective type (max or min) is chosen\n");break;
  }
  
  if (lp->objective==LPmax||lp->objective==LPmin){
    fprintf(f,"* Objective function is\n");  
    for (j=0; j<lp->d; j++){
      if (j>0 && lp->A[lp->objrow-1][j]>=0 ) fprintf(f," +");
      if (j>0 && (j % 5) == 0) fprintf(f,"\n");
      dd_WriteNumber(f,lp->A[lp->objrow-1][j]);
      if (j>0) fprintf(f," X[%3ld]",j);
    }
    fprintf(f,"\n");
  }

  switch (lp->LPS){
  case Optimal:
    fprintf(f,"* LP status: a dual pair (x,y) of optimal solutions found.\n");
    fprintf(f,"begin\n");
    fprintf(f,"  primal_solution\n");
    for (j=1; j<lp->d; j++) {
      fprintf(f,"  %3ld : ",j);
      dd_WriteNumber(f,lp->sol[j]);
      fprintf(f,"\n");
    }
    fprintf(f,"  dual_solution\n");
    for (j=1; j<lp->d; j++){
      if (lp->nbindex[j+1]>0) {
        fprintf(f,"  %3ld : ",lp->nbindex[j+1]);
        dd_WriteNumber(f,lp->dsol[j]); fprintf(f,"\n");
      }
    }
    fprintf(f,"  optimal_value : "); dd_WriteNumber(f,lp->optvalue);
    fprintf(f,"\nend\n");
    break;

  case Inconsistent:
    fprintf(f,"* LP status: LP is inconsistent.\n");
    fprintf(f,"* The positive combination of original inequalities with\n");
    fprintf(f,"* the following coefficients will prove the inconsistency.\n");
    fprintf(f,"begin\n");
    fprintf(f,"  dual_direction\n");
    fprintf(f,"  %3ld : ",lp->re);
    dd_WriteNumber(f,dd_one);  fprintf(f,"\n");
    for (j=1; j<lp->d; j++){
      if (lp->nbindex[j+1]>0) {
        fprintf(f,"  %3ld : ",lp->nbindex[j+1]);
        dd_WriteNumber(f,lp->dsol[j]); fprintf(f,"\n");
      }
    }
    fprintf(f,"end\n");
    break;

  case DualInconsistent: case StrucDualInconsistent:
    fprintf(f,"* LP status: LP is dual inconsistent.\n");
    fprintf(f,"* The linear combination of columns with\n");
    fprintf(f,"* the following coefficients will prove the dual inconsistency.\n");
    fprintf(f,"* (It is also an unbounded direction for the primal LP.)\n");
    fprintf(f,"begin\n");
    fprintf(f,"  primal_direction\n");
    for (j=1; j<lp->d; j++) {
      fprintf(f,"  %3ld : ",j);
      dd_WriteNumber(f,lp->sol[j]);
      fprintf(f,"\n");
    }
    fprintf(f,"end\n");
    break;

  default:
    break;
  }
  fprintf(f,"* number of pivot operations = %ld (ph0 = %ld, ph1 = %ld, ph2 = %ld, ph3 = %ld)\n",lp->total_pivots,lp->pivots[0],lp->pivots[1],lp->pivots[2],lp->pivots[3]);
  dd_WriteLPTimes(f, lp);
_L99:;
}


/* end of dplex.c */


