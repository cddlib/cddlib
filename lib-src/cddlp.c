/* cddlp.c:  dual simplex method c-code
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.93, July 18, 2003
*/

/* cddlp.c : C-Implementation of the dual simplex method for
   solving an LP: max/min  A_(m-1).x subject to  x in P, where
   P= {x :  A_i.x >= 0, i=0,...,m-2, and  x_0=1}, and
   A_i is the i-th row of an m x n matrix A.
   Please read COPYING (GNU General Public Licence) and
   the manual cddlibman.tex for detail.
*/

#include "setoper.h"  /* set operation library header (Ver. May 18, 2000 or later) */
#include "cdd.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#if defined GMPRATIONAL
#include "cdd_f.h"
#endif

#define dd_CDDLPVERSION  "Version 0.93 (July 18, 2003)"

#define dd_FALSE 0
#define dd_TRUE 1

typedef set_type rowset;  /* set_type defined in setoper.h */
typedef set_type colset;

void dd_CrissCrossSolve(dd_LPPtr lp,dd_ErrorType *);
void dd_DualSimplexSolve(dd_LPPtr lp,dd_ErrorType *);
void dd_CrissCrossMinimize(dd_rowrange,dd_colrange,dd_Amatrix,dd_Bmatrix T,dd_rowset,
  dd_rowrange,dd_colrange,
  dd_LPStatusType *,mytype *optvalue,dd_Arow,dd_Arow,dd_colindex,
  dd_rowrange *,dd_colrange *,long *,dd_ErrorType *);
void dd_CrissCrossMaximize(dd_rowrange,dd_colrange,dd_Amatrix,dd_Bmatrix T,dd_rowset,
  dd_rowrange,dd_colrange,
  dd_LPStatusType *,mytype *optvalue,dd_Arow,dd_Arow,dd_colindex,
  dd_rowrange *,dd_colrange *,long *,dd_ErrorType *);
void dd_DualSimplexMinimize(dd_rowrange,dd_colrange, dd_Amatrix,dd_Bmatrix,dd_rowset,
  dd_rowrange,dd_colrange,
  dd_LPStatusType *,mytype *,dd_Arow,dd_Arow,dd_colindex,
  dd_rowrange *,dd_colrange *,long *,dd_ErrorType *);
void dd_DualSimplexMaximize(dd_rowrange,dd_colrange,dd_Amatrix,dd_Bmatrix,dd_rowset,
  dd_rowrange,dd_colrange,
  dd_LPStatusType *,mytype *,dd_Arow,dd_Arow,dd_colindex,
  dd_rowrange *,dd_colrange *,long *,dd_ErrorType *);
void dd_FindLPBasis(dd_rowrange,dd_colrange,dd_Amatrix,dd_Bmatrix,dd_rowindex,dd_rowset,
    dd_colindex,dd_rowindex,dd_rowrange,dd_colrange,
    dd_colrange *,int *,dd_LPStatusType *,long *);
void dd_FindDualFeasibleBasis(dd_rowrange,dd_colrange,dd_Amatrix,dd_Bmatrix,dd_rowindex,
    dd_colindex,long *,dd_rowrange,dd_colrange,
    dd_colrange *,dd_ErrorType *,dd_LPStatusType *,long *, long maxpivots);


#ifdef GMPRATIONAL
void dd_BasisStatus(ddf_LPPtr lpf, dd_LPPtr lp, dd_boolean*);
void dd_BasisStatusMinimize(dd_rowrange,dd_colrange, dd_Amatrix,dd_Bmatrix,dd_rowset,
    dd_rowrange,dd_colrange,ddf_LPStatusType, mytype *,dd_Arow,dd_Arow,ddf_colindex,
    ddf_rowrange,ddf_colrange,long *, int *, int *);
void dd_BasisStatusMaximize(dd_rowrange,dd_colrange, dd_Amatrix,dd_Bmatrix,dd_rowset,
    dd_rowrange,dd_colrange,ddf_LPStatusType, mytype *,dd_Arow,dd_Arow,ddf_colindex,
    ddf_rowrange,ddf_colrange,long *, int *, int *);
#endif

void dd_WriteBmatrix(FILE *f,dd_colrange d_size,dd_Bmatrix T);
void dd_SetNumberType(char *line,dd_NumberType *number,dd_ErrorType *Error);
void dd_ComputeRowOrderVector2(dd_rowrange m_size,dd_colrange d_size,dd_Amatrix A,
    dd_rowindex OV,dd_RowOrderType ho,unsigned int rseed);
void dd_SelectPreorderedNext2(dd_rowrange m_size,dd_colrange d_size,
    rowset excluded,dd_rowindex OV,dd_rowrange *hnext);
void dd_SetSolutions(dd_rowrange,dd_colrange,
   dd_Amatrix,dd_Bmatrix,dd_rowrange,dd_colrange,dd_LPStatusType,
   mytype *,dd_Arow,dd_Arow,dd_colindex,dd_rowrange,dd_colrange);


dd_LPSolutionPtr dd_CopyLPSolution(dd_LPPtr lp)
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
  dd_InitializeArow(lp->d+1,&(lps->sol));
  dd_InitializeArow(lp->d+1,&(lps->dsol));
  lps->nbindex=(long*) calloc((lp->d)+1,sizeof(long));  /* dual solution */
  for (j=0; j<=lp->d; j++){
    dd_set(lps->sol[j],lp->sol[j]);
    dd_set(lps->dsol[j],lp->dsol[j]);
    lps->nbindex[j]=lp->nbindex[j];
  }
  lps->pivots[0]=lp->pivots[0];
  lps->pivots[1]=lp->pivots[1];
  lps->pivots[2]=lp->pivots[2];
  lps->pivots[3]=lp->pivots[3];
  lps->pivots[4]=lp->pivots[4];
  lps->total_pivots=lp->total_pivots;

  return lps;
}


dd_LPPtr dd_CreateLPData(dd_LPObjectiveType obj,
   dd_NumberType nt,dd_rowrange m,dd_colrange d)
{
  dd_LPType *lp;

  lp=(dd_LPPtr) calloc(1,sizeof(dd_LPType));
  lp->solver=dd_DualSimplex;  /* set the default lp solver */
  lp->d=d;
  lp->m=m;
  lp->numbtype=nt;
  lp->objrow=m;
  lp->rhscol=1L;
  lp->objective=dd_LPnone;
  lp->LPS=dd_LPSundecided;
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
  dd_boolean localdebug=dd_FALSE;

  *err=dd_NoError;
  linc=set_card(M->linset);
  m=M->rowsize+1+linc; 
     /* We represent each equation by two inequalities.
        This is not the best way but makes the code simple. */
  d=M->colsize;
  if (localdebug) fprintf(stderr,"number of equalities = %ld\n", linc);
  
  lp=dd_CreateLPData(M->objective, M->numbtype, m, d);
  lp->Homogeneous = dd_TRUE;
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
      if (localdebug) fprintf(stderr,"equality row %ld generates the reverse row %ld.\n",i,irev);
    }
    for (j = 1; j <= M->colsize; j++) {
      dd_set(lp->A[i-1][j-1],M->matrix[i-1][j-1]);
      if (j==1 && i<M->rowsize && dd_Nonzero(M->matrix[i-1][j-1])) lp->Homogeneous = dd_FALSE;
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
    dd_clear(lp->optvalue);
    dd_FreeArow(lp->d_alloc,lp->dsol);
    dd_FreeArow(lp->d_alloc,lp->sol);
    dd_FreeBmatrix(lp->d_alloc,lp->B);
    dd_FreeAmatrix(lp->m_alloc,lp->d_alloc,lp->A);
    set_free(lp->equalityset);
    free(lp->nbindex);
    free(lp->given_nbindex);
    free(lp);
  }
}

void dd_FreeLPSolution(dd_LPSolutionPtr lps)
{
  if (lps!=NULL){
    free(lps->nbindex);
    dd_FreeArow(lps->d+1,lps->dsol);
    dd_FreeArow(lps->d+1,lps->sol);
    dd_clear(lps->optvalue);
    
    free(lps);
  }
}

int dd_LPReverseRow(dd_LPPtr lp, dd_rowrange i)
{
  dd_colrange j;
  int success=0;

  if (i>=1 && i<=lp->m){
    lp->LPS=dd_LPSundecided;
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
    lp->LPS=dd_LPSundecided;
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
    *number = dd_Integer;
    return;
  }
  else if (strncmp(line,"rational",8)==0) {
    *number = dd_Rational;
    return;
  }
  else if (strncmp(line,"real",4)==0) {
    *number = dd_Real;
    return;
  }
  else { 
    *number=dd_Unknown;
    *Error=dd_ImproperInputFormat;
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
  fprintf(f," %ld  %ld  real\n",m_size,d_size);
  fprintf(f,"          |");
  for (j=1; j<= d_size; j++) {
    fprintf(f," %ld",nbindex[j]);
  } fprintf(f,"\n");
  for (j=1; j<= d_size+1; j++) {
    fprintf(f," ----");
  } fprintf(f,"\n");
  for (i=1; i<= m_size; i++) {
    fprintf(f," %3ld(%3ld) |",i,bflag[i]);  
    for (j=1; j<= d_size; j++) {
      dd_TableauEntry(&x,m_size,d_size,A,T,i,j);
      dd_WriteNumber(f,x);
    }
    fprintf(f,"\n");
  }
  fprintf(f,"end\n");
  dd_clear(x);
}


void dd_SelectDualSimplexPivot(dd_rowrange m_size,dd_colrange d_size,
    int Phase1,dd_Amatrix A,dd_Bmatrix T,dd_rowindex OV,
    dd_colindex nbindex,dd_rowindex bflag,
    dd_rowrange objrow,dd_colrange rhscol,
    dd_rowrange *r,dd_colrange *s,int *selected,dd_LPStatusType *lps)
{ 
  /* selects a dual simplex pivot (*r,*s) if the current
     basis is dual feasible and not optimal. If not dual feasible,
     the procedure returns *selected=dd_FALSE and *lps=LPSundecided.
     If Phase1=dd_TRUE, the RHS column will be considered as the negative
     of the column of the largest variable (==m_size).  For this case, it is assumed
     that the caller used the auxiliary row (with variable m_size) to make the current
     dictionary dual feasible before calling this routine so that the nonbasic
     column for m_size corresponds to the auxiliary variable.
  */
  int colselected=dd_FALSE,rowselected=dd_FALSE,
    dualfeasible=dd_TRUE,localdebug=dd_FALSE;
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
  *selected=dd_FALSE;
  *lps=dd_LPSundecided;
  for (j=1; j<=d_size; j++){
    if (j!=rhscol){
      dd_TableauEntry(&(rcost[j-1]),m_size,d_size,A,T,objrow,j);
      if (dd_Positive(rcost[j-1])) { 
        dualfeasible=dd_FALSE;
      }
    }
  }
  if (dualfeasible){
    while ((*lps==dd_LPSundecided) && (!rowselected) && (!colselected)) {
      for (i=1; i<=m_size; i++) {
        if (i!=objrow && bflag[i]==-1) {  /* i is a basic variable */
          if (Phase1){
            dd_TableauEntry(&val, m_size,d_size,A,T,i,bflag[m_size]);
            dd_neg(val,val);
            /* for dual Phase I */
          } 
          else {dd_TableauEntry(&val,m_size,d_size,A,T,i,rhscol);}
          if (dd_Smaller(val,minval)) {
            *r=i;
            dd_set(minval,val);
          }
        }
      }
      if (dd_Nonnegative(minval)) {
        *lps=dd_Optimal;
      }
      else {
        rowselected=dd_TRUE;
        for (j=1; j<=d_size; j++){
          dd_TableauEntry(&val,m_size,d_size,A,T,*r,j);
          if (j!=rhscol && dd_Positive(val)) {
            dd_div(rat,rcost[j-1],val);
            dd_neg(rat,rat);
            if (*s==0 || dd_Smaller(rat,minrat)){
              dd_set(minrat,rat);
              *s=j;
            }
          }
        }
        if (*s>0) {colselected=dd_TRUE; *selected=dd_TRUE;}
        else *lps=dd_Inconsistent;
      }
    } /* end of while */
  }
  if (localdebug) {
     if (Phase1) fprintf(stderr,"Phase 1 : select %ld,%ld\n",*r,*s);
     else fprintf(stderr,"Phase 2 : select %ld,%ld\n",*r,*s);
  }
  dd_clear(val); dd_clear(minval); dd_clear(rat); dd_clear(minrat);
}

void dd_TableauEntry(mytype *x,dd_rowrange m_size, dd_colrange d_size, dd_Amatrix X, dd_Bmatrix T,
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

void dd_SelectPivot2(dd_rowrange m_size,dd_colrange d_size,dd_Amatrix A,dd_Bmatrix T,
            dd_RowOrderType roworder,dd_rowindex ordervec, rowset equalityset,
            dd_rowrange rowmax,rowset NopivotRow,
            colset NopivotCol,dd_rowrange *r,dd_colrange *s,
            dd_boolean *selected)
/* Select a position (*r,*s) in the matrix A.T such that (A.T)[*r][*s] is nonzero
   The choice is feasible, i.e., not on NopivotRow and NopivotCol, and
   best with respect to the specified roworder 
 */
{
  int stop;
  dd_rowrange i,rtemp;
  rowset rowexcluded;
  mytype Xtemp;
  dd_boolean localdebug=dd_FALSE;

  stop = dd_FALSE;
  localdebug=dd_debug;
  dd_init(Xtemp);
  set_initialize(&rowexcluded,m_size);
  set_copy(rowexcluded,NopivotRow);
  for (i=rowmax+1;i<=m_size;i++) {
    set_addelem(rowexcluded,i);   /* cannot pivot on any row > rmax */
  }
  *selected = dd_FALSE;
  do {
    rtemp=0; i=1;
    while (i<=m_size && rtemp==0) {  /* equalityset vars have highest priorities */
      if (set_member(i,equalityset) && !set_member(i,rowexcluded)){
        if (localdebug) fprintf(stderr,"marked set %ld chosen as a candidate\n",i);
        rtemp=i;
      }
      i++;
    }
    if (rtemp==0) dd_SelectPreorderedNext2(m_size,d_size,rowexcluded,ordervec,&rtemp);;
    if (rtemp>=1) {
      *r=rtemp;
      *s=1;
      while (*s <= d_size && !*selected) {
        dd_TableauEntry(&Xtemp,m_size,d_size,A,T,*r,*s);
        if (!set_member(*s,NopivotCol) && dd_Nonzero(Xtemp)) {
          *selected = dd_TRUE;
          stop = dd_TRUE;
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
      stop = dd_TRUE;
    }
  } while (!stop);
  set_free(rowexcluded); dd_clear(Xtemp);
}

void dd_GaussianColumnPivot(dd_rowrange m_size, dd_colrange d_size, 
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
  dd_boolean localdebug=dd_debug;

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
    dd_TableauEntry(&(Rtemp[j-1]), m_size, d_size, X, T, r,j);
  }
  dd_set(Xtemp0,Rtemp[s-1]);
  if (localdebug) {
    fprintf(stderr,"Gaussian Pivot: pivot entry = "); dd_WriteNumber(stderr,Xtemp0);
    fprintf(stderr,"\n");
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
          dd_WriteNumber(stderr, T[j1-1][j-1]);
        }
      }
      if (localdebug) fprintf(stderr,"\n");
    }
    if (localdebug) fprintf(stderr,"\n");
  }
  for (j = 1; j <= d_size; j++)
    dd_div(T[j-1][s - 1],T[j-1][s - 1],Xtemp0);

  dd_clear(Xtemp0); dd_clear(Xtemp1); dd_clear(Xtemp);
}

void dd_GaussianColumnPivot2(dd_rowrange m_size,dd_colrange d_size,
    dd_Amatrix A,dd_Bmatrix T,dd_colindex nbindex,dd_rowindex bflag,dd_rowrange r,dd_colrange s)
/* Update the Transformation matrix T with the pivot operation on (r,s) 
   This procedure performs a implicit pivot operation on the matrix A by
   updating the dual basis inverse  T.
 */
{
  int localdebug=dd_FALSE;
  long entering;

  dd_GaussianColumnPivot(m_size,d_size,A,T,r,s);
  entering=nbindex[s];
  bflag[r]=s;     /* the nonbasic variable r corresponds to column s */
  nbindex[s]=r;   /* the nonbasic variable on s column is r */
  if (localdebug) {
    fprintf(stderr,"Column pivot: (leaving, entering) = (%ld, %ld)\n", r,entering);
    fprintf(stderr, "bflag[%ld] is set to %ld\n", r, s);
    if (d_size <=20 && m_size <40){
      dd_WriteBmatrix(stderr,d_size,T);
      dd_WriteTableau(stderr,m_size,d_size,A,T,nbindex,bflag);
    }
  }

  if (entering>0) bflag[entering]=-1;
     /* original variables have negative index and should not affect the row index */
}


void dd_ResetTableau(dd_rowrange m_size,dd_colrange d_size,dd_Bmatrix T,
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
  int colselected=dd_FALSE,rowselected=dd_FALSE;
  dd_rowrange i;
  mytype val;
  
  dd_init(val);
  *selected=dd_FALSE;
  *lps=dd_LPSundecided;
  while ((*lps==dd_LPSundecided) && (!rowselected) && (!colselected)) {
    for (i=1; i<=m_size; i++) {
      if (i!=objrow && bflag[i]==-1) {  /* i is a basic variable */
        dd_TableauEntry(&val,m_size,d_size,A,T,i,rhscol);
        if (dd_Negative(val)) {
          rowselected=dd_TRUE;
          *r=i;
          break;
        }
      }
      else if (bflag[i] >0) { /* i is nonbasic variable */
        dd_TableauEntry(&val,m_size,d_size,A,T,objrow,bflag[i]);
        if (dd_Positive(val)) {
          colselected=dd_TRUE;
          *s=bflag[i];
          break;
        }
      }
    }
    if  ((!rowselected) && (!colselected)) {
      *lps=dd_Optimal;
      return;
    }
    else if (rowselected) {
     for (i=1; i<=m_size; i++) {
       if (bflag[i] >0) { /* i is nonbasic variable */
          dd_TableauEntry(&val,m_size,d_size,A,T,*r,bflag[i]);
          if (dd_Positive(val)) {
            colselected=dd_TRUE;
            *s=bflag[i];
            *selected=dd_TRUE;
            break;
          }
        }
      }
    }
    else if (colselected) {
      for (i=1; i<=m_size; i++) {
        if (i!=objrow && bflag[i]==-1) {  /* i is a basic variable */
          dd_TableauEntry(&val,m_size,d_size,A,T,i,*s);
          if (dd_Negative(val)) {
            rowselected=dd_TRUE;
            *r=i;
            *selected=dd_TRUE;
            break;
          }
        }
      }
    }
    if (!rowselected) {
      *lps=dd_DualInconsistent;
    }
    else if (!colselected) {
      *lps=dd_Inconsistent;
    }
  }
  dd_clear(val);
}

void dd_CrissCrossSolve(dd_LPPtr lp, dd_ErrorType *err)
{
  switch (lp->objective) {
    case dd_LPmax:
         dd_CrissCrossMaximize(lp->m,lp->d,lp->A,lp->B,lp->equalityset,lp->objrow,lp->rhscol,
           &(lp->LPS),&(lp->optvalue),lp->sol,lp->dsol,lp->nbindex,&(lp->re),&(lp->se),lp->pivots,err);
      break;
      
    case dd_LPmin:
         dd_CrissCrossMinimize(lp->m,lp->d,lp->A,lp->B,lp->equalityset,lp->objrow,lp->rhscol,
           &(lp->LPS),&(lp->optvalue),lp->sol,lp->dsol,lp->nbindex,&(lp->re),&(lp->se),lp->pivots,err);
      break;

    case dd_LPnone: *err=dd_NoLPObjective; break;
  }

}

void dd_DualSimplexSolve(dd_LPPtr lp, dd_ErrorType *err)
{
  switch (lp->objective) {
    case dd_LPmax:
         dd_DualSimplexMaximize(lp->m,lp->d,lp->A,lp->B,lp->equalityset,lp->objrow,lp->rhscol,
           &(lp->LPS),&(lp->optvalue),lp->sol,lp->dsol,lp->nbindex,&(lp->re),&(lp->se),lp->pivots,err);
      break;
      
    case dd_LPmin:
         dd_DualSimplexMinimize(lp->m,lp->d,lp->A,lp->B,lp->equalityset,lp->objrow,lp->rhscol,
           &(lp->LPS),&(lp->optvalue),lp->sol,lp->dsol,lp->nbindex,&(lp->re),&(lp->se),lp->pivots,err);
      break;

    case dd_LPnone: *err=dd_NoLPObjective; break;
  }
}

#ifdef GMPRATIONAL

dd_LPStatusType LPSf2LPS(ddf_LPStatusType lpsf)
{
   dd_LPStatusType lps=dd_LPSundecided;

   switch (lpsf) {
   case ddf_LPSundecided: lps=dd_LPSundecided; break;
   case ddf_Optimal: lps=dd_Optimal; break;
   case ddf_Inconsistent: lps=dd_Inconsistent; break; 
   case ddf_DualInconsistent: lps=dd_DualInconsistent; break;
   case ddf_StrucInconsistent: lps=dd_StrucInconsistent; break; 
   case ddf_StrucDualInconsistent: lps=dd_StrucDualInconsistent; break;
   case ddf_Unbounded: lps=dd_Unbounded; break;
   case ddf_DualUnbounded: lps=dd_DualUnbounded; break;
   }
   return lps;
}


void dd_BasisStatus(ddf_LPPtr lpf, dd_LPPtr lp, dd_boolean *LPScorrect)
{
  int i;
  dd_colrange j;
  dd_boolean basisfound; 
 
  switch (lp->objective) {
    case dd_LPmax:
      dd_BasisStatusMaximize(lp->m,lp->d,lp->A,lp->B,lp->equalityset,lp->objrow,lp->rhscol,
           lpf->LPS,&(lp->optvalue),lp->sol,lp->dsol,lpf->nbindex,lpf->re,lpf->se,lp->pivots, 
           &basisfound, LPScorrect);
      if (*LPScorrect) {
         /* printf("BasisStatus Check: the current basis is verified with GMP\n"); */
         lp->LPS=LPSf2LPS(lpf->LPS);
         lp->re=lpf->re;
         lp->se=lpf->se;
         for (i=1; i<=5; i++) lp->pivots[i-1]+=lpf->pivots[i-1]; 
         for (j=1; j<=lp->d; j++) lp->nbindex[j]=lpf->nbindex[j]; 
      }
      break;
    case dd_LPmin:
      dd_BasisStatusMinimize(lp->m,lp->d,lp->A,lp->B,lp->equalityset,lp->objrow,lp->rhscol,
           lpf->LPS,&(lp->optvalue),lp->sol,lp->dsol,lpf->nbindex,lpf->re,lpf->se,lp->pivots, 
           &basisfound, LPScorrect);
      if (*LPScorrect) {
         /* printf("BasisStatus Check: the current basis is verified with GMP\n"); */
         lp->LPS=LPSf2LPS(lpf->LPS);
         lp->re=lpf->re;
         lp->se=lpf->se;
         for (i=1; i<=4; i++) lp->pivots[i-1]+=lpf->pivots[i-1]; 
         for (j=1; j<=lp->d; j++) lp->nbindex[j]=lpf->nbindex[j]; 
      }
      break;
    case dd_LPnone:  break;
   }      
}
#endif

void dd_CrissCrossMinimize(dd_rowrange m_size,dd_colrange d_size,
    dd_Amatrix A,dd_Bmatrix T,dd_rowset equalityset,
    dd_rowrange objrow,dd_colrange rhscol,dd_LPStatusType *LPS,
    mytype *optvalue,dd_Arow sol,dd_Arow dsol,dd_colindex nbindex,
    dd_rowrange *re,dd_colrange *se,long *pivots,dd_ErrorType *err)
{
   dd_colrange j;

   *err=dd_NoError;
   for (j=1; j<=d_size; j++)
     dd_neg(A[objrow-1][j-1],A[objrow-1][j-1]);
   dd_CrissCrossMaximize(m_size,d_size,A,T,equalityset, objrow,rhscol,
     LPS,optvalue,sol,dsol,nbindex,re,se,pivots,err);
   dd_neg(*optvalue,*optvalue);
   for (j=1; j<=d_size; j++){
     dd_neg(dsol[j-1],dsol[j-1]);
     dd_neg(A[objrow-1][j-1],A[objrow-1][j-1]);
   }
}

void dd_CrissCrossMaximize(dd_rowrange m_size,dd_colrange d_size,
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
  dd_colindex nbtemp;

  *err=dd_NoError;
  nbtemp=(long *) calloc(d_size+1,sizeof(long*));
  for (i=0; i<= 4; i++) pivots[i]=0;
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
  dd_ComputeRowOrderVector2(m_size,d_size,A,OrderVector,dd_MinIndex,rseed);

  *re=0; *se=0; pivots1=0;

  dd_ResetTableau(m_size,d_size,T,nbindex,bflag,objrow,rhscol);

  dd_FindLPBasis(m_size,d_size,A,T,OrderVector, equalityset,
      nbindex,bflag, objrow,rhscol,&s,&found,LPS,&pivots0);
  pivots[0]=pivots0;

  if (!found){
     *se=s;
     goto _L99;
     /* No LP basis is found, and thus Inconsistent.  
     Output the evidence column. */
  }

  stop=dd_FALSE;
  do {   /* Criss-Cross Method */
    dd_SelectCrissCrossPivot(m_size,d_size,A,T,bflag,
       objrow,rhscol,&r,&s,&chosen,LPS);
    if (chosen) {
      dd_GaussianColumnPivot2(m_size,d_size,A,T,nbindex,bflag,r,s);
      pivots1++;
    } else {
      switch (*LPS){
        case dd_Inconsistent: *re=r;
        case dd_DualInconsistent: *se=s;
        default: break;
      }
      stop=dd_TRUE;
    }
  } while(!stop);
  pivots[1]=pivots1;
  
_L99:

  dd_SetSolutions(m_size,d_size,A,T,
   objrow,rhscol,*LPS,optvalue,sol,dsol,nbindex,*re,*se);
  free(nbtemp);

}

void dd_FindLPBasis(dd_rowrange m_size,dd_colrange d_size,
    dd_Amatrix A, dd_Bmatrix T,dd_rowindex OV,dd_rowset equalityset, dd_colindex nbindex,
    dd_rowindex bflag,dd_rowrange objrow,dd_colrange rhscol,
    dd_colrange *cs,int *found,dd_LPStatusType *lps,long *pivot_no)
{ 
  /* Find a LP basis using Gaussian pivots.
     If the problem has an LP basis,
     the procedure returns *found=dd_TRUE,*lps=LPSundecided and an LP basis.
     If the constraint matrix A (excluding the rhs and objective) is not
     column indepent, there are two cases.  If the dependency gives a dual
     inconsistency, this returns *found=dd_FALSE, *lps=dd_StrucDualInconsistent and 
     the evidence column *s.  Otherwise, this returns *found=dd_TRUE, 
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
  *found=dd_FALSE; *cs=0; rank=0;
  *lps=dd_LPSundecided;

  set_initialize(&RowSelected,m_size);
  set_initialize(&ColSelected,d_size);
  set_addelem(RowSelected,objrow);
  set_addelem(ColSelected,rhscol);

  stop=dd_FALSE;
  do {   /* Find a LP basis */
    dd_SelectPivot2(m_size,d_size,A,T,dd_MinIndex,OV,equalityset,
      m_size,RowSelected,ColSelected,&r,&s,&chosen);
    if (chosen) {
      set_addelem(RowSelected,r);
      set_addelem(ColSelected,s);
      dd_GaussianColumnPivot2(m_size,d_size,A,T,nbindex,bflag,r,s);
      pivots_p0++;
      rank++;
    } else {
      for (j=1;j<=d_size  && *lps==dd_LPSundecided; j++) {
        if (j!=rhscol && nbindex[j]<0){
          dd_TableauEntry(&val,m_size,d_size,A,T,objrow,j);
          if (dd_Nonzero(val)){  /* dual inconsistent */
            *lps=dd_StrucDualInconsistent;
            *cs=j;
            /* dual inconsistent because the nonzero reduced cost */
          }
        }
      }
      if (*lps==dd_LPSundecided) *found=dd_TRUE;  
         /* dependent columns but not dual inconsistent. */
      stop=dd_TRUE;
    }
    if (rank==d_size-1) {
      stop = dd_TRUE;
      *found=dd_TRUE;
    }
  } while (!stop);

  *pivot_no=pivots_p0;
  set_free(RowSelected);
  set_free(ColSelected);
  dd_clear(val);
}


void dd_FindLPBasis2(dd_rowrange m_size,dd_colrange d_size,
    dd_Amatrix A, dd_Bmatrix T,dd_rowindex OV,dd_rowset equalityset, dd_colindex nbindex,
    dd_rowindex bflag,dd_rowrange objrow,dd_colrange rhscol,
    dd_colrange *cs,int *found,long *pivot_no)
{ 
  /* Similar to dd_FindLPBasis but it is much simpler.  This tries to recompute T for
  the specified basis given by nbindex.  It will return *found=dd_FALSE if the specified
  basis is not a basis.
  */
  int chosen,stop;
  long pivots_p0=0,rank;
  dd_colset ColSelected;
  dd_rowset RowSelected, NopivotRow;
  mytype val;
  dd_boolean localdebug=dd_FALSE;

  dd_rowrange r,negcount=0;
  dd_colrange j,s;

  dd_init(val);
  *found=dd_FALSE; *cs=0; rank=0;

  set_initialize(&RowSelected,m_size);
  set_initialize(&ColSelected,d_size);
  set_initialize(&NopivotRow,m_size);
  set_addelem(RowSelected,objrow);
  set_addelem(ColSelected,rhscol);
  set_compl(NopivotRow, NopivotRow);  /* set NopivotRow to be the groundset */
  
  for (j=2; j<=d_size; j++) 
    if (nbindex[j]>0) 
       set_delelem(NopivotRow, nbindex[j]);
    else if (nbindex[j]<0) 
       negcount++;       
     
  set_uni(RowSelected, RowSelected, NopivotRow);  /* RowSelected is the set of rows not allowed to poviot on */

  stop=dd_FALSE;
  do {   /* Find a LP basis */
    dd_SelectPivot2(m_size,d_size,A,T,dd_MinIndex,OV,equalityset,
      m_size,RowSelected,ColSelected,&r,&s,&chosen);
    if (chosen) {
      set_addelem(RowSelected,r);
      set_addelem(ColSelected,s);

      dd_GaussianColumnPivot2(m_size,d_size,A,T,nbindex,bflag,r,s);
      if (localdebug && m_size <=10){
        dd_WriteBmatrix(stderr,d_size,T);
        dd_WriteTableau(stderr,m_size,d_size,A,T,nbindex,bflag);
      }
      pivots_p0++;
      rank++;
    } else{
      *found=dd_FALSE;   /* cannot pivot on any of the spacified positions. */
      stop=dd_TRUE;
    }
    if (rank==d_size-1-negcount) {
      stop = dd_TRUE;
      *found=dd_TRUE;
    }
  } while (!stop);

  for (j=1; j<=d_size; j++) if (nbindex[j]>0) bflag[nbindex[j]]=j;
  *pivot_no=pivots_p0;
  set_free(RowSelected);
  set_free(ColSelected);
  set_free(NopivotRow);
  dd_clear(val);
}

void dd_FindDualFeasibleBasis(dd_rowrange m_size,dd_colrange d_size,
    dd_Amatrix A,dd_Bmatrix T,dd_rowindex OV,
    dd_colindex nbindex,dd_rowindex bflag,dd_rowrange objrow,dd_colrange rhscol,
    dd_colrange *s,dd_ErrorType *err,dd_LPStatusType *lps,long *pivot_no, long maxpivots)
{ 
  /* Find a dual feasible basis using Phase I of Dual Simplex method.
     If the problem is dual feasible,
     the procedure returns *err=NoError, *lps=LPSundecided and a dual feasible
     basis.   If the problem is dual infeasible, this returns
     *err=NoError, *lps=DualInconsistent and the evidence column *s.
     Caution: matrix A must have at least one extra row:  the row space A[m_size] must
     have been allocated.
  */
  int phase1,dualfeasible=dd_TRUE,localdebug=dd_FALSE,chosen,stop;
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

  *err=dd_NoError; *lps=dd_LPSundecided; *s=0;
  local_m_size=m_size+1;  /* increase m_size by 1 */

  ms=0;  /* ms will be the index of column which has the largest reduced cost */
  for (j=1; j<=d_size; j++){
    if (j!=rhscol){
      if (localdebug) fprintf(stderr,"checking the column %ld var %ld\n",j,nbindex[j]); 
      dd_TableauEntry(&(rcost[j-1]),local_m_size,d_size,A,T,objrow,j);
      if (localdebug) {fprintf(stderr,"reduced cost = "); dd_WriteNumber(stderr, rcost[j-1]); }
      if (dd_Larger(rcost[j-1],maxcost)) {dd_set(maxcost,rcost[j-1]); ms = j;}
    }
  }
  if (dd_Positive(maxcost)) dualfeasible=dd_FALSE;

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
      fprintf(stderr,"dd_FindDualFeasibleBasis: curruent basis is not dual feasible.\n");
      fprintf(stderr,"because of the column %ld assoc. with var %ld   dual cost =",
       ms,nbindex[ms]);
      dd_WriteNumber(stderr, maxcost);
    }

    /* Pivot on (local_m_size,ms) so that the dual basic solution becomes feasible */
    dd_GaussianColumnPivot2(local_m_size,d_size,A,T,nbindex,bflag,local_m_size,ms);
    pivots_p1=pivots_p1+1;

    phase1=dd_TRUE; stop=dd_FALSE;
    do {   /* Dual Simplex Phase I */
      chosen=dd_FALSE; LPSphase1=dd_LPSundecided;
      if (pivots_p1>maxpivots) {
        *err=dd_LPCycling;
        fprintf(stderr,"max number %ld of pivots performed in Phase I. Switch to the anticycling phase.\n", maxpivots);
        goto _L99;  /* failure due to max no. of pivots performed */
      }
      dd_SelectDualSimplexPivot(local_m_size,d_size,phase1,A,T,OV,nbindex,bflag,
        objrow,rhscol,&r_val,&s_val,&chosen,&LPSphase1);
      if (!chosen) {
        /* The current dictionary is terminal.  There are two cases:
           dd_TableauEntry(local_m_size,d_size,A,T,objrow,ms) is negative or zero.
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
            dd_TableauEntry(&val,local_m_size,d_size,A,T,i,ms);  /* auxiliary column*/
            if (dd_Smaller(val, minval)) {
              r_val=i;
              dd_set(minval,val);
              if (localdebug) {
                fprintf(stderr,"update minval with = ");
                dd_WriteNumber(stderr, minval);
                fprintf(stderr,"  r_val = %ld\n",r_val);
              }
            }
          }
        }
        dd_clear(minval);

        dd_GaussianColumnPivot2(local_m_size,d_size,A,T,nbindex,bflag,r_val,ms);
        pivots_p1=pivots_p1+1;

        dd_TableauEntry(&x,local_m_size,d_size,A,T,objrow,ms);
        if (dd_Negative(x)){
          *err=dd_NoError; *lps=dd_DualInconsistent;  *s=ms;
        }
        stop=dd_TRUE;
      } else {
        dd_GaussianColumnPivot2(local_m_size,d_size,A,T,nbindex,bflag,r_val,s_val);
        pivots_p1=pivots_p1+1;
        if (bflag[local_m_size]<0) {
          stop=dd_TRUE; 
          if (localdebug) 
            fprintf(stderr,"Dual Phase I: the auxiliary variable entered the basis, go to phase II\n");
        }
      }
    } while(!stop);
  }
_L99:
  *pivot_no=pivots_p1;
  dd_clear(x); dd_clear(val); dd_clear(maxcost);
}

void dd_DualSimplexMinimize(dd_rowrange m_size,dd_colrange d_size,
   dd_Amatrix A,dd_Bmatrix T,dd_rowset equalityset,
   dd_rowrange objrow,dd_colrange rhscol,dd_LPStatusType *LPS,
   mytype *optvalue,dd_Arow sol,dd_Arow dsol,dd_colindex nbindex,
   dd_rowrange *re,dd_colrange *se,long *pivots,dd_ErrorType *err)
{
   dd_colrange j;

   *err=dd_NoError;
   for (j=1; j<=d_size; j++)
     dd_neg(A[objrow-1][j-1],A[objrow-1][j-1]);
   dd_DualSimplexMaximize(m_size,d_size,A,T,equalityset, objrow,rhscol,
     LPS,optvalue,sol,dsol,nbindex,re,se,pivots,err);
   dd_neg(*optvalue,*optvalue);
   for (j=1; j<=d_size; j++){
     dd_neg(dsol[j-1],dsol[j-1]);
     dd_neg(A[objrow-1][j-1],A[objrow-1][j-1]);
   }
}

void dd_DualSimplexMaximize(dd_rowrange m_size,dd_colrange d_size,
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
  long pivots_ds=0,pivots_p0=0,pivots_p1=0,pivots_pc=0,maxpivots,maxpivfactor=40;
  dd_rowrange i,r;
  dd_colrange s;
  static dd_rowindex bflag;
  static long mlast=0,nlast=0;
  static dd_rowindex OrderVector;  /* the permutation vector to store a preordered row indeces */
  unsigned int rseed=1;
  
  /* *err=dd_NoError; */
  for (i=0; i<= 4; i++) pivots[i]=0;
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
  dd_ComputeRowOrderVector2(m_size,d_size,A,OrderVector,dd_MinIndex,rseed);

  *re=0; *se=0;
  
  dd_ResetTableau(m_size,d_size,T,nbindex,bflag,objrow,rhscol);
   
  dd_FindLPBasis(m_size,d_size,A,T,OrderVector,equalityset,nbindex,bflag,
      objrow,rhscol,&s,&found,LPS,&pivots_p0);
  pivots[0]=pivots_p0;

  if (!found){
     *se=s;
     goto _L99;
     /* No LP basis is found, and thus Inconsistent.  
     Output the evidence column. */
  }

  dd_FindDualFeasibleBasis(m_size,d_size,A,T,OrderVector,nbindex,bflag,
      objrow,rhscol,&s, err, LPS,&pivots_p1, maxpivots);
  pivots[1]=pivots_p1;

  if (*err==dd_LPCycling){
    if (dd_debug) fprintf(stderr, "Phase I failed and thus switch to the Criss-Cross method\n");
    dd_CrissCrossMaximize(m_size,d_size,A,T,equalityset,
    objrow,rhscol,LPS,optvalue,sol,dsol,nbindex,re,se,pivots,err);
    return;
  }
  if (*LPS==dd_DualInconsistent){
     *se=s;
     goto _L99;
     /* No dual feasible basis is found, and thus DualInconsistent.  
     Output the evidence column. */
  }

  /* Dual Simplex Method */
  stop=dd_FALSE;
  do {
    chosen=dd_FALSE; *LPS=dd_LPSundecided; phase1=dd_FALSE;
    if (pivots_ds<maxpivots) {
      dd_SelectDualSimplexPivot(m_size,d_size,phase1,A,T,OrderVector,nbindex,bflag,
        objrow,rhscol,&r,&s,&chosen,LPS);
    }
    if (chosen) pivots_ds=pivots_ds+1;
    if (!chosen && *LPS==dd_LPSundecided) {  
      if (dd_debug) fprintf(stderr,"Warning: an emergency CC pivot in Phase II is performed\n");
      /* In principle this should not be executed because we already have dual feasibility
         attained and dual simplex pivot should have been chosen.  This might occur
         under floating point computation, or the case of cycling.
      */
      dd_SelectCrissCrossPivot(m_size,d_size,A,T,bflag,
        objrow,rhscol,&r,&s,&chosen,LPS);
      if (chosen) pivots_pc=pivots_pc+1;
    }
    if (chosen) {
      dd_GaussianColumnPivot2(m_size,d_size,A,T,nbindex,bflag,r,s);
    } else {
      switch (*LPS){
        case dd_Inconsistent: *re=r;
        case dd_DualInconsistent: *se=s;
        default: break;
      }
      stop=dd_TRUE;
    }
  } while(!stop);
  pivots[2]=pivots_ds;
  pivots[3]=pivots_pc;

_L99: 
  dd_SetSolutions(m_size,d_size,A,T,objrow,rhscol,*LPS,optvalue,sol,dsol,nbindex,*re,*se);
}

void dd_SetSolutions(dd_rowrange m_size,dd_colrange d_size,
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
  int localdebug=dd_FALSE;
  
  dd_init(x); dd_init(sw);
  switch (LPS){
  case dd_Optimal:
    for (j=1;j<=d_size; j++) {
      dd_set(sol[j-1],T[j-1][rhscol-1]);
      dd_TableauEntry(&x,m_size,d_size,A,T,objrow,j);
      dd_neg(dsol[j-1],x);
      dd_TableauEntry(optvalue,m_size,d_size,A,T,objrow,rhscol);
      if (localdebug) {fprintf(stderr,"dsol[%ld]= ",nbindex[j]); dd_WriteNumber(stderr, dsol[j-1]); }
    }
    break;
  case dd_Inconsistent:
    if (localdebug) fprintf(stderr,"DualSimplexSolve: LP is inconsistent.\n");
    for (j=1;j<=d_size; j++) {
      dd_set(sol[j-1],T[j-1][rhscol-1]);
      dd_TableauEntry(&x,m_size,d_size,A,T,re,j);
      dd_neg(dsol[j-1],x);
      if (localdebug) {fprintf(stderr,"dsol[%ld]= ",nbindex[j]); dd_WriteNumber(stderr,dsol[j-1]);}
    }
    break;
  case dd_DualInconsistent:
    for (j=1;j<=d_size; j++) {
      dd_set(sol[j-1],T[j-1][se-1]);
      dd_TableauEntry(&x,m_size,d_size,A,T,objrow,j);
      dd_neg(dsol[j-1],x);
      if (localdebug) {fprintf(stderr,"dsol[%ld]= \n",nbindex[j]);dd_WriteNumber(stderr,dsol[j-1]);}
    }
    if (localdebug) printf( "DualSimplexSolve: LP is dual inconsistent.\n");
    break;

  case dd_StrucDualInconsistent:
    dd_TableauEntry(&x,m_size,d_size,A,T,objrow,se);
    if (dd_Positive(x)) dd_set(sw,dd_one);
    else dd_neg(sw,dd_one);
    for (j=1;j<=d_size; j++) {
      dd_mul(sol[j-1],sw,T[j-1][se-1]);
      dd_TableauEntry(&x,m_size,d_size,A,T,objrow,j);
      dd_neg(dsol[j-1],x);
      if (localdebug) {fprintf(stderr,"dsol[%ld]= ",nbindex[j]);dd_WriteNumber(stderr,dsol[j-1]);}
    }
    if (localdebug) fprintf(stderr,"DualSimplexSolve: LP is dual inconsistent.\n");
    break;

  default:break;
  }
  dd_clear(x); dd_clear(sw);
}


void dd_RandomPermutation2(dd_rowindex OV,long t,unsigned int seed)
{
  long k,j,ovj;
  double u,xk,r,rand_max=(double) RAND_MAX;
  int localdebug=dd_FALSE;

  srand(seed);
  for (j=t; j>1 ; j--) {
    r=rand();
    u=r/rand_max;
    xk=(double)(j*u +1);
    k=(long)xk;
    if (localdebug) fprintf(stderr,"u=%g, k=%ld, r=%g, randmax= %g\n",u,k,r,rand_max);
    ovj=OV[j];
    OV[j]=OV[k];
    OV[k]=ovj;
    if (localdebug) fprintf(stderr,"row %ld is exchanged with %ld\n",j,k); 
  }
}

void dd_ComputeRowOrderVector2(dd_rowrange m_size,dd_colrange d_size,dd_Amatrix A,
    dd_rowindex OV,dd_RowOrderType ho,unsigned int rseed)
{
  long i,itemp;
  
  OV[0]=0;
  switch (ho){
  case dd_MaxIndex:
    for(i=1; i<=m_size; i++) OV[i]=m_size-i+1;
    break;

  case dd_LexMin:
    for(i=1; i<=m_size; i++) OV[i]=i;
    dd_QuickSort(OV,1,m_size,A,d_size);
   break;

  case dd_LexMax:
    for(i=1; i<=m_size; i++) OV[i]=i;
    dd_QuickSort(OV,1,m_size,A,d_size);
    for(i=1; i<=m_size/2;i++){   /* just reverse the order */
      itemp=OV[i];
      OV[i]=OV[m_size-i+1];
      OV[m_size-i+1]=itemp;
    }
    break;

  case dd_RandomRow:
    for(i=1; i<=m_size; i++) OV[i]=i;
    if (rseed<=0) rseed=1;
    dd_RandomPermutation2(OV,m_size,rseed);
    break;

  case dd_MinIndex: 
    for(i=1; i<=m_size; i++) OV[i]=i;
    break;

  default: 
    for(i=1; i<=m_size; i++) OV[i]=i;
    break;
 }
}

void dd_SelectPreorderedNext2(dd_rowrange m_size,dd_colrange d_size,
    rowset excluded,dd_rowindex OV,dd_rowrange *hnext)
{
  dd_rowrange i,k;
  
  *hnext=0;
  for (i=1; i<=m_size && *hnext==0; i++){
    k=OV[i];
    if (!set_member(k,excluded)) *hnext=k ;
  }
}

#ifdef GMPRATIONAL

ddf_LPObjectiveType Obj2Obj(dd_LPObjectiveType obj)
{
   ddf_LPObjectiveType objf=ddf_LPnone;

   switch (obj) {
   case dd_LPnone: objf=ddf_LPnone; break;
   case dd_LPmax: objf=ddf_LPmax; break;
   case dd_LPmin: objf=ddf_LPmin; break;
   }
   return objf;
}

ddf_LPPtr dd_LPgmp2LPf(dd_LPPtr lp)
{
  dd_rowrange i;
  dd_colrange j;
  ddf_LPType *lpf;
  double val;
  dd_boolean localdebug=dd_FALSE;

  if (localdebug) fprintf(stderr,"Converting a GMP-LP to a float-LP.\n");
  
  lpf=ddf_CreateLPData(Obj2Obj(lp->objective), ddf_Real, lp->m, lp->d);
  lpf->Homogeneous = lp->Homogeneous;
  lpf->eqnumber=lp->eqnumber;  /* this records the number of equations */

  for (i = 1; i <= lp->m; i++) {
    if (set_member(i, lp->equalityset)) set_addelem(lpf->equalityset,i);    
          /* it is equality. Its reversed row will not be in this set */
      for (j = 1; j <= lp->d; j++) {
        val=mpq_get_d(lp->A[i-1][j-1]);
        ddf_set_d(lpf->A[i-1][j-1],val);
      }  /*of j*/
  }  /*of i*/

  return lpf;
}


#endif


dd_boolean dd_LPSolve(dd_LPPtr lp,dd_LPSolverType solver,dd_ErrorType *err)
/* 
The current version of dd_LPSolve that solves an LP with floating-arithmetics first
and then with the specified arithimetics if it is GMP.

When LP is inconsistent then *re returns the evidence row.
When LP is dual-inconsistent then *se returns the evidence column.
*/
{
  int i;
  dd_boolean found=dd_FALSE;
#ifdef GMPRATIONAL
  ddf_LPPtr lpf;
  ddf_ErrorType errf;
  dd_boolean LPScorrect=dd_FALSE;
  dd_boolean localdebug=dd_FALSE;
#endif

  *err=dd_NoError;
  lp->solver=solver;
  time(&lp->starttime);

#ifndef GMPRATIONAL
  switch (lp->solver) {
    case dd_CrissCross:
      dd_CrissCrossSolve(lp,err);
      break;
    case dd_DualSimplex:
      dd_DualSimplexSolve(lp,err);
      break;
  }
#else
  lpf=dd_LPgmp2LPf(lp);
  switch (lp->solver) {
    case dd_CrissCross:
      ddf_CrissCrossSolve(lpf,&errf);                /* First, run with double float. */
      dd_BasisStatus(lpf,lp, &LPScorrect);    /* Check the basis. */
      if (!LPScorrect) {
         if (localdebug) printf("BasisStatus: the current basis is NOT verified with GMP. Rerun with GMP.\n");
         dd_CrissCrossSolve(lp,err);  /* Rerun with GMP if fails. */
      } else {
         if (localdebug) printf("BasisStatus: the current basis is verified with GMP. The LP Solved.\n");
      }
      break;
    case dd_DualSimplex:
      ddf_DualSimplexSolve(lpf,&errf);                /* First, run with double float. */
      dd_BasisStatus(lpf,lp, &LPScorrect);    /* Check the basis. */
      if (!LPScorrect){
         if (localdebug) printf("BasisStatus: the current basis is NOT verified with GMP. Rerun with GMP.\n");
         dd_DualSimplexSolve(lp,err);  /* Rerun with GMP if fails. */
         if (localdebug){
            printf("*total number pivots = %ld (ph0 = %ld, ph1 = %ld, ph2 = %ld, ph3 = %ld, ph4 = %ld)\n",
               lp->total_pivots,lp->pivots[0],lp->pivots[1],lp->pivots[2],lp->pivots[3],lp->pivots[4]);
            ddf_WriteLPResult(stdout, lpf, errf);
            dd_WriteLP(stdout, lp);
         }
      } else {
         if (localdebug) printf("BasisStatus: the current basis is verified with GMP. The LP Solved.\n");
      }
      break;
  }
  ddf_FreeLPData(lpf);
#endif

  time(&lp->endtime);
  lp->total_pivots=0;
  for (i=0; i<=4; i++) lp->total_pivots+=lp->pivots[i];
  if (*err==dd_NoError) found=dd_TRUE;
  return found;
}


dd_boolean dd_LPSolve0(dd_LPPtr lp,dd_LPSolverType solver,dd_ErrorType *err)
/* 
The original version of dd_LPSolve that solves an LP with specified arithimetics.

When LP is inconsistent then *re returns the evidence row.
When LP is dual-inconsistent then *se returns the evidence column.
*/
{
  int i;
  dd_boolean found=dd_FALSE;

  *err=dd_NoError;
  lp->solver=solver;
  time(&lp->starttime);

  switch (lp->solver) {
    case dd_CrissCross:
      dd_CrissCrossSolve(lp,err);
      break;
    case dd_DualSimplex:
      dd_DualSimplexSolve(lp,err);
      break;
  }

  time(&lp->endtime);
  lp->total_pivots=0;
  for (i=0; i<=4; i++) lp->total_pivots+=lp->pivots[i];
  if (*err==dd_NoError) found=dd_TRUE;
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
  int localdebug=dd_FALSE;

  dd_init(bm); dd_init(bmax); dd_init(bceil);
  dd_add(bm,dd_one,dd_one); dd_set(bmax,dd_one);
  numbtype=lp->numbtype;
  m=lp->m+1;
  d=lp->d+1;
  obj=dd_LPmax;

  lpnew=dd_CreateLPData(obj, numbtype, m, d);

  for (i=1; i<=lp->m; i++) {
    if (dd_Larger(lp->A[i-1][lp->rhscol-1],bmax)) 
      dd_set(bmax,lp->A[i-1][lp->rhscol-1]);
  }
  dd_mul(bceil,bm,bmax);
  if (localdebug) {fprintf(stderr,"bceil is set to "); dd_WriteNumber(stderr, bceil);}
  
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
 
  if (localdebug) dd_WriteAmatrix(stderr, lp->A, lp->m, lp->d);
  if (localdebug) dd_WriteAmatrix(stderr, lpnew->A, lpnew->m, lpnew->d);
  dd_clear(bm); dd_clear(bmax); dd_clear(bceil);

  return lpnew;
}


void dd_WriteLPResult(FILE *f,dd_LPPtr lp,dd_ErrorType err)
{
  long j;

  fprintf(f,"* cdd LP solver result\n");
  
  if (err!=dd_NoError) {
    dd_WriteErrorMessages(f,err);
    goto _L99;
  }

  dd_WriteProgramDescription(f);

  fprintf(f,"* #constraints = %ld\n",lp->m-1);
  fprintf(f,"* #variables   = %ld\n",lp->d-1);

  switch (lp->solver) {
    case dd_DualSimplex:
      fprintf(f,"* Algorithm: dual simplex algorithm\n");break; 
    case dd_CrissCross:
      fprintf(f,"* Algorithm: criss-cross method\n");break;
  }

  switch (lp->objective) {
    case dd_LPmax:
      fprintf(f,"* maximization is chosen\n");break; 
    case dd_LPmin:
      fprintf(f,"* minimization is chosen\n");break;
    case dd_LPnone:
      fprintf(f,"* no objective type (max or min) is chosen\n");break;
  }
  
  if (lp->objective==dd_LPmax||lp->objective==dd_LPmin){
    fprintf(f,"* Objective function is\n");  
    for (j=0; j<lp->d; j++){
      if (j>0 && dd_Nonnegative(lp->A[lp->objrow-1][j]) ) fprintf(f," +");
      if (j>0 && (j % 5) == 0) fprintf(f,"\n");
      dd_WriteNumber(f,lp->A[lp->objrow-1][j]);
      if (j>0) fprintf(f," X[%3ld]",j);
    }
    fprintf(f,"\n");
  }

  switch (lp->LPS){
  case dd_Optimal:
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

  case dd_Inconsistent:
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

  case dd_DualInconsistent: case dd_StrucDualInconsistent:
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
  fprintf(f,"* number of pivot operations = %ld (ph0 = %ld, ph1 = %ld, ph2 = %ld, ph3 = %ld, ph4 = %ld)\n",lp->total_pivots,lp->pivots[0],lp->pivots[1],lp->pivots[2],lp->pivots[3],lp->pivots[4]);
  dd_WriteLPTimes(f, lp);
_L99:;
}

dd_LPPtr dd_CreateLP_H_Redundancy(dd_MatrixPtr M, dd_rowrange itest)
{
  dd_rowrange m, i, irev, linc;
  dd_colrange d, j;
  dd_LPPtr lp;
  dd_boolean localdebug=dd_FALSE;

  linc=set_card(M->linset);
  m=M->rowsize+1+linc; 
     /* We represent each equation by two inequalities.
        This is not the best way but makes the code simple. */
  d=M->colsize;
  
  lp=dd_CreateLPData(M->objective, M->numbtype, m, d);
  lp->Homogeneous = dd_TRUE;
  lp->objective = dd_LPmin;
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
      if (localdebug) fprintf(stderr,"equality row %ld generates the reverse row %ld.\n",i,irev);
    }
    for (j = 1; j <= M->colsize; j++) {
      dd_set(lp->A[i-1][j-1],M->matrix[i-1][j-1]);
      if (j==1 && i<M->rowsize && dd_Nonzero(M->matrix[i-1][j-1])) lp->Homogeneous = dd_FALSE;
    }  /*of j*/
  }  /*of i*/
  for (j = 1; j <= M->colsize; j++) {
    dd_set(lp->A[m-1][j-1],M->matrix[itest-1][j-1]);
      /* objective is to violate the inequality in question.  */
  }  /*of j*/
  dd_add(lp->A[itest-1][0],lp->A[itest-1][0],dd_one); /* relax the original inequality by one */

  return lp;
}

dd_LPPtr dd_CreateLP_V_Redundancy(dd_MatrixPtr M, dd_rowrange itest)
{
  dd_rowrange m, i, irev, linc;
  dd_colrange d, j;
  dd_LPPtr lp;
  dd_boolean localdebug=dd_FALSE;

  linc=set_card(M->linset);
  m=M->rowsize+1+linc; 
     /* We represent each equation by two inequalities.
        This is not the best way but makes the code simple. */
  d=(M->colsize)+1;  
     /* One more column.  This is different from the H-reprentation case */
  
/* The below must be modified for V-representation!!!  */

  lp=dd_CreateLPData(M->objective, M->numbtype, m, d);
  lp->Homogeneous = dd_FALSE;
  lp->objective = dd_LPmin;
  lp->eqnumber=linc;  /* this records the number of equations */

  irev=M->rowsize; /* the first row of the linc reversed inequalities. */
  for (i = 1; i <= M->rowsize; i++) {
    if (i==itest){
      dd_set(lp->A[i-1][0],dd_one); /* this is to make the LP bounded, ie. the min >= -1 */
    } else {
      dd_set(lp->A[i-1][0],dd_purezero);  /* It is almost completely degerate LP */
    }
    if (set_member(i, M->linset)) {
      irev=irev+1;
      set_addelem(lp->equalityset,i);    /* it is equality. */
            /* the reversed row irev is not in the equality set. */
      for (j = 2; j <= (M->colsize)+1; j++) {
        dd_neg(lp->A[irev-1][j-1],M->matrix[i-1][j-2]);
      }  /*of j*/
      if (localdebug) fprintf(stderr,"equality row %ld generates the reverse row %ld.\n",i,irev);
    }
    for (j = 2; j <= (M->colsize)+1; j++) {
      dd_set(lp->A[i-1][j-1],M->matrix[i-1][j-2]);
    }  /*of j*/
  }  /*of i*/
  for (j = 2; j <= (M->colsize)+1; j++) {
    dd_set(lp->A[m-1][j-1],M->matrix[itest-1][j-2]);
      /* objective is to violate the inequality in question.  */
  }  /*of j*/
  dd_set(lp->A[m-1][0],dd_purezero);   /* the constant term for the objective is zero */

  if (localdebug) dd_WriteLP(stdout, lp);

  return lp;
}

dd_boolean dd_Redundant(dd_MatrixPtr M, dd_rowrange itest, dd_Arow certificate, dd_ErrorType *error)  
  /* 092 */
{
  /* Checks whether the row itest is redundant for the representation.
     All linearity rows are not checked and considered NONredundant. 
     This code works for both H- and V-representations.  A certificate is
     given in the case of non-redundancy, showing a solution x violating only the itest
     inequality for H-representation, a hyperplane RHS and normal (x_0, x) that
     separates the itest from the rest.  More explicitly, the LP to be setup is

     H-representation
       f* = minimize  
         b_itest     + A_itest x
       subject to
         b_itest + 1 + A_itest x     >= 0 (relaxed inequality to make an LP bounded)
         b_{I-itest} + A_{I-itest} x >= 0 (all inequalities except for itest)
         b_L         + A_L x = 0.  (linearity)

     V-representation (=separation problem)
       f* = minimize  
         b_itest x_0     + A_itest x
       subject to
         b_itest x_0     + A_itest x     >= -1 (to make an LP bounded)
         b_{I-itest} x_0 + A_{I-itest} x >=  0 (all nonlinearity generators except for itest in one side)
         b_L x_0         + A_L x = 0.  (linearity generators)
    
    Here, the input matrix is considered as (b, A), i.e. b corresponds to the first column of input
    and the row indices of input is partitioned into I and L where L is the set of linearity.
    In both cases, the itest data is nonredundant if and only if the optimal value f* is negative.
    The certificate has dimension one more for V-representation case.
  */

  dd_colrange j;
  dd_LPPtr lp;
  dd_LPSolutionPtr lps;
  dd_ErrorType err=dd_NoError;
  dd_boolean answer=dd_FALSE,localdebug=dd_FALSE;

  *error=dd_NoError;
  if (set_member(itest, M->linset)){
    if (localdebug) printf("The %ld th row is linearity and redundancy checking is skipped.\n",itest);
    goto _L99;
  }
  
  /* Create an LP data for redundancy checking */
  if (M->representation==dd_Generator){
    lp=dd_CreateLP_V_Redundancy(M, itest);
  } else {
    lp=dd_CreateLP_H_Redundancy(M, itest);
  }

  dd_LPSolve(lp,dd_DualSimplex,&err);
  if (err!=dd_NoError){
    *error=err;
    goto _L999;
  } else {
    lps=dd_CopyLPSolution(lp);

    for (j=0; j<lps->d; j++) {
      dd_set(certificate[j], lps->sol[j]);
    }

    if (dd_Negative(lps->optvalue)){
      answer=dd_FALSE;
      if (localdebug) fprintf(stderr,"==> %ld th inequality is nonredundant.\n",itest);
    } else {
      answer=dd_TRUE;
      if (localdebug) fprintf(stderr,"==> %ld th inequality is redundant.\n",itest);
    }
    dd_FreeLPSolution(lps);
  }
  _L999:
  dd_FreeLPData(lp);
_L99:
  return answer;
}

dd_rowset dd_RedundantRows(dd_MatrixPtr M, dd_ErrorType *error)  /* 092 */
{
  dd_rowrange i,m;
  dd_colrange d;
  dd_rowset redset;
  dd_MatrixPtr Mcopy;
  dd_Arow cvec; /* certificate */  
  dd_boolean localdebug=dd_FALSE;

  m=M->rowsize;
  if (M->representation==dd_Generator){
    d=(M->colsize)+1;
  } else {
    d=M->colsize;
  }
  Mcopy=dd_MatrixCopy(M);
  dd_InitializeArow(d,&cvec); 
  set_initialize(&redset, m);
  for (i=m; i>=1; i--) {
    if (dd_Redundant(Mcopy, i, cvec, error)) {
      if (localdebug) printf("dd_RedundantRows: the row %ld is redundant.\n", i);
      set_addelem(redset, i);
      dd_MatrixRowRemove(&Mcopy, i);
    } else {
      if (localdebug) printf("dd_RedundantRows: the row %ld is essential.\n", i);
    }
    if (*error!=dd_NoError) goto _L99;
  }
_L99:
  dd_FreeMatrix(Mcopy);
  dd_FreeArow(d, cvec);
  return redset;
}

dd_rowset dd_RedundantRowsViaShooting(dd_MatrixPtr M, dd_ErrorType *error)  /* 092 */
{
  /* 
     For H-representation only and not quite reliable,
     especially when floating-point arithmetic is used.
     Use the ordinary (slower) method dd_RedundantRows.
  */

  dd_rowrange i,m, ired, irow=0;
  dd_colrange j,k,d;
  dd_rowset redset;
  dd_rowindex rowflag; 
    /* ith comp is negative if the ith inequality (i-1 st row) is redundant.
                   zero     if it is not decided.
                   k > 0    if it is nonredundant and assigned to the (k-1)th row of M1.
    */
  dd_MatrixPtr M1;
  dd_Arow shootdir, cvec=NULL;
  dd_LPPtr lp0, lp;
  dd_LPSolutionPtr lps; 
  dd_ErrorType err;
  dd_LPSolverType solver=dd_DualSimplex; 
  dd_boolean localdebug=dd_FALSE;

  m=M->rowsize;
  d=M->colsize;
  M1=dd_CreateMatrix(m,d);
  M1->rowsize=0;  /* cheat the rowsize so that smaller matrix can be stored */
  set_initialize(&redset, m);
  dd_InitializeArow(d, &shootdir);
  dd_InitializeArow(d, &cvec);

  rowflag=(long *)calloc(m+1, sizeof(long)); 

  /* First find some (likely) nonredundant inequalities by Interior Point Find. */
  lp0=dd_Matrix2LP(M, &err);
  lp=dd_MakeLPforInteriorFinding(lp0);
  dd_FreeLPData(lp0); 
  dd_LPSolve(lp, solver, &err);  /* Solve the LP */
  lps=dd_CopyLPSolution(lp);

  if (dd_Positive(lps->optvalue)){
    /* An interior point is found.  Use rayshooting to find some nonredundant
       inequalities. */
    for (j=1; j<d; j++){
      for (k=1; k<=d; k++) dd_set(shootdir[k-1], dd_purezero);
      dd_set(shootdir[j], dd_one);  /* j-th unit vector */
      ired=dd_RayShooting(M, lps->sol, shootdir);
      if (localdebug) printf("nonredundant row %3ld found by shooting.\n", ired);
      if (ired>0 && rowflag[ired]<=0) {
        irow++;
        rowflag[ired]=irow;
        for (k=1; k<=d; k++) dd_set(M1->matrix[irow-1][k-1], M->matrix[ired-1][k-1]); 
      }
        
      dd_neg(shootdir[j], dd_one);  /* negative of the j-th unit vector */
      ired=dd_RayShooting(M, lps->sol, shootdir);
      if (localdebug) printf("nonredundant row %3ld found by shooting.\n", ired);
      if (ired>0 && rowflag[ired]<=0) {
        irow++;
        rowflag[ired]=irow;
        for (k=1; k<=d; k++) dd_set(M1->matrix[irow-1][k-1], M->matrix[ired-1][k-1]); 
      }
    }

    M1->rowsize=irow;
    if (localdebug) {
      printf("The initial nonredundant set is:");
      for (i=1; i<=m; i++) if (rowflag[i]>0) printf(" %ld", i);
      printf("\n");
    }
    
    i=1;
    while(i<=m){
      if (rowflag[i]==0){ /* the ith inequality is not yet checked */
        if (localdebug) fprintf(stderr, "Checking redundancy of %ld th inequality\n", i);
        irow++;  M1->rowsize=irow;
        for (k=1; k<=d; k++) dd_set(M1->matrix[irow-1][k-1], M->matrix[i-1][k-1]);
        if (!dd_Redundant(M1, irow, cvec, &err)){
          for (k=1; k<=d; k++) dd_sub(shootdir[k-1], cvec[k-1], lps->sol[k-1]); 
          ired=dd_RayShooting(M, lps->sol, shootdir);
          rowflag[ired]=irow;
          for (k=1; k<=d; k++) dd_set(M1->matrix[irow-1][k-1], M->matrix[ired-1][k-1]);
          if (localdebug) {
            fprintf(stderr, "The %ld th inequality is nonredundant for the subsystem\n", i);
            fprintf(stderr, "The nonredundancy of %ld th inequality is found by shooting.\n", ired);
          }
        } else {
          if (localdebug) fprintf(stderr, "The %ld th inequality is redundant for the subsystem and thus for the whole.\n", i);
          rowflag[i]=-1;
          set_addelem(redset, i);
          i++;
        }
      } else {
        i++;
      }
    } /* endwhile */
  } else {
    /* No interior point is found.  Apply the standard LP technique.  */
    redset=dd_RedundantRows(M, error);
  }

  dd_FreeLPData(lp);
  dd_FreeLPSolution(lps);

  M1->rowsize=m; M1->colsize=d;  /* recover the original sizes */
  dd_FreeMatrix(M1);
  dd_FreeArow(d, shootdir);
  dd_FreeArow(d, cvec);
  free(rowflag);
  return redset;
}

dd_SetFamilyPtr dd_Matrix2Adjacency(dd_MatrixPtr M, dd_ErrorType *error)  /* 093 */
{
  /* This is to generate the graph of a polyheron representation given by M using LPs.
     Since it does not use the representation conversion, it should work for a large
     scale problem.
  */
  dd_rowrange i,m;
  dd_colrange d;
  dd_rowset redset;
  dd_MatrixPtr Mcopy;
  dd_SetFamilyPtr F=NULL;

  m=M->rowsize;
  d=M->colsize;
  if (m<=0 ||d<=0) {
    *error=dd_EmptyRepresentation;
    goto _L999;
  }
  Mcopy=dd_MatrixCopy(M);
  F=dd_CreateSetFamily(m, m);
  for (i=1; i<=m; i++) {
    if (!set_member(i, M->linset)){
      set_addelem(Mcopy->linset, i);
      redset=dd_RedundantRows(Mcopy, error);  /* redset should contain all nonadjacent ones */
      set_uni(redset, redset, Mcopy->linset); /* all linearity elements should be nonadjacent */
      set_compl(F->set[i-1], redset); /* set the adjacency list of vertex i */
      set_delelem(Mcopy->linset, i);
      set_free(redset);
      if (*error!=dd_NoError) goto _L99;
    }
  }
_L99:
  dd_FreeMatrix(Mcopy);
_L999:
  return F;
}


dd_boolean dd_ImplicitLinearity(dd_MatrixPtr M, dd_rowrange itest, dd_Arow certificate, dd_ErrorType *error)  
  /* 092 */
{
  /* Checks whether the row itest is implicit linearity for the representation.
     All linearity rows are not checked and considered non implicit linearity (dd_FALSE). 
     This code works for both H- and V-representations.  A certificate is
     given in the case of dd_FALSE, showing a feasible solution x satisfying the itest
     strict inequality for H-representation, a hyperplane RHS and normal (x_0, x) that
     separates the itest from the rest.  More explicitly, the LP to be setup is
     the same thing as redundancy case but with maximization:

     H-representation
       f* = maximize  
         b_itest     + A_itest x
       subject to
         b_itest + 1 + A_itest x     >= 0 (relaxed inequality. This is not necessary but kept for simplicity of the code)
         b_{I-itest} + A_{I-itest} x >= 0 (all inequalities except for itest)
         b_L         + A_L x = 0.  (linearity)

     V-representation (=separation problem)
       f* = maximize  
         b_itest x_0     + A_itest x
       subject to
         b_itest x_0     + A_itest x     >= -1 (again, this is not necessary but kept for simplicity.)
         b_{I-itest} x_0 + A_{I-itest} x >=  0 (all nonlinearity generators except for itest in one side)
         b_L x_0         + A_L x = 0.  (linearity generators)
    
    Here, the input matrix is considered as (b, A), i.e. b corresponds to the first column of input
    and the row indices of input is partitioned into I and L where L is the set of linearity.
    In both cases, the itest data is implicit linearity if and only if the optimal value f* is nonpositive.
    The certificate has dimension one more for V-representation case.
  */

  dd_colrange j;
  dd_LPPtr lp;
  dd_LPSolutionPtr lps;
  dd_ErrorType err=dd_NoError;
  dd_boolean answer=dd_FALSE,localdebug=dd_FALSE;

  *error=dd_NoError;
  if (set_member(itest, M->linset)){
    if (localdebug) printf("The %ld th row is linearity and redundancy checking is skipped.\n",itest);
    goto _L99;
  }
  
  /* Create an LP data for redundancy checking */
  if (M->representation==dd_Generator){
    lp=dd_CreateLP_V_Redundancy(M, itest);
  } else {
    lp=dd_CreateLP_H_Redundancy(M, itest);
  }

  lp->objective = dd_LPmax;  /* the lp->objective is set by CreateLP* to LPmin */
  dd_LPSolve(lp,dd_DualSimplex,&err);
  if (err!=dd_NoError){
    *error=err;
    goto _L999;
  } else {
    lps=dd_CopyLPSolution(lp);

    for (j=0; j<lps->d; j++) {
      dd_set(certificate[j], lps->sol[j]);
    }

    if (lps->LPS==dd_Optimal && dd_EqualToZero(lps->optvalue)){
      answer=dd_TRUE;
      if (localdebug) fprintf(stderr,"==> %ld th data is an implicit linearity.\n",itest);
    } else {
      answer=dd_FALSE;
      if (localdebug) fprintf(stderr,"==> %ld th data is not an implicit linearity.\n",itest);
    }
    dd_FreeLPSolution(lps);
  }
  _L999:
  dd_FreeLPData(lp);
_L99:
  return answer;
}


dd_rowset dd_ImplicitLinearityRows(dd_MatrixPtr M, dd_ErrorType *error)  /* 092 */
{
  dd_rowrange i,m;
  dd_colrange d;
  dd_rowset eqset;
  dd_Arow cvec; /* certificate */

  if (M->representation==dd_Generator){
    d=(M->colsize)+1;
  } else {
    d=M->colsize;
  }

  dd_InitializeArow(d,&cvec);
  m=M->rowsize;
  set_initialize(&eqset, m);
  for (i=m; i>=1; i--) {
    if (dd_ImplicitLinearity(M, i, cvec, error)) {
      set_addelem(eqset, i);
    }
    if (*error!=dd_NoError) goto _L99;
  }
_L99:
  dd_FreeArow(d, cvec);
  return eqset;
}

dd_rowrange dd_RayShooting(dd_MatrixPtr M, dd_Arow p, dd_Arow r)
{
/* 092, find the first inequality "hit" by a ray from an intpt.  */
  dd_rowrange imin=-1,i,m;
  dd_colrange j, d;
  dd_Arow vecmin, vec;
  mytype min,t1,t2,alpha, t1min;  
  dd_boolean started=dd_FALSE;
  dd_boolean localdebug=dd_FALSE;

  m=M->rowsize;
  d=M->colsize;
  if (!dd_Equal(dd_one, p[0])){
    fprintf(stderr, "Warning: RayShooting is called with a point with first coordinate not 1.\n");
    dd_set(p[0],dd_one);
  }
  if (!dd_EqualToZero(r[0])){
    fprintf(stderr, "Warning: RayShooting is called with a direction with first coordinate not 0.\n");
    dd_set(r[0],dd_purezero);
  }

  dd_init(alpha); dd_init(min); dd_init(t1); dd_init(t2); dd_init(t1min);
  dd_InitializeArow(d,&vecmin);
  dd_InitializeArow(d,&vec);

  for (i=1; i<=m; i++){
    dd_InnerProduct(t1, d, M->matrix[i-1], p);
    if (dd_Positive(t1)) {
      dd_InnerProduct(t2, d, M->matrix[i-1], r);
      dd_div(alpha, t2, t1);
      if (!started){
        imin=i;  dd_set(min, alpha);
        dd_set(t1min, t1);  /* store the denominator. */
        started=dd_TRUE;
        if (localdebug) {
          fprintf(stderr," Level 1: imin = %ld and min = ", imin);
          dd_WriteNumber(stderr, min);
          fprintf(stderr,"\n");
        }
      } else {
        if (dd_Smaller(alpha, min)){
          imin=i;  dd_set(min, alpha);
          dd_set(t1min, t1);  /* store the denominator. */
          if (localdebug) {
            fprintf(stderr," Level 2: imin = %ld and min = ", imin);
            dd_WriteNumber(stderr, min);
            fprintf(stderr,"\n");
          }
        } else {
          if (dd_Equal(alpha, min)) { /* tie break */
            for (j=1; j<= d; j++){
              dd_div(vecmin[j-1], M->matrix[imin-1][j-1], t1min);
              dd_div(vec[j-1], M->matrix[i-1][j-1], t1);
            }
            if (dd_LexSmaller(vec,vecmin, d)){
              imin=i;  dd_set(min, alpha);
              dd_set(t1min, t1);  /* store the denominator. */
              if (localdebug) {
                fprintf(stderr," Level 3: imin = %ld and min = ", imin);
                dd_WriteNumber(stderr, min);
                fprintf(stderr,"\n");
              }
            }
          }
        }
      }       
    }
  }

  dd_clear(alpha); dd_clear(min); dd_clear(t1); dd_clear(t2); dd_clear(t1min);
  dd_FreeArow(d, vecmin);
  dd_FreeArow(d, vec);
  return imin;
}

#ifdef GMPRATIONAL
void dd_BasisStatusMaximize(dd_rowrange m_size,dd_colrange d_size,
    dd_Amatrix A,dd_Bmatrix T,dd_rowset equalityset,
    dd_rowrange objrow,dd_colrange rhscol,ddf_LPStatusType LPS,
    mytype *optvalue,dd_Arow sol,dd_Arow dsol,ddf_colindex nbindex,
    ddf_rowrange re,ddf_colrange se,long *pivots, int *found, int *LPScorrect)
/*  This is just to check whether the status LPS of the basis given by 
nbindex with extra certificates se or re is correct.  It is done
by recomputing the basis inverse matrix T.  It does not solve the LP
when the status *LPS is undecided.  Thus the input is
m_size, d_size, A, equalityset, LPS, nbindex, re and se.
Other values will be recomputed from scratch.

The main purpose of the function is to verify the correctness
of the result of floating point computation with the GMP rational
arithmetics.
*/
{
  long pivots0,pivots1;
  dd_rowrange i,is;
  dd_colrange s,j;
  static dd_rowindex bflag;
  static long mlast=0;
  static dd_rowindex OrderVector;  /* the permutation vector to store a preordered row indeces */
  unsigned int rseed=1;
  mytype val;
  dd_colindex nbtemp;
  dd_LPStatusType ddlps;
  dd_boolean localdebug=dd_FALSE;


  dd_init(val);
  nbtemp=(long *) calloc(d_size+1,sizeof(long*));
  for (i=0; i<= 4; i++) pivots[i]=0;
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
  dd_ComputeRowOrderVector2(m_size,d_size,A,OrderVector,dd_MinIndex,rseed);

  pivots1=0;

  dd_ResetTableau(m_size,d_size,T,nbtemp,bflag,objrow,rhscol);

  if (localdebug){
     printf("\nnbindex:");
     for (j=1; j<=d_size; j++) printf(" %ld", nbindex[j]);
     printf("\n");
     printf("re = %ld,   se=%ld\n", re, se);
  }
  
  is=nbindex[se];
  if (localdebug) printf("se=%ld,  is=%ld\n", se, is);
  
  dd_FindLPBasis2(m_size,d_size,A,T,OrderVector, equalityset,nbindex,bflag,
      objrow,rhscol,&s,found,&pivots0);
      
  pivots[4]=pivots0;  /*GMP postopt pivots */

  if (!(*found)){
    if (localdebug) {
       printf("dd_BasisStatusMaximize: a specified basis DOES NOT exist.\n");
    }

       goto _L99;
     /* No speficied LP basis is found. */
  }

  if (localdebug) {
    printf("dd_BasisStatusMaximize: a specified basis exists.\n");
    dd_WriteTableau(stdout,m_size,d_size,A,T,nbindex,bflag);
  }

  /* Check whether a recomputed basis is of the type specified by LPS */
  *LPScorrect=dd_TRUE;
  switch (LPS){
     case dd_Optimal: 
       for (i=1; i<=m_size; i++) {
         if (i!=objrow && bflag[i]==-1) {  /* i is a basic variable */
            dd_TableauEntry(&val,m_size,d_size,A,T,i,rhscol);
            if (dd_Negative(val)) {
               if (localdebug) printf("RHS entry for %ld is negative\n", i);
               *LPScorrect=dd_FALSE;
               break;
            }
          } else if (bflag[i] >0) { /* i is nonbasic variable */
            dd_TableauEntry(&val,m_size,d_size,A,T,objrow,bflag[i]);
            if (dd_Positive(val)) {
               if (localdebug) printf("Reduced cost entry for %ld is positive\n", i);
               *LPScorrect=dd_FALSE;
               break;
            }
          }
       };
       break;
     case dd_Inconsistent: 
       for (j=1; j<=d_size; j++){
          dd_TableauEntry(&val,m_size,d_size,A,T,re,j);
          if (j==rhscol){
             if (dd_Nonnegative(val)){
               if (localdebug) printf("RHS entry for %ld is nonnegative\n", re);
               *LPScorrect=dd_FALSE;
               break;             
             }
           } else if (dd_Positive(val)){
               if (localdebug) printf("the row entry for(%ld, %ld) is positive\n", re, j);
               *LPScorrect=dd_FALSE;
               break;             
           }
       };
       break;
     case dd_DualInconsistent:
        for (i=1; i<=m_size; i++){
          dd_TableauEntry(&val,m_size,d_size,A,T,i,bflag[is]);
          if (i==objrow){
             if (dd_Nonpositive(val)){
               if (localdebug) printf("Reduced cost entry for %ld is nonpositive\n", bflag[is]);
               *LPScorrect=dd_FALSE;
               break;             
             }
           } else if (dd_Negative(val)){
               if (localdebug) printf("the column entry for(%ld, %ld) is positive\n", i, bflag[is]);
               *LPScorrect=dd_FALSE;
               break;             
           }
       };
       break;
;
     default: break;
  }

  ddlps=LPSf2LPS(LPS);

  dd_SetSolutions(m_size,d_size,A,T,
   objrow,rhscol,ddlps,optvalue,sol,dsol,nbindex,re,se);

  
_L99:
  dd_clear(val);
  free(nbtemp);
}

void dd_BasisStatusMinimize(dd_rowrange m_size,dd_colrange d_size,
    dd_Amatrix A,dd_Bmatrix T,dd_rowset equalityset,
    dd_rowrange objrow,dd_colrange rhscol,ddf_LPStatusType LPS,
    mytype *optvalue,dd_Arow sol,dd_Arow dsol,ddf_colindex nbindex,
    ddf_rowrange re,ddf_colrange se,long *pivots, int *found, int *LPScorrect)
{
   dd_colrange j;
   
   for (j=1; j<=d_size; j++) dd_neg(A[objrow-1][j-1],A[objrow-1][j-1]);
   dd_BasisStatusMaximize(m_size,d_size,A,T,equalityset, objrow,rhscol,
     LPS,optvalue,sol,dsol,nbindex,re,se,pivots,found,LPScorrect);
   dd_neg(*optvalue,*optvalue);
   for (j=1; j<=d_size; j++){
     dd_neg(dsol[j-1],dsol[j-1]);
     dd_neg(A[objrow-1][j-1],A[objrow-1][j-1]);
   }
}
#endif

/* end of cddlp.c */

