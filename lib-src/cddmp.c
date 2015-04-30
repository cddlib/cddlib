/* cddmp.c       (cddlib arithmetic operations using gmp)
   written by Komei Fukuda, fukuda@math.ethz.ch
   Version 0.94h, April 30, 2015
*/
/* This program is free software; you can redistribute it and/or modify
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

#include "setoper.h"  /* set operation library header (Ver. March 16,1995 or later) */
#include "cdd.h"

void dd_set_global_constants()
{
 dd_init(dd_zero);
 dd_init(dd_minuszero);
 dd_init(dd_one);
 dd_init(dd_minusone);
 dd_init(dd_purezero);
  
 time(&dd_statStartTime); /* cddlib starting time */
 dd_statBApivots=0;  /* basis finding pivots */
 dd_statCCpivots=0;  /* criss-cross pivots */
 dd_statDS1pivots=0; /* phase 1 pivots */
 dd_statDS2pivots=0; /* phase 2 pivots */
 dd_statACpivots=0;  /* anticycling (cc) pivots */

 dd_choiceLPSolverDefault=dd_DualSimplex;  /* Default LP solver Algorithm */
 dd_choiceRedcheckAlgorithm=dd_DualSimplex;  /* Redundancy Checking Algorithm */
 dd_choiceLexicoPivotQ=dd_TRUE;    /* whether to use the lexicographic pivot */
 
#if defined GMPRATIONAL
 dd_statBSpivots=0;  /* basis status checking pivots */
 mpq_set_ui(dd_zero,0U,1U);
 mpq_set_ui(dd_purezero,0U,1U);
 mpq_set_ui(dd_one,1U,1U);
 mpq_set_si(dd_minusone,-1L,1U);
 ddf_set_global_constants();
#elif defined GMPFLOAT
 mpf_set_d(dd_zero,dd_almostzero);
 mpf_set_ui(dd_purezero,0U);
 mpf_set_ui(dd_one,1U);
 mpf_set_si(dd_minusone,-1L,1U);
#else
 dd_zero[0]= dd_almostzero;  /*real zero */
 dd_purezero[0]= 0.0;
 dd_one[0]= 1L;
 dd_minusone[0]= -1L;
#endif
 dd_neg(dd_minuszero,dd_zero);
}

void dd_free_global_constants()
{
 dd_clear(dd_zero);
 dd_clear(dd_minuszero);
 dd_clear(dd_one);
 dd_clear(dd_minusone);
 dd_clear(dd_purezero);
  
 time(&dd_statStartTime); /* cddlib starting time */
 dd_statBApivots=0;  /* basis finding pivots */
 dd_statCCpivots=0;  /* criss-cross pivots */
 dd_statDS1pivots=0; /* phase 1 pivots */
 dd_statDS2pivots=0; /* phase 2 pivots */
 dd_statACpivots=0;  /* anticycling (cc) pivots */

 dd_choiceLPSolverDefault=dd_DualSimplex;  /* Default LP solver Algorithm */
 dd_choiceRedcheckAlgorithm=dd_DualSimplex;  /* Redundancy Checking Algorithm */
 dd_choiceLexicoPivotQ=dd_TRUE;    /* whether to use the lexicographic pivot */
 
#if defined GMPRATIONAL
 dd_statBSpivots=0;  /* basis status checking pivots */
 ddf_free_global_constants();
#endif
}


#if defined GMPRATIONAL
void ddd_mpq_set_si(mytype a,signed long b)
{
  mpz_t nz, dz;

  mpz_init(nz); mpz_init(dz);

  mpz_set_si(nz, b);
  mpz_set_ui(dz, 1U);
  mpq_set_num(a, nz);
  mpq_set_den(a, dz);
  mpz_clear(nz);  mpz_clear(dz);
}
#endif

#if defined dd_CDOUBLE
void ddd_init(mytype a)   
{
  a[0]=0L;
}
  
void ddd_clear(mytype a)
{
  /* a[0]=0L;  */
}

void ddd_set(mytype a,mytype b)
{
  a[0]=b[0];
}

void ddd_set_d(mytype a,double b)
{
  a[0]=b;
}

void ddd_set_si(mytype a,signed long b)
{
  a[0]=(double)b;
}

void ddd_set_si2(mytype a,signed long b, unsigned long c)
{
  a[0]=(double)b/(double)c;
}

void ddd_add(mytype a,mytype b,mytype c)
{
  a[0]=b[0]+c[0];
}

void ddd_sub(mytype a,mytype b,mytype c)
{
  a[0]=b[0]-c[0];
}

void ddd_mul(mytype a,mytype b,mytype c)
{
  a[0]=b[0]*c[0];
}

void ddd_div(mytype a,mytype b,mytype c)
{
  a[0]=b[0]/c[0];
}

void ddd_neg(mytype a,mytype b)
{
  a[0]=-b[0];
}

void ddd_inv(mytype a,mytype b)
{
  a[0]=1/b[0];
}

int ddd_cmp(mytype a,mytype b)
{
  if (a[0]-b[0]>0) return 1;
  else if (a[0]-b[0]>=0) return 0;
  else return -1;
}

int ddd_sgn(mytype a)
{
  if (a[0]>0) return 1;
  else if (a[0]>=0) return 0;
  else return -1;
}

double ddd_get_d(mytype a)
{
  return a[0];
}
#endif

/* end of  cddmp.h  */
