/* dplex.h: Header file for dplex.c 
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.8, March 13, 1999
*/

/* dplex.c : C-Implementation of the dual simplex method for
   solving an LP: max/min  A_(m-1).x subject to  x in P, where
   P= {x :  A_i.x >= 0, i=0,...,m-2, and  x_0=1}, and
   A_i is the i-th row of an m x n matrix A.
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#include <time.h>

#define dp_zero      1.0E-6    /*real zero */
#define dp_magfac    1000000   /*reciprocal of dp_zero*/
#define dp_wordlenmax  32
#define dp_linelenmax  255
#define dp_filenamelen 255
typedef long dp_rowrange;
typedef long dp_colrange;
typedef long *dp_rowindex;   
typedef long *dp_colindex;
typedef double dp_mynumber;
typedef double **dp_Amatrix;
typedef double *dp_Arow;
typedef double **dp_Bmatrix;

typedef char dp_FilenameType[dp_filenamelen];

typedef enum {
  dp_DimensionTooLarge, dp_LowColumnRank, dp_ImproperInputFormat, 
  dp_FileNotFound, 
  dp_LPLoadAmatrixError, 
  dp_None
} dp_ErrorType;

typedef enum {
  dp_Real, dp_Rational, dp_Integer, dp_Unknown
} dp_NumberType;

typedef enum {
  dp_LPmax, dp_LPmin
} dp_LPConversionType;

typedef enum {
  dp_CrissCross, dp_DualSimplex
} dp_LPSolverType;

typedef enum {
  dp_LPSundecided, dp_Optimal, dp_Inconsistent, dp_DualInconsistent,
  dp_StrucInconsistent, dp_StrucDualInconsistent,
  dp_Unbounded, dp_DualUnbounded
} dp_LPStatusType;


typedef struct {
  dp_FilenameType filename;
  dp_LPConversionType conv;
  dp_LPSolverType solver; 
  dp_rowrange m;
  dp_colrange d;
  dp_NumberType number;

  dp_LPStatusType LPS;  /* the current solution status */
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
} dp_LPSolutionType ;

typedef dp_LPSolutionType *dp_LPSolutionPtr;

typedef struct lpdata *dp_LPPtr;

dp_LPPtr dp_LPInput(FILE **f, dp_ErrorType *err);  

dp_LPPtr dp_LPLoad(dp_LPConversionType,
   dp_NumberType, dp_rowrange m, dp_colrange d, dp_Amatrix, 
   dp_rowrange OBJrow, dp_colrange RHScol, dp_ErrorType *err);  
   /* 
      Load an LP safely.  This creates a copy of LP data,
      and returns a pointer to the LPDataType.  
   */
dp_LPPtr dp_LPDirectLoad(dp_LPConversionType,
   dp_NumberType, dp_rowrange m, dp_colrange d, dp_Amatrix*, 
   dp_rowrange OBJrow, dp_colrange RHScol, dp_ErrorType *err);  
   /* 
      Load an LP quickly.  This creates a direct links to A
      and a copy of other data, and returns a pointer to the LPDataType.
      This deletes the pointer A so that the loaded LPdata won't be
      affected by the user.  Use this function only when saving memory
      or time is extremely important.
   */
void dp_LPSolve(dp_LPPtr, dp_ErrorType *);
dp_LPPtr dp_MakeLPforInteriorFinding(dp_LPPtr);  
dp_LPSolutionPtr dp_LPSolutionLoad(dp_LPPtr lp);

int dp_LPReverseRow(dp_LPPtr, dp_rowrange);
    /* reverse the i-th row (1 <= i <= no. of rows) */
int dp_LPReplaceRow(dp_LPPtr, dp_rowrange, dp_Arow);
    /* replace the i-th row (1 <= i <= no. of rows) */
dp_LPPtr dp_LPCopy(dp_LPPtr);
    /* copy an LP data */
dp_Arow dp_LPCopyRow(dp_LPPtr, dp_rowrange);
    /* copy the i-th row (1 <= i <= no. of rows) */

void dp_FreeLPData(dp_LPPtr*);
void dp_FreeLPSolution(dp_LPSolutionPtr*);

void dp_InitializeArow(dp_colrange,dp_Arow *);
void dp_InitializeAmatrix(dp_rowrange,dp_colrange,dp_Amatrix *);

void dp_WriteLPResult(FILE *, dp_LPPtr, dp_ErrorType);
void dp_WriteErrorMessages(FILE *, dp_ErrorType);
void dp_WriteReal(FILE *, double);

int dp_Nonnegative(double);
int dp_Nonpositive(double);
int dp_Positive(double);
int dp_Negative(double);
int dp_Zero(double);
int dp_Nonzero(double);


/* end of dplex.h */
