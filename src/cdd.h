/* cddlib.h: Header file for cddlib.c 
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.90, May 18, 2000
*/

/* cddlib.c : C-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#ifndef  __CDDLIB_H
#define  __CDDLIB_H
#endif  /* __CDDLIB_H */

#ifndef  __CDDMP_H
#include "cddmp.h"
#endif  /* __CDDMP_H */

#ifndef  __CDDTYPES_H
#include "cddtypes.h"
#endif  /* __CDDTYPES_H */


/* GLOBAL CONSTANTS (to be set by dd_set_global_constants() */
extern mytype dd_zero;
extern mytype dd_one;
extern mytype dd_purezero;
extern mytype dd_minuszero;
   /* to be used to avoid creating temporary spaces for mytype */
#define dd_almostzero  1.0E-6

/* ---------- FUNCTIONS MEANT TO BE PUBLIC ---------- */

/* basic matrix manipulations */
void dd_InitializeArow(dd_colrange,dd_Arow *);
void dd_InitializeAmatrix(dd_rowrange,dd_colrange,dd_Amatrix *);
void dd_InitializeBmatrix(dd_colrange, dd_Bmatrix *);
dd_SetFamilyPtr dd_CreateSetFamily(dd_bigrange,dd_bigrange);
void dd_FreeSetFamily(dd_SetFamilyPtr);
dd_MatrixPtr dd_CreateMatrix(dd_rowrange,dd_colrange);
void dd_FreeAmatrix(dd_rowrange,dd_colrange,dd_Amatrix);
void dd_FreeArow(dd_colrange, dd_Arow);
void dd_FreeBmatrix(dd_colrange,dd_Bmatrix);
void dd_FreeDDMemory(dd_PolyhedraPtr);
void dd_FreePolyhedra(dd_PolyhedraPtr);
void dd_FreeMatrix(dd_MatrixPtr);
void dd_SetToIdentity(dd_colrange, dd_Bmatrix);

/* sign recognitions */
boolean dd_Nonnegative(mytype);
boolean dd_Nonpositive(mytype);
boolean dd_Positive(mytype);
boolean dd_Negative(mytype);
boolean dd_EqualToZero(mytype);
boolean dd_Nonzero(mytype);
boolean dd_Equal(mytype,mytype);
boolean dd_Larger(mytype,mytype);
boolean dd_Smaller(mytype,mytype);
void dd_abs(mytype, mytype);
void dd_lincomb(mytype, mytype, mytype, mytype, mytype);


/* major cddlib operations */
dd_MatrixPtr dd_CopyInequalities(dd_PolyhedraPtr);
dd_MatrixPtr dd_CopyGenerators(dd_PolyhedraPtr);
dd_SetFamilyPtr dd_CopyIncidence(dd_PolyhedraPtr);
dd_SetFamilyPtr dd_CopyAdjacency(dd_PolyhedraPtr);
dd_SetFamilyPtr dd_CopyInputIncidence(dd_PolyhedraPtr);
dd_SetFamilyPtr dd_CopyInputAdjacency(dd_PolyhedraPtr);
boolean dd_DoubleDescription(dd_PolyhedraPtr, dd_ErrorType*);
boolean dd_DDFile2File(char *ifile, char *ofile, dd_ErrorType *err);
boolean dd_DDAddInequalities(dd_PolyhedraPtr, dd_MatrixPtr, dd_ErrorType*);
dd_MatrixPtr dd_PolyFile2Matrix (FILE *f, dd_ErrorType *);
dd_PolyhedraPtr dd_Matrix2Poly(dd_MatrixPtr, dd_ErrorType *);

/* input/output */
void dd_SetInputFile(FILE **f,dd_DataFileType inputfile, dd_ErrorType *);
void dd_SetWriteFileName(dd_DataFileType, dd_DataFileType, char, dd_RepresentationType);

void dd_WriteAmatrix(FILE *, dd_Amatrix, dd_rowrange, dd_colrange);
void dd_WriteBmatrix(FILE *, dd_colrange, dd_Bmatrix T);
void dd_WriteMatrix(FILE *, dd_MatrixPtr);
void dd_MatrixIntegerFilter(dd_MatrixPtr);
void dd_WriteReal(FILE *, mytype);
void dd_WriteNumber(FILE *f, mytype x); 
    /* write a number depending on the arithmetic used.  */
void dd_WritePolyFile(FILE *, dd_PolyhedraPtr);
void dd_WriteRunningMode(FILE *, dd_ConePtr);
void dd_WriteErrorMessages(FILE *, dd_ErrorType);
void dd_WriteSetFamily(FILE *, dd_SetFamilyPtr);
void dd_WriteSetFamilyCompressed(FILE *, dd_SetFamilyPtr);
void dd_WriteProgramDescription(FILE *);
void dd_WriteDDTimes(FILE *, dd_PolyhedraPtr);
void dd_WriteTimes(FILE *, time_t, time_t);


/* ---------- FUNCTIONS MEANT TO BE NON-PUBLIC ---------- */

void fread_rational_value (FILE *f, mytype value);
void AddNewHalfspace1(dd_ConePtr, dd_rowrange);
void AddNewHalfspace2(dd_ConePtr, dd_rowrange);
void AddRay(dd_ConePtr, mytype *);
void AddArtificialRay(dd_ConePtr);
void AValue(mytype*,dd_colrange, dd_Amatrix, mytype *, dd_rowrange);
void CheckAdjacency(dd_ConePtr,
    dd_RayPtr*, dd_RayPtr*, boolean *);
void CheckEquality(dd_colrange, dd_RayPtr *, dd_RayPtr *, boolean *);
void ComputeRowOrderVector(dd_ConePtr);
void ConditionalAddEdge(dd_ConePtr,dd_RayPtr, dd_RayPtr, dd_RayPtr);
void CopyArow(mytype *, mytype *, dd_colrange);
void CopyAmatrix(mytype **, mytype **, dd_rowrange, dd_colrange);
void CopyBmatrix(dd_colrange, dd_Bmatrix T, dd_Bmatrix TCOPY);
void CopyRay(mytype *, dd_colrange, dd_RayPtr,
   dd_RepresentationType, dd_colindex);
void CreateInitialEdges(dd_ConePtr);
void CreateNewRay(dd_ConePtr, dd_RayPtr, dd_RayPtr, dd_rowrange);
void Eliminate(dd_ConePtr, dd_RayPtr*);
void EvaluateARay1(dd_rowrange, dd_ConePtr);
void EvaluateARay2(dd_rowrange, dd_ConePtr);
void FeasibilityIndices(long *, long *, dd_rowrange, dd_ConePtr);
void FindBasis(dd_ConePtr, long *rank);
void FindInitialRays(dd_ConePtr, boolean *);
void ColumnReduce(dd_ConePtr);
void GaussianColumnPivot(dd_rowrange, dd_colrange, dd_Amatrix, dd_Bmatrix,  dd_rowrange, dd_colrange);
boolean LexSmaller(mytype *, mytype *, long);
boolean LexLarger(mytype *, mytype *, long);
void Normalize(dd_colrange, mytype *);
void MatrixIntegerFilter(dd_MatrixPtr);
void ProcessCommandLine(FILE*,dd_MatrixPtr, char *);
void SelectNextHalfspace(dd_ConePtr, dd_rowset, dd_rowrange *);
void SelectPivot2(dd_rowrange,dd_colrange,dd_Amatrix,
dd_Bmatrix,dd_RowOrderType,dd_rowindex, dd_rowset,dd_rowrange,dd_rowset,
dd_colset,dd_rowrange *,dd_colrange *,boolean *);
void SelectPreorderedNext(dd_ConePtr, dd_rowset, dd_rowrange *);
void SetInequalitySets(dd_ConePtr);
void SnapToInteger(mytype, mytype);
void StoreRay1(dd_ConePtr, mytype *, boolean *);
void StoreRay2(dd_ConePtr, mytype *, boolean *, boolean *);
void TableauEntry(mytype *, dd_rowrange, dd_colrange, dd_Amatrix, dd_Bmatrix T, dd_rowrange, dd_colrange);
void UpdateEdges(dd_ConePtr, dd_RayPtr, dd_RayPtr);
void UpdateRowOrderVector(dd_ConePtr, dd_rowset PriorityRows);
void dd_WriteIncidence(FILE *, dd_PolyhedraPtr);
void dd_WriteAdjacency(FILE *, dd_PolyhedraPtr);
void dd_WriteInputAdjacency(FILE *, dd_PolyhedraPtr);
void dd_WriteInputIncidence(FILE *, dd_PolyhedraPtr);
void WriteArow(FILE *f, dd_colrange d_origsize, dd_Arow a);
void WriteRay(FILE *, dd_colrange, dd_RayPtr,
   dd_RepresentationType, dd_colindex);
void ZeroIndexSet(dd_rowrange, dd_colrange, dd_Amatrix, mytype *, dd_rowset);

/* New functions to handle data loading, NON-PUBLIC */
dd_NumberType GetNumberType(char *);
dd_ConePtr ConeDataLoad(dd_PolyhedraPtr);
dd_PolyhedraPtr CreatePolyhedraData(dd_rowrange, dd_colrange);
boolean InitializeConeData(dd_rowrange, dd_colrange, dd_ConePtr*);
void AddInequalities(dd_PolyhedraPtr, dd_MatrixPtr);


/* functions and types for LP solving */

dd_LPPtr dd_Matrix2LP(dd_MatrixPtr, dd_ErrorType *);
  /* a new way to load a matrix to create an LP object. */

boolean dd_LPSolve(dd_LPPtr,dd_LPSolverType,dd_ErrorType *);
dd_LPPtr dd_MakeLPforInteriorFinding(dd_LPPtr);  
dd_LPSolutionPtr dd_LPSolutionLoad(dd_LPPtr lp);

int dd_LPReverseRow(dd_LPPtr, dd_rowrange);
    /* reverse the i-th row (1 <= i <= no. of rows) */
int dd_LPReplaceRow(dd_LPPtr, dd_rowrange, dd_Arow);
    /* replace the i-th row (1 <= i <= no. of rows) */
dd_Arow dd_LPCopyRow(dd_LPPtr, dd_rowrange);
    /* copy the i-th row (1 <= i <= no. of rows) */

void dd_FreeLPData(dd_LPPtr);
void dd_FreeLPSolution(dd_LPSolutionPtr);

void dd_WriteLPResult(FILE *, dd_LPPtr, dd_ErrorType);
void dd_WriteLPErrorMessages(FILE *, dd_ErrorType);
void dd_WriteLPTimes(FILE *, dd_LPPtr);


/* end of cddlib.h */
