/* cddlib.h: Header file for cddlib.c 
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.8alpha, June 1998
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

#ifndef  __CDDTYPES_H
#include "cddtypes.h"
#endif  /* __CDDTYPES_H */

/* ---------- FUNCTIONS MEANT TO BE PUBLIC ---------- */

/* basic matrix manipulations */
void dd_InitializeArow(dd_colrange,dd_Arow *);
void dd_InitializeAmatrix(dd_rowrange,dd_colrange,dd_Amatrix *);
void dd_InitializeBmatrix(dd_colrange, dd_Bmatrix *);
dd_SetFamilyPtr dd_CreateSetFamily(dd_bigrange,dd_bigrange);
void dd_FreeSetFamily(dd_SetFamilyPtr *);
dd_MatrixPtr dd_CreateMatrix(dd_rowrange,dd_colrange);
void dd_FreeAmatrix(dd_rowrange,dd_Amatrix *);
void dd_FreeBmatrix(dd_colrange,dd_Bmatrix *);
void dd_FreeDDMemory(dd_PolyhedraPtr);
void dd_FreeMatrix(dd_MatrixPtr *);
void dd_SetToIdentity(dd_colrange, dd_Bmatrix);

/* sign recognitions */
boolean dd_Nonnegative(double);
boolean dd_Nonpositive(double);
boolean dd_Positive(double);
boolean dd_Negative(double);
boolean dd_Zero(double);
boolean dd_Nonzero(double);

/* major cddlib operations */
dd_MatrixPtr dd_CopyInequalities(dd_PolyhedraPtr);
dd_MatrixPtr dd_CopyGenerators(dd_PolyhedraPtr);
dd_SetFamilyPtr dd_CopyIncidence(dd_PolyhedraPtr);
dd_SetFamilyPtr dd_CopyAdjacency(dd_PolyhedraPtr);
boolean dd_DoubleDescription(dd_PolyhedraPtr);
boolean dd_DDAddInequalities(dd_PolyhedraPtr, dd_MatrixPtr);
boolean dd_PolyhedraInput(dd_ErrorType*, dd_PolyhedraPtr *);
void dd_PolyhedraLoadMatrix(dd_PolyhedraPtr *, 
  dd_RepresentationType, dd_MatrixPtr );

/* output */
void dd_WriteAmatrix(FILE *, dd_Amatrix, dd_rowrange, dd_colrange);
void dd_WriteBmatrix(FILE *, dd_colrange, dd_Bmatrix T);
void dd_WriteMatrix(FILE *, dd_MatrixPtr);
void dd_MatrixIntegerFilter(dd_MatrixPtr);
void dd_WriteReal(FILE *, double);
void dd_WritePolyhedraFile(FILE *f, dd_PolyhedraPtr poly);
void dd_WriteRunningMode(FILE *, dd_PolyhedraPtr poly);
void dd_WriteErrorMessages(FILE *, dd_ErrorType);
void dd_WriteSetFamily(FILE *, dd_SetFamilyPtr);
void dd_WriteSetFamilyWithNumbers(FILE *, dd_SetFamilyPtr);


/* ---------- FUNCTIONS MEANT TO BE NON-PUBLIC ---------- */

void AddNewHalfspace1(dd_ConePtr, dd_rowrange);
void AddNewHalfspace2(dd_ConePtr, dd_rowrange);
void AddRay(dd_ConePtr, double *);
void AddArtificialRay(dd_ConePtr);
double AValue(dd_colrange, dd_Amatrix, double *, dd_rowrange);
void CheckAdjacency(dd_ConePtr,
    dd_RayPtr*, dd_RayPtr*, boolean *);
void CheckEquality(dd_colrange, dd_RayPtr *, dd_RayPtr *, boolean *);
void ComputeRowOrderVector(dd_ConePtr);
void ConditionalAddEdge(dd_ConePtr,dd_RayPtr, dd_RayPtr, dd_RayPtr);
void CopyArow(double *, double *, dd_colrange);
void CopyAmatrix(double **, double **, dd_rowrange, dd_colrange);
void CopyBmatrix(dd_colrange, dd_Bmatrix T, dd_Bmatrix TCOPY);
void CopyRay(double *, dd_colrange, dd_RayPtr,
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
boolean LexSmaller(double *, double *, long);
boolean LexLarger(double *, double *, long);
void Normalize(dd_colrange, double *);
void MatrixIntegerFilter(dd_MatrixPtr);
void ProcessCommandLine(dd_PolyhedraPtr , char *);
void SelectNextHalfspace(dd_ConePtr, dd_rowset, dd_rowrange *);
void SelectPivot(dd_ConePtr, dd_rowrange, dd_rowset,
   dd_colset, dd_rowrange *, dd_colrange *, boolean *selected);
void SelectPreorderedNext(dd_ConePtr, dd_rowset, dd_rowrange *);
void SetInequalitySets(dd_ConePtr);
double SnapToInteger(double);
void StoreRay1(dd_ConePtr, double *, boolean *);
void StoreRay2(dd_ConePtr, double *, boolean *, boolean *);
double TableauEntry(dd_rowrange, dd_colrange, dd_Amatrix, dd_Bmatrix T, dd_rowrange, dd_colrange);
void UpdateEdges(dd_ConePtr, dd_RayPtr, dd_RayPtr);
void UpdateRowOrderVector(dd_ConePtr, dd_rowset PriorityRows);
void WriteIncidence(FILE *, dd_ConePtr, dd_RayPtr);
void WriteRay(FILE *, dd_colrange, dd_RayPtr,
   dd_RepresentationType, dd_colindex);
void ZeroIndexSet(dd_rowrange, dd_colrange, dd_Amatrix, double *, dd_rowset);

/* New functions to handle data loading, NON-PUBLIC */
dd_NumberType GetNumberType(char *);
dd_ConePtr ConeDataLoad(dd_PolyhedraPtr);
dd_PolyhedraPtr CreatePolyhedraData(dd_rowrange, dd_colrange);
boolean InitializeConeData(dd_rowrange, dd_colrange, dd_ConePtr*);
void AddInequalities(dd_PolyhedraPtr, dd_Matrix*);


/* end of cddlib.h */
