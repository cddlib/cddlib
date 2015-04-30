/* cdd.h: Header file for cddlib.c 
   written by Komei Fukuda, fukuda@math.ethz.ch
   Version 0.94h, April 30, 2015
*/

/* cddlib.c : C-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddlibman.tex for detail.
*/

#ifndef  __CDD_H
#define  __CDD_H
#endif  /* __CDD_H */

#ifndef  __CDDMP_H
#include "cddmp.h"
#endif  /* __CDDMP_H */

#ifndef  __CDDTYPES_H
#include "cddtypes.h"
#endif  /* __CDDTYPES_H */

#ifdef GMPRATIONAL
#ifndef __CDD_HF
#include "cdd_f.h"
#endif
#endif

/* GLOBAL CONSTANTS and STATISTICS VARIABLES (to be set by dd_set_global_constants() */
extern mytype dd_zero;
extern mytype dd_one;
extern mytype dd_purezero;
extern mytype dd_minuszero;
extern mytype dd_minusone;

extern time_t dd_statStartTime; /* cddlib starting time */
extern long dd_statBApivots;  /* basis finding pivots */
extern long dd_statCCpivots;  /* criss-cross pivots */
extern long dd_statDS1pivots; /* phase 1 pivots */
extern long dd_statDS2pivots; /* phase 2 pivots */
extern long dd_statACpivots;  /* anticycling (cc) pivots */
#ifdef GMPRATIONAL
extern long dd_statBSpivots;  /* basis status checking pivots */
#endif
extern dd_LPSolverType dd_choiceLPSolverDefault;  /* Default LP solver Algorithm */
extern dd_LPSolverType dd_choiceRedcheckAlgorithm;  /* Redundancy Checking Algorithm */
extern dd_boolean dd_choiceLexicoPivotQ;    /* whether to use the lexicographic pivot */

   /* to be used to avoid creating temporary spaces for mytype */
#define dd_almostzero  1.0E-7

/* ---------- FUNCTIONS MEANT TO BE PUBLIC ---------- */

#if defined(__cplusplus)
extern "C" {
#endif

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
dd_boolean dd_Nonnegative(mytype);
dd_boolean dd_Nonpositive(mytype);
dd_boolean dd_Positive(mytype);
dd_boolean dd_Negative(mytype);
dd_boolean dd_EqualToZero(mytype);
dd_boolean dd_Nonzero(mytype);
dd_boolean dd_Equal(mytype,mytype);
dd_boolean dd_Larger(mytype,mytype);
dd_boolean dd_Smaller(mytype,mytype);
void dd_abs(mytype, mytype);
void dd_LinearComb(mytype, mytype, mytype, mytype, mytype);
void dd_InnerProduct(mytype, dd_colrange, dd_Arow, dd_Arow);

/* major cddlib operations */
dd_MatrixPtr dd_CopyInput(dd_PolyhedraPtr);
dd_MatrixPtr dd_CopyOutput(dd_PolyhedraPtr);
dd_MatrixPtr dd_CopyInequalities(dd_PolyhedraPtr);
dd_MatrixPtr dd_CopyGenerators(dd_PolyhedraPtr);
dd_SetFamilyPtr dd_CopyIncidence(dd_PolyhedraPtr);
dd_SetFamilyPtr dd_CopyAdjacency(dd_PolyhedraPtr);
dd_SetFamilyPtr dd_CopyInputIncidence(dd_PolyhedraPtr);
dd_SetFamilyPtr dd_CopyInputAdjacency(dd_PolyhedraPtr);
dd_boolean dd_DDFile2File(char *ifile, char *ofile, dd_ErrorType *err);
dd_boolean dd_DDInputAppend(dd_PolyhedraPtr*, dd_MatrixPtr, dd_ErrorType*);
dd_MatrixPtr dd_PolyFile2Matrix(FILE *f, dd_ErrorType *);

dd_PolyhedraPtr dd_DDMatrix2Poly(dd_MatrixPtr, dd_ErrorType *);
dd_PolyhedraPtr dd_DDMatrix2Poly2(dd_MatrixPtr, dd_RowOrderType, dd_ErrorType *);
dd_boolean dd_Redundant(dd_MatrixPtr, dd_rowrange, dd_Arow, dd_ErrorType *);  /* 092 */
dd_rowset dd_RedundantRows(dd_MatrixPtr, dd_ErrorType *);  /* 092 */
dd_boolean dd_SRedundant(dd_MatrixPtr, dd_rowrange, dd_Arow, dd_ErrorType *);  /* 093a */
dd_rowset dd_SRedundantRows(dd_MatrixPtr, dd_ErrorType *);  /* 093a */
dd_rowset dd_RedundantRowsViaShooting(dd_MatrixPtr, dd_ErrorType *); /* 092 */
dd_rowrange dd_RayShooting(dd_MatrixPtr, dd_Arow intpt, dd_Arow direction);  /* 092 */ 
 /* 092, find the first inequality "hit" by a ray from an intpt.  */
dd_boolean dd_ImplicitLinearity(dd_MatrixPtr, dd_rowrange, dd_Arow, dd_ErrorType *);  /* 092 */
dd_rowset dd_ImplicitLinearityRows(dd_MatrixPtr, dd_ErrorType *);  /* 092  */
int dd_FreeOfImplicitLinearity(dd_MatrixPtr, dd_Arow, dd_rowset *, dd_ErrorType *) ; /* 094 */
dd_boolean dd_MatrixCanonicalizeLinearity(dd_MatrixPtr *, dd_rowset *,dd_rowindex *, dd_ErrorType *); /* 094 */
dd_boolean dd_MatrixCanonicalize(dd_MatrixPtr *, dd_rowset *, dd_rowset *, dd_rowindex *, dd_ErrorType *); /* 094 */
dd_boolean dd_MatrixRedundancyRemove(dd_MatrixPtr *M, dd_rowset *redset,dd_rowindex *newpos, dd_ErrorType *); /* 094 */
dd_boolean dd_FindRelativeInterior(dd_MatrixPtr, dd_rowset *, dd_rowset *, dd_LPSolutionPtr *, dd_ErrorType *);  /* 094 */
dd_boolean dd_ExistsRestrictedFace(dd_MatrixPtr, dd_rowset, dd_rowset, dd_ErrorType *);  /* 0.94 */
dd_boolean dd_ExistsRestrictedFace2(dd_MatrixPtr, dd_rowset, dd_rowset, dd_LPSolutionPtr *, dd_ErrorType *); /* 0.94 */

dd_SetFamilyPtr dd_Matrix2Adjacency(dd_MatrixPtr, dd_ErrorType *);  /* 093 */
dd_SetFamilyPtr dd_Matrix2WeakAdjacency(dd_MatrixPtr, dd_ErrorType *);  /* 093a */
long dd_MatrixRank(dd_MatrixPtr, dd_rowset, dd_colset, dd_rowset *, dd_colset *);

/* Matrix Basic Operations */
dd_MatrixPtr dd_MatrixCopy(dd_MatrixPtr); /* a new name for dd_CopyMatrix */
dd_MatrixPtr dd_CopyMatrix(dd_MatrixPtr); /* 090c, kept for compatibility */
dd_MatrixPtr dd_MatrixNormalizedCopy(dd_MatrixPtr); /* 094 */
dd_MatrixPtr dd_MatrixNormalizedSortedCopy(dd_MatrixPtr,dd_rowindex*); /* 094 */
dd_MatrixPtr dd_MatrixUniqueCopy(dd_MatrixPtr,dd_rowindex*); /* 094 */
dd_MatrixPtr dd_MatrixNormalizedSortedUniqueCopy(dd_MatrixPtr,dd_rowindex*); /* 094 */
dd_MatrixPtr dd_MatrixSortedUniqueCopy(dd_MatrixPtr,dd_rowindex*); /* 094 */

dd_MatrixPtr dd_MatrixAppend(dd_MatrixPtr, dd_MatrixPtr);  /* a name for dd_AppendMatrix */
dd_MatrixPtr dd_AppendMatrix(dd_MatrixPtr, dd_MatrixPtr);  /* 090c, kept for compatibility */

int dd_MatrixAppendTo(dd_MatrixPtr*, dd_MatrixPtr);  /* 092 */
int dd_Remove(dd_MatrixPtr*, dd_rowrange);  /* 092 */
dd_MatrixPtr dd_MatrixSubmatrix(dd_MatrixPtr, dd_rowset delset); /* 092 */
dd_MatrixPtr dd_MatrixSubmatrix2(dd_MatrixPtr, dd_rowset delset,dd_rowindex*); /* 094.  It returns new row positions. */
dd_MatrixPtr dd_MatrixSubmatrix2L(dd_MatrixPtr, dd_rowset delset,dd_rowindex*); /* 094.  Linearity shifted up. */
int dd_MatrixShiftupLinearity(dd_MatrixPtr *,dd_rowindex *); /* 094 */
int dd_MatrixRowRemove(dd_MatrixPtr *M, dd_rowrange r); /* 092 */
int dd_MatrixRowRemove2(dd_MatrixPtr *M, dd_rowrange r,dd_rowindex*); /* 094*/
int dd_MatrixRowsRemove(dd_MatrixPtr *M, dd_rowset delset); /* 094 */
int dd_MatrixRowsRemove2(dd_MatrixPtr *M, dd_rowset delset,dd_rowindex*); /* 094 */

/* input/output */
void dd_SetInputFile(FILE **f,dd_DataFileType inputfile, dd_ErrorType *);
void dd_SetWriteFileName(dd_DataFileType, dd_DataFileType, char, dd_RepresentationType);

void dd_WriteAmatrix(FILE *, dd_Amatrix, dd_rowrange, dd_colrange);
void dd_WriteArow(FILE *f, dd_Arow a, dd_colrange);
void dd_WriteBmatrix(FILE *, dd_colrange, dd_Bmatrix T);
void dd_WriteMatrix(FILE *, dd_MatrixPtr);
void dd_MatrixIntegerFilter(dd_MatrixPtr);
void dd_WriteReal(FILE *, mytype);
void dd_WriteNumber(FILE *f, mytype x); 
    /* write a number depending on the arithmetic used.  */
void dd_WritePolyFile(FILE *, dd_PolyhedraPtr);
void dd_WriteRunningMode(FILE *, dd_PolyhedraPtr);
void dd_WriteErrorMessages(FILE *, dd_ErrorType);
void dd_WriteSetFamily(FILE *, dd_SetFamilyPtr);
void dd_WriteSetFamilyCompressed(FILE *, dd_SetFamilyPtr);
void dd_WriteProgramDescription(FILE *);
void dd_WriteDDTimes(FILE *, dd_PolyhedraPtr);
void dd_WriteTimes(FILE *, time_t, time_t);
void dd_WriteIncidence(FILE *, dd_PolyhedraPtr);
void dd_WriteAdjacency(FILE *, dd_PolyhedraPtr);
void dd_WriteInputAdjacency(FILE *, dd_PolyhedraPtr);
void dd_WriteInputIncidence(FILE *, dd_PolyhedraPtr);

/* functions and types for LP solving */

dd_LPPtr dd_Matrix2LP(dd_MatrixPtr, dd_ErrorType *);
  /* Load a matrix to create an LP object. */
  
dd_LPPtr dd_Matrix2Feasibility(dd_MatrixPtr, dd_ErrorType *);
  /* Load a matrix to create an LP object for feasibility (obj == 0) .*/  /*  094 */
  
dd_LPPtr dd_Matrix2Feasibility2(dd_MatrixPtr, dd_rowset, dd_rowset, dd_ErrorType *);
  /* Load a matrix to create an LP object for feasibility with additional equality and
   strict inequality constraints. */  /*  094 */

dd_boolean dd_LPSolve(dd_LPPtr,dd_LPSolverType,dd_ErrorType *);
dd_boolean dd_LPSolve0(dd_LPPtr,dd_LPSolverType,dd_ErrorType *);
void dd_CrissCrossSolve(dd_LPPtr lp,dd_ErrorType *);
void dd_DualSimplexSolve(dd_LPPtr lp,dd_ErrorType *);

dd_LPPtr dd_MakeLPforInteriorFinding(dd_LPPtr);  
dd_LPSolutionPtr dd_CopyLPSolution(dd_LPPtr);  /* 0.90c */
void dd_WriteLP(FILE *, dd_LPPtr); /* 092 */

dd_LPPtr dd_CreateLPData(dd_LPObjectiveType,dd_NumberType,dd_rowrange,dd_colrange);
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
void dd_WriteLPStats(FILE *f);
void dd_WriteLPMode(FILE *f);

dd_MatrixPtr dd_FourierElimination(dd_MatrixPtr,dd_ErrorType *);
dd_MatrixPtr dd_BlockElimination(dd_MatrixPtr, dd_colset, dd_ErrorType *);

#if defined(__cplusplus)
}
#endif

/* ---------- FUNCTIONS MEANT TO BE NON-PUBLIC ---------- */
void dd_QuickSort(dd_rowindex, long, long, dd_Amatrix, long);
void dd_RandomPermutation(dd_rowindex, long, unsigned int seed);
void dd_UniqueRows(dd_rowindex, long, long, dd_Amatrix, long, dd_rowset, long *);

dd_boolean dd_DoubleDescription(dd_PolyhedraPtr, dd_ErrorType*);
dd_boolean dd_DoubleDescription2(dd_PolyhedraPtr, dd_RowOrderType, dd_ErrorType *);

void dd_FreeDDMemory0(dd_ConePtr);
void dd_fread_rational_value (FILE *f, mytype value);
void dd_sread_rational_value (const char *s, mytype value);
void dd_AddNewHalfspace1(dd_ConePtr, dd_rowrange);
void dd_AddNewHalfspace2(dd_ConePtr, dd_rowrange);
void dd_AddRay(dd_ConePtr, mytype *);
void dd_AddArtificialRay(dd_ConePtr);
void dd_AValue(mytype*,dd_colrange, dd_Amatrix, mytype *, dd_rowrange);
void dd_CheckAdjacency(dd_ConePtr, dd_RayPtr*, dd_RayPtr*, dd_boolean *);
void dd_CheckEquality(dd_colrange, dd_RayPtr *, dd_RayPtr *, dd_boolean *);
void dd_ComputeRowOrderVector(dd_ConePtr);
void dd_ConditionalAddEdge(dd_ConePtr,dd_RayPtr, dd_RayPtr, dd_RayPtr);
void dd_CopyArow(mytype *, mytype *, dd_colrange);
void dd_CopyNormalizedAmatrix(mytype **, mytype **, dd_rowrange, dd_colrange);
void dd_CopyNormalizedArow(mytype *, mytype *, dd_colrange);
void dd_CopyAmatrix(mytype **, mytype **, dd_rowrange, dd_colrange);
void dd_PermuteCopyAmatrix(mytype **, mytype **, dd_rowrange, dd_colrange, dd_rowindex);
void dd_PermutePartialCopyAmatrix(mytype **, mytype **, dd_rowrange, dd_colrange, dd_rowindex,dd_rowrange, dd_rowrange);
void dd_CopyBmatrix(dd_colrange, dd_Bmatrix T, dd_Bmatrix TCOPY);
void dd_CopyRay(mytype *, dd_colrange, dd_RayPtr,
   dd_RepresentationType, dd_colindex);
void dd_CreateInitialEdges(dd_ConePtr);
void dd_CreateNewRay(dd_ConePtr, dd_RayPtr, dd_RayPtr, dd_rowrange);
void dd_Eliminate(dd_ConePtr, dd_RayPtr*);
void dd_EvaluateARay1(dd_rowrange, dd_ConePtr);
void dd_EvaluateARay2(dd_rowrange, dd_ConePtr);
void dd_FeasibilityIndices(long *, long *, dd_rowrange, dd_ConePtr);
void dd_FindBasis(dd_ConePtr, long *rank);
void dd_FindInitialRays(dd_ConePtr, dd_boolean *);
void dd_ColumnReduce(dd_ConePtr);
void dd_GaussianColumnPivot(dd_rowrange, dd_colrange, dd_Amatrix, dd_Bmatrix,  dd_rowrange, dd_colrange);
dd_boolean dd_LexSmaller(mytype *, mytype *, long);
dd_boolean dd_LexLarger(mytype *, mytype *, long);
dd_boolean dd_LexEqual(mytype *, mytype *, long);
void dd_Normalize(dd_colrange, mytype *);
void dd_MatrixIntegerFilter(dd_MatrixPtr);
void dd_ProcessCommandLine(FILE*,dd_MatrixPtr, const char *);
void dd_SelectNextHalfspace(dd_ConePtr, dd_rowset, dd_rowrange *);
void dd_SelectPivot2(dd_rowrange,dd_colrange,dd_Amatrix,
dd_Bmatrix,dd_RowOrderType,dd_rowindex, dd_rowset,dd_rowrange,dd_rowset,
dd_colset,dd_rowrange *,dd_colrange *,dd_boolean *);
void dd_SelectPreorderedNext(dd_ConePtr, dd_rowset, dd_rowrange *);
void dd_SetInequalitySets(dd_ConePtr);
void dd_SnapToInteger(mytype, mytype);
void dd_StoreRay1(dd_ConePtr, mytype *, dd_boolean *);
void dd_StoreRay2(dd_ConePtr, mytype *, dd_boolean *, dd_boolean *);
void dd_TableauEntry(mytype *, dd_rowrange, dd_colrange, dd_Amatrix, dd_Bmatrix T, dd_rowrange, dd_colrange);
void dd_UpdateEdges(dd_ConePtr, dd_RayPtr, dd_RayPtr);
void dd_UpdateRowOrderVector(dd_ConePtr, dd_rowset PriorityRows);
void dd_WriteRay(FILE *, dd_colrange, dd_RayPtr,
   dd_RepresentationType, dd_colindex);
void dd_ZeroIndexSet(dd_rowrange, dd_colrange, dd_Amatrix, mytype *, dd_rowset);

/* New functions to handle data loading, NON-PUBLIC */
dd_NumberType dd_GetNumberType(const char *);
dd_ConePtr dd_ConeDataLoad(dd_PolyhedraPtr);
dd_PolyhedraPtr dd_CreatePolyhedraData(dd_rowrange, dd_colrange);
dd_boolean dd_InitializeConeData(dd_rowrange, dd_colrange, dd_ConePtr*);
dd_boolean dd_AppendMatrix2Poly(dd_PolyhedraPtr*, dd_MatrixPtr);





/* end of cddlib.h */
