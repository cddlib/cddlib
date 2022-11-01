/* cdd.h: Header file for cddlib.c 
   written by Komei Fukuda, fukuda@math.ethz.ch
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

#if defined(__cplusplus)
extern "C" {
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

/* basic matrix manipulations */
void dd_InitializeArow(dd_colrange d,dd_Arow * a);
void dd_InitializeAmatrix(dd_rowrange m,dd_colrange d,dd_Amatrix * A);
void dd_InitializeBmatrix(dd_colrange d, dd_Bmatrix * B);
dd_SetFamilyPtr dd_CreateSetFamily(dd_bigrange fsize,dd_bigrange ssize);
void dd_FreeSetFamily(dd_SetFamilyPtr F);
dd_MatrixPtr dd_CreateMatrix(dd_rowrange m_size,dd_colrange d_size);
void dd_FreeAmatrix(dd_rowrange m,dd_colrange d,dd_Amatrix A);
void dd_FreeArow(dd_colrange d, dd_Arow a);
void dd_FreeBmatrix(dd_colrange d,dd_Bmatrix B);
void dd_FreeDDMemory(dd_PolyhedraPtr poly);
void dd_FreePolyhedra(dd_PolyhedraPtr poly);
void dd_FreeMatrix(dd_MatrixPtr M);
void dd_SetToIdentity(dd_colrange d_size, dd_Bmatrix T);

/* sign recognitions */
dd_boolean dd_Nonnegative(mytype val);
dd_boolean dd_Nonpositive(mytype val);
dd_boolean dd_Positive(mytype val);
dd_boolean dd_Negative(mytype val);
dd_boolean dd_EqualToZero(mytype val);
dd_boolean dd_Nonzero(mytype val);
dd_boolean dd_Equal(mytype val1,mytype val2);
dd_boolean dd_Larger(mytype val1,mytype val2);
dd_boolean dd_Smaller(mytype val1,mytype val2);
void dd_abs(mytype absval, mytype val);
void dd_LinearComb(mytype lc, mytype v1, mytype c1, mytype v2, mytype c2);
void dd_InnerProduct(mytype prod, dd_colrange d, dd_Arow v1, dd_Arow v2);

/* major cddlib operations */
dd_MatrixPtr dd_CopyInput(dd_PolyhedraPtr poly);
dd_MatrixPtr dd_CopyOutput(dd_PolyhedraPtr poly);
dd_MatrixPtr dd_CopyInequalities(dd_PolyhedraPtr poly);
dd_MatrixPtr dd_CopyGenerators(dd_PolyhedraPtr poly);
dd_SetFamilyPtr dd_CopyIncidence(dd_PolyhedraPtr poly);
dd_SetFamilyPtr dd_CopyAdjacency(dd_PolyhedraPtr poly);
dd_SetFamilyPtr dd_CopyInputIncidence(dd_PolyhedraPtr poly);
dd_SetFamilyPtr dd_CopyInputAdjacency(dd_PolyhedraPtr poly);
dd_boolean dd_DDFile2File(char *ifile, char *ofile, dd_ErrorType *err);
dd_boolean dd_DDInputAppend(dd_PolyhedraPtr* poly, dd_MatrixPtr M, dd_ErrorType* err);
dd_MatrixPtr dd_PolyFile2Matrix(FILE *f, dd_ErrorType * Error);

dd_PolyhedraPtr dd_DDMatrix2Poly(dd_MatrixPtr M, dd_ErrorType * err);
dd_PolyhedraPtr dd_DDMatrix2Poly2(dd_MatrixPtr M, dd_RowOrderType horder, dd_ErrorType * err);
dd_boolean dd_Redundant(dd_MatrixPtr M, dd_rowrange itest, dd_Arow certificate, dd_ErrorType * error);  /* 092 */
dd_rowset dd_RedundantRows(dd_MatrixPtr M, dd_ErrorType * error);  /* 092 */
dd_boolean dd_SRedundant(dd_MatrixPtr M, dd_rowrange itest, dd_Arow certificate, dd_ErrorType * error);  /* 093a */
dd_rowset dd_SRedundantRows(dd_MatrixPtr M, dd_ErrorType * error);  /* 093a */
dd_rowset dd_RedundantRowsViaShooting(dd_MatrixPtr M, dd_ErrorType * error); /* 092 */
dd_rowrange dd_RayShooting(dd_MatrixPtr M, dd_Arow p, dd_Arow r);  /* 092 */ 
 /* 092, find the first inequality "hit" by a ray from an intpt.  */
dd_boolean dd_ImplicitLinearity(dd_MatrixPtr M, dd_rowrange itest, dd_Arow certificate, dd_ErrorType * error);  /* 092 */
dd_rowset dd_ImplicitLinearityRows(dd_MatrixPtr M, dd_ErrorType * error);  /* 092  */
int dd_FreeOfImplicitLinearity(dd_MatrixPtr M, dd_Arow certificate, dd_rowset * imp_linrows, dd_ErrorType * error) ; /* 094 */
dd_boolean dd_MatrixCanonicalizeLinearity(dd_MatrixPtr * M, dd_rowset * impl_linset,dd_rowindex * newpos, dd_ErrorType * error); /* 094 */
dd_boolean dd_MatrixCanonicalize(dd_MatrixPtr * M, dd_rowset * impl_linset, dd_rowset * redset, dd_rowindex * newpos, dd_ErrorType * error); /* 094 */
dd_boolean dd_MatrixRedundancyRemove(dd_MatrixPtr *M, dd_rowset *redset,dd_rowindex *newpos, dd_ErrorType * error); /* 094 */
dd_boolean dd_FindRelativeInterior(dd_MatrixPtr M, dd_rowset * ImL, dd_rowset * Lbasis, dd_LPSolutionPtr * lps, dd_ErrorType * err);  /* 094 */
dd_boolean dd_ExistsRestrictedFace(dd_MatrixPtr M, dd_rowset R, dd_rowset S, dd_ErrorType * err);  /* 0.94 */
dd_boolean dd_ExistsRestrictedFace2(dd_MatrixPtr M, dd_rowset R, dd_rowset S, dd_LPSolutionPtr * lps, dd_ErrorType * err); /* 0.94 */

dd_SetFamilyPtr dd_Matrix2Adjacency(dd_MatrixPtr M, dd_ErrorType * error);  /* 093 */
dd_SetFamilyPtr dd_Matrix2WeakAdjacency(dd_MatrixPtr M, dd_ErrorType * error);  /* 093a */
long dd_MatrixRank(dd_MatrixPtr M, dd_rowset ignoredrows, dd_colset ignoredcols, dd_rowset * rowbasis, dd_colset * colbasis);

/* Matrix Basic Operations */
dd_MatrixPtr dd_MatrixCopy(dd_MatrixPtr M); /* a new name for dd_CopyMatrix */
dd_MatrixPtr dd_CopyMatrix(dd_MatrixPtr M); /* 090c, kept for compatibility */
dd_MatrixPtr dd_MatrixNormalizedCopy(dd_MatrixPtr M); /* 094 */
dd_MatrixPtr dd_MatrixNormalizedSortedCopy(dd_MatrixPtr M,dd_rowindex* newpos); /* 094 */
dd_MatrixPtr dd_MatrixUniqueCopy(dd_MatrixPtr M,dd_rowindex* newpos); /* 094 */
dd_MatrixPtr dd_MatrixNormalizedSortedUniqueCopy(dd_MatrixPtr M,dd_rowindex* newpos); /* 094 */
dd_MatrixPtr dd_MatrixSortedUniqueCopy(dd_MatrixPtr M,dd_rowindex* newpos); /* 094 */

dd_MatrixPtr dd_MatrixAppend(dd_MatrixPtr M1, dd_MatrixPtr M2);  /* a name for dd_AppendMatrix */
dd_MatrixPtr dd_AppendMatrix(dd_MatrixPtr M1, dd_MatrixPtr M2);  /* 090c, kept for compatibility */

int dd_MatrixAppendTo(dd_MatrixPtr* M1, dd_MatrixPtr M2);  /* 092 */
int dd_Remove(dd_MatrixPtr*, dd_rowrange);  /* 092 */
dd_MatrixPtr dd_MatrixSubmatrix(dd_MatrixPtr M, dd_rowset delset); /* 092 */
dd_MatrixPtr dd_MatrixSubmatrix2(dd_MatrixPtr M, dd_rowset delset,dd_rowindex* newpos); /* 094.  It returns new row positions. */
dd_MatrixPtr dd_MatrixSubmatrix2L(dd_MatrixPtr M, dd_rowset delset,dd_rowindex* newpos); /* 094.  Linearity shifted up. */
int dd_MatrixShiftupLinearity(dd_MatrixPtr * M,dd_rowindex * newpos); /* 094 */
int dd_MatrixRowRemove(dd_MatrixPtr *M, dd_rowrange r); /* 092 */
int dd_MatrixRowRemove2(dd_MatrixPtr *M, dd_rowrange r,dd_rowindex* newpos); /* 094*/
int dd_MatrixRowsRemove(dd_MatrixPtr *M, dd_rowset delset); /* 094 */
int dd_MatrixRowsRemove2(dd_MatrixPtr *M, dd_rowset delset,dd_rowindex* newpos); /* 094 */

/* input/output */
void dd_SetInputFile(FILE **f,dd_DataFileType inputfile, dd_ErrorType * Error);
void dd_SetWriteFileName(dd_DataFileType inputfile, dd_DataFileType outfile, char cflag, dd_RepresentationType rep);

void dd_WriteAmatrix(FILE * f, dd_Amatrix A, dd_rowrange rowmax, dd_colrange colmax);
void dd_WriteArow(FILE *f, dd_Arow a, dd_colrange d);
void dd_WriteBmatrix(FILE * f, dd_colrange d_size, dd_Bmatrix B);
void dd_WriteMatrix(FILE * f, dd_MatrixPtr M);
void dd_MatrixIntegerFilter(dd_MatrixPtr M);
void dd_WriteReal(FILE * f, mytype x);
void dd_WriteNumber(FILE *f, mytype x); 
    /* write a number depending on the arithmetic used.  */
void dd_WritePolyFile(FILE * f, dd_PolyhedraPtr poly);
void dd_WriteRunningMode(FILE * f, dd_PolyhedraPtr poly);
void dd_WriteErrorMessages(FILE * f, dd_ErrorType Error);
void dd_WriteSetFamily(FILE * f, dd_SetFamilyPtr F);
void dd_WriteSetFamilyCompressed(FILE * f, dd_SetFamilyPtr F);
void dd_WriteProgramDescription(FILE * f);
void dd_WriteDDTimes(FILE * f, dd_PolyhedraPtr poly);
void dd_WriteTimes(FILE * f, time_t starttime, time_t endtime);
void dd_WriteIncidence(FILE * f, dd_PolyhedraPtr poly);
void dd_WriteAdjacency(FILE * f, dd_PolyhedraPtr poly);
void dd_WriteInputAdjacency(FILE * f, dd_PolyhedraPtr poly);
void dd_WriteInputIncidence(FILE * f, dd_PolyhedraPtr poly);

/* functions and types for LP solving */

dd_LPPtr dd_Matrix2LP(dd_MatrixPtr M, dd_ErrorType * err);
  /* Load a matrix to create an LP object. */
  
dd_LPPtr dd_Matrix2Feasibility(dd_MatrixPtr M, dd_ErrorType * err);
  /* Load a matrix to create an LP object for feasibility (obj == 0) .*/  /*  094 */
  
dd_LPPtr dd_Matrix2Feasibility2(dd_MatrixPtr M, dd_rowset R, dd_rowset S, dd_ErrorType * err);
  /* Load a matrix to create an LP object for feasibility with additional equality and
   strict inequality constraints. */  /*  094 */

dd_boolean dd_LPSolve(dd_LPPtr lp,dd_LPSolverType solver,dd_ErrorType * err);
dd_boolean dd_LPSolve0(dd_LPPtr lp,dd_LPSolverType solver,dd_ErrorType * err);
void dd_CrissCrossSolve(dd_LPPtr lp,dd_ErrorType * err);
void dd_DualSimplexSolve(dd_LPPtr lp,dd_ErrorType * err);

dd_LPPtr dd_MakeLPforInteriorFinding(dd_LPPtr lp);  
dd_LPSolutionPtr dd_CopyLPSolution(dd_LPPtr lp);  /* 0.90c */
void dd_WriteLP(FILE * f, dd_LPPtr lp); /* 092 */

dd_LPPtr dd_CreateLPData(dd_LPObjectiveType obj,dd_NumberType nt,dd_rowrange m,dd_colrange d);
int dd_LPReverseRow(dd_LPPtr lp, dd_rowrange i);
    /* reverse the i-th row (1 <= i <= no. of rows) */
int dd_LPReplaceRow(dd_LPPtr lp, dd_rowrange i, dd_Arow a);
    /* replace the i-th row (1 <= i <= no. of rows) */
dd_Arow dd_LPCopyRow(dd_LPPtr lp, dd_rowrange i);
    /* copy the i-th row (1 <= i <= no. of rows) */

void dd_FreeLPData(dd_LPPtr lp);
void dd_FreeLPSolution(dd_LPSolutionPtr lps);

void dd_WriteLPResult(FILE * f, dd_LPPtr lp, dd_ErrorType err);
void dd_WriteLPErrorMessages(FILE *, dd_ErrorType);
void dd_WriteLPTimes(FILE * f, dd_LPPtr lp);
void dd_WriteLPStats(FILE *f);
void dd_WriteLPMode(FILE *f);

dd_MatrixPtr dd_FourierElimination(dd_MatrixPtr M,dd_ErrorType * error);
dd_MatrixPtr dd_BlockElimination(dd_MatrixPtr M, dd_colset delset, dd_ErrorType * error);

/* ---------- FUNCTIONS MEANT TO BE NON-PUBLIC ---------- */
void dd_QuickSort(dd_rowindex OV, long p, long r, dd_Amatrix A, long dmax);
void dd_RandomPermutation(dd_rowindex OV, long t, unsigned int seed);
void dd_UniqueRows(dd_rowindex OV, long p, long r, dd_Amatrix A, long dmax, dd_rowset preferred, long * uniqrows);

dd_boolean dd_DoubleDescription(dd_PolyhedraPtr poly, dd_ErrorType* err);
dd_boolean dd_DoubleDescription2(dd_PolyhedraPtr poly, dd_RowOrderType horder, dd_ErrorType * err);

void dd_FreeDDMemory0(dd_ConePtr cone);
void dd_fread_rational_value (FILE *f, mytype value);
void dd_sread_rational_value (const char *s, mytype value);
void dd_AddNewHalfspace1(dd_ConePtr cone, dd_rowrange hnew);
void dd_AddNewHalfspace2(dd_ConePtr cone, dd_rowrange hnew);
void dd_AddRay(dd_ConePtr cone, mytype * p);
void dd_AddArtificialRay(dd_ConePtr cone);
void dd_AValue(mytype* val,dd_colrange d_size, dd_Amatrix A, mytype * p, dd_rowrange i);
void dd_CheckAdjacency(dd_ConePtr cone, dd_RayPtr* RP1, dd_RayPtr* RP2, dd_boolean * adjacent);
void dd_CheckEquality(dd_colrange d_size, dd_RayPtr * RP1, dd_RayPtr * RP2, dd_boolean * equal);
void dd_ComputeRowOrderVector(dd_ConePtr cone);
void dd_ConditionalAddEdge(dd_ConePtr cone,dd_RayPtr Ray1, dd_RayPtr Ray2, dd_RayPtr ValidFirstRay);
void dd_CopyArow(mytype * acopy, mytype * a, dd_colrange d);
void dd_CopyNormalizedAmatrix(mytype ** Acopy, mytype ** A, dd_rowrange m, dd_colrange d);
void dd_CopyNormalizedArow(mytype * acopy, mytype * a, dd_colrange d);
void dd_CopyAmatrix(mytype ** Acopy, mytype ** A, dd_rowrange m, dd_colrange d);
void dd_PermuteCopyAmatrix(mytype ** Acopy, mytype ** A, dd_rowrange m, dd_colrange d, dd_rowindex roworder);
void dd_PermutePartialCopyAmatrix(mytype ** Acopy, mytype ** A, dd_rowrange m, dd_colrange d, dd_rowindex roworder, dd_rowrange p, dd_rowrange q);
void dd_SetMatrixObjective(dd_MatrixPtr M, dd_LPObjectiveType objective);
void dd_SetMatrixNumberType(dd_MatrixPtr M, dd_NumberType numbtype);
void dd_SetMatrixRepresentationType(dd_MatrixPtr M, dd_RepresentationType representation);
void dd_CopyBmatrix(dd_colrange d_size, dd_Bmatrix T, dd_Bmatrix TCOPY);
void dd_CopyRay(mytype * a, dd_colrange d_origsize, dd_RayPtr RR, dd_RepresentationType rep, dd_colindex reducedcol);
void dd_CreateInitialEdges(dd_ConePtr cone);
void dd_CreateNewRay(dd_ConePtr cone, dd_RayPtr Ptr1, dd_RayPtr Ptr2, dd_rowrange ii);
void dd_Eliminate(dd_ConePtr cone, dd_RayPtr* Ptr);
void dd_EvaluateARay1(dd_rowrange i, dd_ConePtr cone);
void dd_EvaluateARay2(dd_rowrange i, dd_ConePtr cone);
void dd_FeasibilityIndices(long * fnum, long * infnum, dd_rowrange i, dd_ConePtr cone);
void dd_FindBasis(dd_ConePtr cone, long *rank);
void dd_FindInitialRays(dd_ConePtr cone, dd_boolean * found);
void dd_ColumnReduce(dd_ConePtr cone);
void dd_GaussianColumnPivot(dd_rowrange m_size, dd_colrange d_size, dd_Amatrix X, dd_Bmatrix T,  dd_rowrange r, dd_colrange s);
dd_boolean dd_LexSmaller(mytype * v1, mytype * v2, long dmax);
dd_boolean dd_LexLarger(mytype * v1, mytype * v2, long dmax);
dd_boolean dd_LexEqual(mytype * v1, mytype * v2, long dmax);
void dd_Normalize(dd_colrange d_size, mytype * V);
void dd_MatrixIntegerFilter(dd_MatrixPtr M);
void dd_ProcessCommandLine(FILE* f,dd_MatrixPtr M, const char * line);
void dd_SelectNextHalfspace(dd_ConePtr cone, dd_rowset excluded, dd_rowrange * hh);
void dd_SelectPivot2(dd_rowrange m_size,dd_colrange d_size,dd_Amatrix A,
dd_Bmatrix T,dd_RowOrderType roworder,dd_rowindex ordervec, dd_rowset equalityset,dd_rowrange rowmax,dd_rowset NopivotRow,
dd_colset NopivotCol,dd_rowrange * r,dd_colrange * s,dd_boolean * selected);
void dd_SelectPreorderedNext(dd_ConePtr cone, dd_rowset excluded, dd_rowrange * hh);
void dd_SetInequalitySets(dd_ConePtr cone);
void dd_SnapToInteger(mytype y, mytype x);
void dd_StoreRay1(dd_ConePtr cone, mytype * p, dd_boolean * feasible);
void dd_StoreRay2(dd_ConePtr cone, mytype * p, dd_boolean * feasible, dd_boolean * weaklyfeasible);
void dd_TableauEntry(mytype * x, dd_rowrange m_size, dd_colrange d_size, dd_Amatrix X, dd_Bmatrix T, dd_rowrange r, dd_colrange s);
void dd_UpdateEdges(dd_ConePtr cone, dd_RayPtr RRbegin, dd_RayPtr RRend);
void dd_UpdateRowOrderVector(dd_ConePtr cone, dd_rowset PriorityRows);
void dd_WriteRay(FILE * f, dd_colrange d_origsize, dd_RayPtr RR,
   dd_RepresentationType rep, dd_colindex reducedcol);
void dd_ZeroIndexSet(dd_rowrange m_size, dd_colrange d_size, dd_Amatrix A, mytype * x, dd_rowset ZS);

/* New functions to handle data loading, NON-PUBLIC */
dd_NumberType dd_GetNumberType(const char * line);
dd_ConePtr dd_ConeDataLoad(dd_PolyhedraPtr poly);
dd_PolyhedraPtr dd_CreatePolyhedraData(dd_rowrange m, dd_colrange d);
dd_boolean dd_InitializeConeData(dd_rowrange m, dd_colrange d, dd_ConePtr* cone);
dd_boolean dd_AppendMatrix2Poly(dd_PolyhedraPtr* poly, dd_MatrixPtr M);

#if defined(__cplusplus)
}
#endif

/* end of cddlib.h */
