/* cddtypes.h: Header file for cddlib.c 
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.8alpha, June 1998
*/

/* cddlib.c : C-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#ifndef  __CDDTYPES_H
#define  __CDDTYPES_H
#endif  /* __CDDTYPES_H */


#define COPYRIGHT   "Copyright (C) 1996, Komei Fukuda, fukuda@ifor.math.ethz.ch"
#define DDVERSION   "Version 0.8alpha (June 1998)"
#include <time.h>
#include "dplex.h"

#define dd_datawidth       10
#define dd_namelenmax      256 
#define dd_wordlenmax      128 
#define dd_linelenmax      256

#define FALSE 0
#define TRUE 1

#define dd_zero            dp_zero     /*real zero*/
#define dd_magfac          dp_magfac   /*reciprocal of dp_zero*/

typedef char boolean;

typedef long dd_rowrange;
typedef long dd_colrange;
typedef long dd_bigrange;

typedef set_type dd_rowset;
typedef set_type dd_colset;
typedef long *dd_rowindex;   
typedef int *dd_rowflag;   
typedef long *dd_colindex;
typedef double **dd_Amatrix;
typedef double *dd_Arow;
typedef set_type *dd_SetVector;
typedef double **dd_Bmatrix;

typedef char dd_DataFileType[dd_namelenmax];
typedef char dd_LineType[dd_linelenmax];
typedef char dd_WordType[dd_wordlenmax];

typedef struct dd_RayData *dd_RayPtr;

typedef struct dd_RayData {
  double *Ray;
  dd_rowset ZeroSet;
  dd_rowrange FirstInfeasIndex;  /* the first inequality the ray violates */
  boolean feasible;  /* flag to store the feasibility */
  double ARay;   /* temporary area to store some row of A*Ray */
  dd_RayPtr Next;
} dd_Ray;

typedef struct dd_AdjacencyData *dd_AdjacencyPtr;

typedef struct dd_AdjacencyData {
  dd_RayPtr Ray1, Ray2;
  dd_AdjacencyPtr Next;
} dd_Adjacency;

typedef enum {
  Combinatorial, Algebraic
} dd_AdjacencyTestType;

typedef enum {
  MaxIndex, MinIndex, MinCutoff, MaxCutoff, MixCutoff,
   LexMin, LexMax, RandomRow, LineShelling
} dd_HalfspaceOrderType;

typedef enum {
  Real, Rational, Integer, Unknown
} dd_NumberType;

typedef enum {
  Inequality, Generator
} dd_RepresentationType;

typedef enum {
  IneToGen, GenToIne, LPmax, LPmin, InteriorFind
} dd_ConversionType;

typedef enum {
  CrissCross,DualSimplex,CombMaxImprove
} dd_LPsolverType;

typedef enum {
  IncOff=0, IncCardinality, IncSet
} dd_IncidenceOutputType;

typedef enum {
  AdjOff=0, AdjacencyList,  AdjacencyDegree
} dd_AdjacencyOutputType;

typedef enum {
  Auto, SemiAuto, Manual
} dd_FileInputModeType;   
   /* Auto if a input filename is specified by command arguments */

typedef enum {
  DimensionTooLarge, LowColumnRank, ImproperInputFormat, EmptyVrepresentation,
  FileNotFound, None
} dd_ErrorType;

typedef enum {
  InProgress, AllFound, RegionEmpty
} dd_CompStatusType;

typedef struct dd_MatrixData *dd_MatrixPtr;

typedef struct dd_MatrixData {
  dd_rowrange rowsize;
  dd_rowset linset; 
    /*  a subset of rows of linearity (ie, generators of
        linearity space for V-representation, and equations
        for H-representation. */
  dd_colrange colsize;
  dd_NumberType number;
  dd_Amatrix matrix;  
} dd_Matrix;

typedef struct dd_SetFamilyData *dd_SetFamilyPtr;

typedef struct dd_SetFamilyData {
  dd_bigrange famsize;
  dd_bigrange setsize;
  dd_SetVector set;  
} dd_SetFamily;

typedef struct dd_NodeData *dd_NodePtr;
typedef struct dd_NodeData {dd_bigrange key; dd_NodePtr next;} dd_Node;

typedef struct dd_GraphData *dd_GraphPtr;

typedef struct dd_GraphData {
  dd_bigrange vsize;
  dd_NodePtr *adjlist;  /* should be initialized to have vsize components */
} dd_Graph;

typedef struct dd_PolyhedraData *dd_PolyhedraPtr;
typedef struct dd_ConeData *dd_ConePtr;

typedef struct dd_PolyhedraData {
  dd_RepresentationType Representation;  /* given representation */
  boolean Homogeneous;
  dd_colrange d;
  dd_rowrange m;
  dd_Amatrix A;   /* Inequality System:  m times d matrix */
  dd_NumberType Number;
  dd_ConePtr child;  /* pointing to the homogenized cone data */
  dd_rowrange m_alloc; /* allocated row size of matrix A */
  dd_colrange d_alloc; /* allocated col size of matrix A */
  dd_Arow c;           /* cost vector */

  dd_rowflag EqualityIndex;  
    /* ith component is 1 if it is equality, -1 if it is strict inequality, 0 otherwise. */

  boolean NondegAssumed;
  boolean InitBasisAtBottom;
  boolean RestrictedEnumeration;
  boolean RelaxedEnumeration;

} dd_Polyhedra;


typedef struct dd_ConeData {
  dd_RepresentationType Representation;
  dd_rowrange m;
  dd_colrange d;
  dd_Amatrix A;
  dd_NumberType Number;
  dd_PolyhedraPtr parent;  /* pointing to the original polyhedra data */
  dd_rowrange m_alloc; /* allocated row size of matrix A */
  dd_colrange d_alloc; /* allocated col size of matrix A */

/* CONTROL: variables to control computation */
  dd_rowrange Iteration;
  dd_HalfspaceOrderType HalfspaceOrder;
  dd_RayPtr FirstRay, LastRay, ArtificialRay; /* The second description: Generator */
  dd_RayPtr PosHead, ZeroHead, NegHead, PosLast, ZeroLast, NegLast;
  dd_Adjacency **Edges;  /* adjacency relation storage for iteration k */
  unsigned int rseed;  /* random seed for random row permutation */

  boolean ColReduced;  /* flag to indicate that a column basis is computed and reduced */
  dd_bigrange LinearityDim;   
           /*  the dimension of the linearity space (when input is H), and
               the size of a minimal system of equations to determine the space (when V). */
  dd_colrange d_orig;  /* the size d of the original matrix A */
  dd_colindex newcol;  /* the size d of the original matrix A */
  
  dd_colindex InitialRayIndex;   /* InitialRayIndex[s] (s>=1) stores the corr. row index */
  dd_rowindex OrderVector;
  boolean RecomputeRowOrder;
  boolean PreOrderedRun;
  dd_rowset GroundSet, EqualitySet, NonequalitySet, 
       AddedHalfspaces, WeaklyAddedHalfspaces, InitialHalfspaces;
  long RayCount, FeasibleRayCount, WeaklyFeasibleRayCount,
       TotalRayCount, ZeroRayCount;
  long EdgeCount, TotalEdgeCount;
  long count_int,count_int_good,count_int_bad; /* no. of intersection operations */

  dd_Bmatrix B;
  dd_Bmatrix Bsave;  /* a copy of the dual basis inverse used to reduce the matrix A */

/* STATES: variables to represent current state. */
  dd_ErrorType Error;
  dd_CompStatusType CompStatus;  /* Computation Status */
  time_t starttime, endtime;
} dd_Cone;

/* Global Variables */
extern boolean debug;

/* end of cddtypes.h */
