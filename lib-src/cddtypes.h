/* cddtypes.h: Header file for cddlib.c 
   written by Komei Fukuda, fukuda@math.ethz.ch
   Version 0.94h, April 30, 2015
*/

/* cddlib.c : C-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddlibman.tex for detail.
*/

#ifndef  __CDDTYPES_H
#define  __CDDTYPES_H
#endif  /* __CDDTYPES_H */

#define dd_COPYRIGHT   "Copyright (C) 1996, Komei Fukuda, fukuda@ifor.math.ethz.ch"
#define dd_DDVERSION   "Version 0.94g (March 23, 2012)"
#include <time.h>

#define dd_wordlenmax    1024 
#define dd_linelenmax    4096
#define dd_datawidth       10
#define dd_filenamelen    255

#define dd_FALSE 0
#define dd_TRUE 1

typedef int dd_boolean;

typedef long dd_rowrange;
typedef long dd_colrange;
typedef long dd_bigrange;

typedef set_type dd_rowset;
typedef set_type dd_colset;
typedef long *dd_rowindex;   
typedef int *dd_rowflag;   
typedef long *dd_colindex;
typedef mytype **dd_Amatrix;
typedef mytype *dd_Arow;
typedef set_type *dd_SetVector;
typedef mytype **dd_Bmatrix;
typedef set_type *dd_Aincidence;

/* typedef char dd_FilenameType[dd_filenamelen]; deleted 000505*/
typedef char dd_DataFileType[dd_filenamelen];
typedef char dd_LineType[dd_linelenmax];
typedef char dd_WordType[dd_wordlenmax];

typedef struct dd_raydata *dd_RayPtr;

typedef struct dd_raydata {
  mytype *Ray;
  dd_rowset ZeroSet;
  dd_rowrange FirstInfeasIndex;  /* the first inequality the ray violates */
  dd_boolean feasible;  /* flag to store the feasibility */
  mytype ARay;   /* temporary area to store some row of A*Ray */
  dd_RayPtr Next;
} dd_RayType;

typedef struct dd_adjacencydata *dd_AdjacencyPtr;
typedef struct dd_adjacencydata {
  dd_RayPtr Ray1, Ray2;
  dd_AdjacencyPtr Next;
} dd_AdjacencyType;

typedef enum {
  dd_Combinatorial, dd_Algebraic
} dd_AdjacencyTestType;

typedef enum {
  dd_MaxIndex, dd_MinIndex, dd_MinCutoff, dd_MaxCutoff, dd_MixCutoff,
   dd_LexMin, dd_LexMax, dd_RandomRow
} dd_RowOrderType;

typedef enum {
  dd_Unknown=0, dd_Real, dd_Rational, dd_Integer
} dd_NumberType;

typedef enum {
  dd_Unspecified=0, dd_Inequality, dd_Generator
} dd_RepresentationType;

typedef enum {
  dd_IneToGen, dd_GenToIne, dd_LPMax, dd_LPMin, dd_InteriorFind
} dd_ConversionType;

typedef enum {
  dd_IncOff=0, dd_IncCardinality, dd_IncSet
} dd_IncidenceOutputType;

typedef enum {
  dd_AdjOff=0, dd_AdjacencyList,  dd_AdjacencyDegree
} dd_AdjacencyOutputType;

typedef enum {
  dd_Auto, dd_SemiAuto, dd_Manual
} dd_FileInputModeType;   
   /* Auto if a input filename is specified by command arguments */

typedef enum {
  dd_DimensionTooLarge, dd_ImproperInputFormat, 
  dd_NegativeMatrixSize, dd_EmptyVrepresentation, dd_EmptyHrepresentation, dd_EmptyRepresentation,
  dd_IFileNotFound, dd_OFileNotOpen, dd_NoLPObjective, dd_NoRealNumberSupport,
  dd_NotAvailForH, dd_NotAvailForV, dd_CannotHandleLinearity,
  dd_RowIndexOutOfRange, dd_ColIndexOutOfRange,
  dd_LPCycling, dd_NumericallyInconsistent,
  dd_NoError
} dd_ErrorType;

typedef enum {
  dd_InProgress, dd_AllFound, dd_RegionEmpty
} dd_CompStatusType;

/* --- LP types ---- */

typedef enum {
  dd_LPnone=0, dd_LPmax, dd_LPmin
} dd_LPObjectiveType;

typedef enum {
  dd_CrissCross, dd_DualSimplex
} dd_LPSolverType;

typedef enum {
  dd_LPSundecided, dd_Optimal, dd_Inconsistent, dd_DualInconsistent,
  dd_StrucInconsistent, dd_StrucDualInconsistent,
  dd_Unbounded, dd_DualUnbounded
} dd_LPStatusType;

typedef struct dd_lpsolution *dd_LPSolutionPtr;
typedef struct dd_lpsolution {
  dd_DataFileType filename;
  dd_LPObjectiveType objective;
  dd_LPSolverType solver; 
  dd_rowrange m;
  dd_colrange d;
  dd_NumberType numbtype;

  dd_LPStatusType LPS;  /* the current solution status */
  mytype optvalue;  /* optimal value */
  dd_Arow sol;   /* primal solution */
  dd_Arow dsol;  /* dual solution */
  dd_colindex nbindex;  /* current basis represented by nonbasic indices */
  dd_rowrange re;  /* row index as a certificate in the case of inconsistency */
  dd_colrange se;  /* col index as a certificate in the case of dual inconsistency */
  long pivots[5]; 
   /* pivots[0]=setup (to find a basis), pivots[1]=PhaseI or Criss-Cross,
      pivots[2]=Phase II, pivots[3]=Anticycling, pivots[4]=GMP postopt. */
  long total_pivots;
} dd_LPSolutionType;


typedef struct dd_lpdata *dd_LPPtr;
typedef struct dd_lpdata {
  dd_DataFileType filename;
  dd_LPObjectiveType objective;
  dd_LPSolverType solver; 
  dd_boolean Homogeneous;  
     /* The first column except for the obj row is all zeros. */
  dd_rowrange m;
  dd_colrange d;
  dd_Amatrix A;
  dd_Bmatrix B;
  dd_rowrange objrow;
  dd_colrange rhscol;
  dd_NumberType numbtype;
  dd_rowrange eqnumber;  /* the number of equalities */
  dd_rowset equalityset;  

  dd_boolean redcheck_extensive;  /* Apply the extensive redundancy check. */
  dd_rowrange ired; /* the row index for the redundancy checking */
  dd_rowset redset_extra;  /* a set of rows that are newly recognized redundan by the extensive search. */
  dd_rowset redset_accum;  /* the accumulated set of rows that are recognized redundant */
  dd_rowset posset_extra;  /* a set of rows that are recognized non-linearity */

  dd_boolean lexicopivot;  /* flag to use the lexicogrphic pivot rule (symbolic perturbation). */

  dd_LPStatusType LPS;  /* the current solution status */
  dd_rowrange m_alloc; /* the allocated row size of matrix A */
  dd_colrange d_alloc; /* the allocated col size of matrix A */
  mytype optvalue;  /* optimal value */
  dd_Arow sol;   /* primal solution */
  dd_Arow dsol;  /* dual solution */
  dd_colindex nbindex;  /* current basis represented by nonbasic indices */
  dd_rowrange re;  /* row index as a certificate in the case of inconsistency */
  dd_colrange se;  /* col index as a certificate in the case of dual inconsistency */
  long pivots[5]; 
   /* pivots[0]=setup (to find a basis), pivots[1]=PhaseI or Criss-Cross,
      pivots[2]=Phase II, pivots[3]=Anticycling, pivots[4]=GMP postopt. */
  long total_pivots;
  int use_given_basis;  /* switch to indicate the use of the given basis */
  dd_colindex given_nbindex;  /* given basis represented by nonbasic indices */
  time_t starttime;
  time_t endtime;
} dd_LPType;


/*----  end of LP Types ----- */


typedef struct  dd_matrixdata *dd_MatrixPtr;
typedef struct  dd_matrixdata {
  dd_rowrange rowsize;
  dd_rowset linset; 
    /*  a subset of rows of linearity (ie, generators of
        linearity space for V-representation, and equations
        for H-representation. */
  dd_colrange colsize;
  dd_RepresentationType representation;
  dd_NumberType numbtype;
  dd_Amatrix matrix;
  dd_LPObjectiveType objective;
  dd_Arow rowvec;
} dd_MatrixType;

typedef struct dd_setfamily *dd_SetFamilyPtr;
typedef struct dd_setfamily {
  dd_bigrange famsize;
  dd_bigrange setsize;
  dd_SetVector set;  
} dd_SetFamilyType;


typedef struct dd_nodedata *dd_NodePtr;
typedef struct dd_nodedata {dd_bigrange key; dd_NodePtr next;} dd_NodeType;

typedef struct dd_graphdata *dd_GraphPtr;
typedef struct dd_graphdata {
  dd_bigrange vsize;
  dd_NodePtr *adjlist;  /* should be initialized to have vsize components */
} dd_GraphType;


typedef struct dd_polyhedradata *dd_PolyhedraPtr;
typedef struct dd_conedata *dd_ConePtr;

typedef struct dd_polyhedradata {
  dd_RepresentationType representation;  /* given representation */
  dd_boolean homogeneous;
  dd_colrange d;
  dd_rowrange m;
  dd_Amatrix A;   /* Inequality System:  m times d matrix */
  dd_NumberType numbtype;
  dd_ConePtr child;  /* pointing to the homogenized cone data */
  dd_rowrange m_alloc; /* allocated row size of matrix A */
  dd_colrange d_alloc; /* allocated col size of matrix A */
  dd_Arow c;           /* cost vector */

  dd_rowflag EqualityIndex;  
    /* ith component is 1 if it is equality, -1 if it is strict inequality, 0 otherwise. */

  dd_boolean IsEmpty;  /* This is to tell whether the set is empty or not */
  
  dd_boolean NondegAssumed;
  dd_boolean InitBasisAtBottom;
  dd_boolean RestrictedEnumeration;
  dd_boolean RelaxedEnumeration;

  dd_rowrange m1; 
    /* = m or m+1 (when representation=Inequality && !homogeneous)
       This data is written after dd_ConeDataLoad is called.  This
       determines the size of Ainc. */
  dd_boolean AincGenerated;
    /* Indicates whether Ainc, Ared, Adom are all computed.
       All the variables below are valid only when this is TRUE */
  dd_colrange ldim;   /* linearity dimension */
  dd_bigrange n; 
    /* the size of output = total number of rays 
       in the computed cone + linearity dimension */
  dd_Aincidence Ainc;
    /* incidence of input and output */
  dd_rowset Ared;  
    /* redundant set of rows whose removal results in a minimal system */
  dd_rowset Adom;  
    /* dominant set of rows (those containing all rays). */

} dd_PolyhedraType;


typedef struct dd_conedata {
  dd_RepresentationType representation;
  dd_rowrange m;
  dd_colrange d;
  dd_Amatrix A;
  dd_NumberType numbtype;
  dd_PolyhedraPtr parent;  /* pointing to the original polyhedra data */
  dd_rowrange m_alloc; /* allocated row size of matrix A */
  dd_colrange d_alloc; /* allocated col size of matrix A */

/* CONTROL: variables to control computation */
  dd_rowrange Iteration;
  dd_RowOrderType HalfspaceOrder;
  dd_RayPtr FirstRay, LastRay, ArtificialRay; /* The second description: Generator */
  dd_RayPtr PosHead, ZeroHead, NegHead, PosLast, ZeroLast, NegLast;
  dd_AdjacencyType **Edges;  /* adjacency relation storage for iteration k */
  unsigned int rseed;  /* random seed for random row permutation */

  dd_boolean ColReduced;  /* flag to indicate that a column basis is computed and reduced */
  dd_bigrange LinearityDim;   
           /*  the dimension of the linearity space (when input is H), and
               the size of a minimal system of equations to determine the space (when V). */
  dd_colrange d_orig;  /* the size d of the original matrix A */
  dd_colindex newcol;  /* the size d of the original matrix A */
  
  dd_colindex InitialRayIndex;   /* InitialRayIndex[s] (s>=1) stores the corr. row index */
  dd_rowindex OrderVector;
  dd_boolean RecomputeRowOrder;
  dd_boolean PreOrderedRun;
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
} dd_ConeType;

/* Global Variables */
extern dd_boolean dd_debug;
extern dd_boolean dd_log;

/* end of cddtypes.h */
