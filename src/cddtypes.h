/* cddtypes.h: Header file for cddlib.c 
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.90e, July 12, 2000
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

#define COPYRIGHT   "Copyright (C) 1996, Komei Fukuda, fukuda@ifor.math.ethz.ch"
#define DDVERSION   "Version 0.90e (July 12, 2000)"
#include <time.h>

#define dd_wordlenmax     127
#define dd_linelenmax     255
#define dd_datawidth       10
#define dd_filenamelen    255

#define FALSE 0
#define TRUE 1

typedef int boolean;

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

typedef struct raydata *dd_RayPtr;

typedef struct raydata {
  mytype *Ray;
  dd_rowset ZeroSet;
  dd_rowrange FirstInfeasIndex;  /* the first inequality the ray violates */
  boolean feasible;  /* flag to store the feasibility */
  mytype ARay;   /* temporary area to store some row of A*Ray */
  dd_RayPtr Next;
} dd_RayType;

typedef struct adjacencydata *dd_AdjacencyPtr;
typedef struct adjacencydata {
  dd_RayPtr Ray1, Ray2;
  dd_AdjacencyPtr Next;
} dd_AdjacencyType;

typedef enum {
  Combinatorial, Algebraic
} dd_AdjacencyTestType;

typedef enum {
  MaxIndex, MinIndex, MinCutoff, MaxCutoff, MixCutoff,
   LexMin, LexMax, RandomRow
} dd_RowOrderType;

typedef enum {
  Real, Rational, Integer, Unknown
} dd_NumberType;

typedef enum {
  Inequality, Generator, Unspecified
} dd_RepresentationType;

typedef enum {
  IneToGen, GenToIne, LPMax, LPMin, InteriorFind
} dd_ConversionType;

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
  DimensionTooLarge, ImproperInputFormat, 
  NegativeMatrixSize, EmptyVrepresentation,
  IFileNotFound, OFileNotOpen, NoLPObjective, NoRealNumberSupport, NoError
} dd_ErrorType;

typedef enum {
  InProgress, AllFound, RegionEmpty
} dd_CompStatusType;

/* --- LP types ---- */

typedef enum {
  LPnone=0, LPmax, LPmin
} dd_LPObjectiveType;

typedef enum {
  CrissCross, DualSimplex
} dd_LPSolverType;

typedef enum {
  LPSundecided, Optimal, Inconsistent, DualInconsistent,
  StrucInconsistent, StrucDualInconsistent,
  Unbounded, DualUnbounded
} dd_LPStatusType;

typedef struct lpsolution *dd_LPSolutionPtr;
typedef struct lpsolution {
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
  long pivots[4]; 
   /* pivots[0]=setup (to find a basis), pivots[1]=PhaseI or Criss-Cross,
      pivots[2]=Phase II, pivots[3]=Anticycling */
  long total_pivots;
} dd_LPSolutionType;


typedef struct lpdata *dd_LPPtr;
typedef struct lpdata {
  dd_DataFileType filename;
  dd_LPObjectiveType objective;
  dd_LPSolverType solver; 
  boolean Homogeneous;  
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

  dd_LPStatusType LPS;  /* the current solution status */
  dd_rowrange m_alloc; /* the allocated row size of matrix A */
  dd_colrange d_alloc; /* the allocated col size of matrix A */
  mytype optvalue;  /* optimal value */
  dd_Arow sol;   /* primal solution */
  dd_Arow dsol;  /* dual solution */
  dd_colindex nbindex;  /* current basis represented by nonbasic indices */
  dd_rowrange re;  /* row index as a certificate in the case of inconsistency */
  dd_colrange se;  /* col index as a certificate in the case of dual inconsistency */
  long pivots[4]; 
   /* pivots[0]=setup (to find a basis), pivots[1]=PhaseI or Criss-Cross,
      pivots[2]=Phase II, pivots[3]=Anticycling */
  long total_pivots;
  int use_given_basis;  /* switch to indicate the use of the given basis */
  dd_colindex given_nbindex;  /* given basis represented by nonbasic indices */
  time_t starttime;
  time_t endtime;
} dd_LPType;


/*----  end of LP Types ----- */


typedef struct matrixdata *dd_MatrixPtr;
typedef struct matrixdata {
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
}  dd_MatrixType;

typedef struct setfamily *dd_SetFamilyPtr;
typedef struct setfamily {
  dd_bigrange famsize;
  dd_bigrange setsize;
  dd_SetVector set;  
} dd_SetFamilyType;


typedef struct nodedata *dd_NodePtr;
typedef struct nodedata {dd_bigrange key; dd_NodePtr next;} dd_NodeType;

typedef struct graphdata *dd_GraphPtr;
typedef struct graphdata {
  dd_bigrange vsize;
  dd_NodePtr *adjlist;  /* should be initialized to have vsize components */
} dd_GraphType;


typedef struct polyhedradata *dd_PolyhedraPtr;
typedef struct conedata *dd_ConePtr;

typedef struct polyhedradata {
  dd_RepresentationType representation;  /* given representation */
  boolean homogeneous;
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

  boolean NondegAssumed;
  boolean InitBasisAtBottom;
  boolean RestrictedEnumeration;
  boolean RelaxedEnumeration;

  dd_rowrange m1; 
    /* = m or m+1 (when representation=Inequality && !homogeneous)
       This data is written after ConeDataLoad is called.  This
       determines the size of Ainc. */
  boolean AincGenerated;
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


typedef struct conedata {
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
} dd_ConeType;

/* Global Variables */
extern boolean debug;

/* end of cddtypes.h */
