/* cddio.c:  Basic Input and Output Procedures for cddlib.c
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.90, May 28, 2000
*/

/* cdd.c : C-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#include "setoper.h"  /* set operation library header (Ver. March 16,1995 or later) */
#include "cdd.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

/* void fread_rational_value (FILE *, mytype *); */
void SetLinearity(dd_MatrixPtr, char *);

void dd_SetInputFile(FILE **f,dd_DataFileType inputfile,dd_ErrorType *Error)
{
  int opened=0,stop,quit=0;
  int i,dotpos=0,trial=0;
  char ch;
  char *tempname;
  
  
  *Error=NoError;
  while (!opened && !quit) {
    printf("\n>> Input file: ");
    scanf("%s",inputfile);
    ch=getchar();
    stop=FALSE;
    for (i=0; i<dd_filenamelen && !stop; i++){
      ch=inputfile[i];
      switch (ch) {
        case '.': 
          dotpos=i+1;
          break;
        case ';':  case ' ':  case '\0':  case '\n':  case '\t':     
          stop=TRUE;
          tempname=(char*)calloc(dd_filenamelen,sizeof(ch));
          strncpy(tempname,inputfile,i);
          strcpy(inputfile,tempname);
          break;
      }
    }
    if ( ( *f = fopen(inputfile,"r") )!= NULL) {
      printf("input file %s is open\n",inputfile);
      opened=1;
      *Error=NoError;
    }
    else{
      printf("The file %s not found\n",inputfile);
      trial++;
      if (trial>5) {
        *Error=IFileNotFound;
        quit=1;
      }
    }
  }
}

void dd_SetWriteFileName(dd_DataFileType inputfile, dd_DataFileType outfile, char cflag, dd_RepresentationType rep)
{
  char *extension;
  dd_DataFileType ifilehead="";
  int i,dotpos;
  
  switch (cflag) {
    case 'o':
      switch (rep) {
        case Generator:
          extension=".ine"; break;     /* output file for ine data */
        case Inequality:
          extension=".ext"; break;     /* output file for ext data */
        default: 
          extension=".xxx";break;
      }
      break;

    case 'a':         /* decide for output adjacence */
      if (rep==Inequality)
        extension=".ead";       /* adjacency file for ext data */
      else
        extension=".iad";       /* adjacency file for ine data */
      break;
    case 'i':         /* decide for output incidence */
      if (rep==Inequality)
        extension=".ecd";       /* ext incidence file */
      else
        extension=".icd";       /* ine incidence file */
      break;
    case 'n':         /* decide for input incidence */
      if (rep==Inequality)
        extension=".icd";       /* ine incidence file */
      else
        extension=".ecd";       /* ext incidence file */
      break;
    case 'j':        /* decide for input adjacency */
      if (rep==Inequality)
        extension=".iad";       /* ine adjacency file */
      else
        extension=".ead";       /* ext adjacency file */
      break;
    case 'l':
      extension=".ddl";break;   /* log file */
    case 'd':
      extension=".dex";break;   /* decomposition output */
    case 'p':
      extension="sub.ine";break;  /* preprojection sub inequality file */
    case 'v':
      extension=".solved";break;  /* verify_input file */
    case 's':
      extension=".lps";break;   /* LP solution file */
    default:
      extension=".xxx";break;
  }
  dotpos=-1;
  for (i=0; i< strlen(inputfile); i++){
    if (inputfile[i]=='.') dotpos=i;
  }
  if (dotpos>1) strncpy(ifilehead, inputfile, dotpos);
  else strcpy(ifilehead,inputfile);
  if (strlen(inputfile)<=0) strcpy(ifilehead,"tempcdd");
  strcpy(outfile,ifilehead); 
  strcat(outfile,extension); 
  if (strcmp(inputfile, outfile)==0) {
    strcpy(outfile,inputfile); 
    strcat(outfile,extension); 
  }
/*  fprintf(stdout,"outfile name = %s\n",outfile);  */
}


dd_NumberType GetNumberType(char *line)
{
  dd_NumberType nt;

  if (strncmp(line, "integer", 7)==0) {
    nt = Integer;
  }
  else if (strncmp(line, "rational", 8)==0) {
    nt = Rational;
  }
  else if (strncmp(line, "real", 4)==0) {
    nt = Real;
  }
  else { 
    nt=Unknown;
  }
  return nt;
}

void ProcessCommandLine(FILE *f, dd_MatrixPtr M, char *line)
{
  char newline[dd_linelenmax];
  dd_colrange j;
  mytype value;

  dd_init(value);
  if (strncmp(line, "hull", 4)==0) {
    M->representation = Generator;
  }
  if (strncmp(line, "debug", 5)==0) {
    debug = TRUE;
  }
  if (strncmp(line, "partial_enum", 12)==0 ||
       strncmp(line, "equality", 8)==0  ||
       strncmp(line, "linearity", 9)==0 ) {
    fgets(newline,dd_linelenmax,f);    
    SetLinearity(M,newline);
  }
  if (strncmp(line, "maximize", 8)==0 ||
      strncmp(line, "minimize", 8)==0) {
    if (strncmp(line, "maximize", 8)==0) M->objective=LPmax;
    else M->objective=LPmin;
    for (j = 1; j <= M->colsize; j++) {
    if (M->numbtype==Real) {
#if !defined(GMPRATIONAL)
        double rvalue;
        fscanf(f, "%lf", &rvalue);
        dd_set_d(value, rvalue);
#endif
      } else {
        fread_rational_value (f, value);
      }
      dd_set(M->rowvec[j - 1],value);
      if (debug) {printf("cost(%5ld) =",j); dd_WriteNumber(stdout,value);}
    }  /*of j*/
  }
  dd_clear(value);
}

void AddInequalities(dd_PolyhedraPtr poly, dd_MatrixPtr M)
{
  dd_Amatrix Anew;
  dd_rowrange i, m_new;
  dd_colrange j, d_new;

  /* poly->child->CompStatus=InProgress;  086 */
  
  m_new=poly->m + M->rowsize;
  d_new=poly->d_alloc;

  if (poly->m_alloc < m_new){
    dd_InitializeAmatrix(m_new,d_new,&(Anew));
    CopyAmatrix(Anew, poly->A, poly->m, poly->d);
    dd_FreeAmatrix(poly->m_alloc,poly->d_alloc,poly->A);
    poly->A=Anew;
    poly->m_alloc=m_new;
  }
  for (i=0; i<M->rowsize; i++){
    for (j=0; j<poly->d; j++) dd_set(poly->A[(poly->m)+i][j],M->matrix[i][j]);
  }
  poly->m=m_new;  
}


dd_PolyhedraPtr CreatePolyhedraData(dd_rowrange m, dd_colrange d)
{
  dd_rowrange i;
  dd_PolyhedraPtr poly=NULL;

  poly=(dd_PolyhedraPtr) malloc (sizeof(dd_PolyhedraType));
  poly->child       =NULL; /* this links the homogenized cone data */
  poly->m           =m;
  poly->d           =d;  
  poly->m_alloc     =m+2; /* the allocated row size of matrix A */
  poly->d_alloc     =d;   /* the allocated col size of matrix A */
  poly->numbtype=Real;
  dd_InitializeAmatrix(poly->m_alloc,poly->d_alloc,&(poly->A));
  dd_InitializeArow(d,&(poly->c));           /* cost vector */
  poly->representation       =Inequality;
  poly->homogeneous =FALSE;

  poly->EqualityIndex=(int *)calloc(m+1, sizeof(int));  
    /* ith component is 1 if it is equality, -1 if it is strict inequality, 0 otherwise. */
  for (i = 0; i <= m; i++) poly->EqualityIndex[i]=0;

  poly->NondegAssumed           = FALSE;
  poly->InitBasisAtBottom       = FALSE;
  poly->RestrictedEnumeration   = FALSE;
  poly->RelaxedEnumeration      = FALSE;

  return poly;
}

boolean InitializeConeData(dd_rowrange m, dd_colrange d, dd_ConePtr *cone)
{
  boolean success=TRUE;
  dd_colrange j;

  (*cone)=(dd_ConePtr)calloc(1, sizeof(dd_ConeType));

/* INPUT: A given representation of a cone: inequality */
  (*cone)->m=m;
  (*cone)->d=d;
  (*cone)->m_alloc=m+2; /* allocated row size of matrix A */
  (*cone)->d_alloc=d;   /* allocated col size of matrix A */
  (*cone)->numbtype=Real;
  (*cone)->parent=NULL;

/* CONTROL: variables to control computation */
  (*cone)->Iteration=0;

  (*cone)->HalfspaceOrder=LexMin;

  (*cone)->ArtificialRay=NULL;
  (*cone)->FirstRay=NULL;
  (*cone)->LastRay=NULL; /* The second description: Generator */
  (*cone)->PosHead=NULL;
  (*cone)->ZeroHead=NULL;
  (*cone)->NegHead=NULL;
  (*cone)->PosLast=NULL;
  (*cone)->ZeroLast=NULL;
  (*cone)->NegLast=NULL;
  (*cone)->RecomputeRowOrder  = TRUE;
  (*cone)->PreOrderedRun      = FALSE;
  set_initialize(&((*cone)->GroundSet),(*cone)->m_alloc);
  set_initialize(&((*cone)->EqualitySet),(*cone)->m_alloc);
  set_initialize(&((*cone)->NonequalitySet),(*cone)->m_alloc);
  set_initialize(&((*cone)->AddedHalfspaces),(*cone)->m_alloc);
  set_initialize(&((*cone)->WeaklyAddedHalfspaces),(*cone)->m_alloc);
  set_initialize(&((*cone)->InitialHalfspaces),(*cone)->m_alloc);
  (*cone)->RayCount=0;
  (*cone)->FeasibleRayCount=0;
  (*cone)->WeaklyFeasibleRayCount=0;
  (*cone)->TotalRayCount=0;
  (*cone)->ZeroRayCount=0;
  (*cone)->EdgeCount=0;
  (*cone)->TotalEdgeCount=0;
  (*cone)->count_int=0;
  (*cone)->count_int_good=0;
  (*cone)->count_int_bad=0;
  (*cone)->rseed=1;  /* random seed for random row permutation */
 
  dd_InitializeBmatrix((*cone)->d, &((*cone)->B));
  dd_InitializeBmatrix((*cone)->d, &((*cone)->Bsave));
  dd_InitializeAmatrix((*cone)->m_alloc,(*cone)->d_alloc,&((*cone)->A));
  (*cone)->AicdGenerated=FALSE;  /* Aicd is a set array to store the input incidence. */

  (*cone)->Edges
     =(dd_AdjacencyType**) calloc((*cone)->m_alloc,sizeof(dd_AdjacencyType*));
  (*cone)->InitialRayIndex=(long*)calloc(d+1,sizeof(long));
  (*cone)->OrderVector=(long*)calloc((*cone)->m_alloc+1,sizeof(long));


  (*cone)->newcol=(long*)calloc(((*cone)->d)+1,sizeof(long));
  for (j=0; j<=(*cone)->d; j++) (*cone)->newcol[j]=j;  /* identity map, initially */
  (*cone)->LinearityDim = -2; /* -2 if it is not computed */
  (*cone)->ColReduced   = FALSE;
  (*cone)->d_orig = d;

/* STATES: variables to represent current state. */
/*(*cone)->Error;
  (*cone)->CompStatus;
  (*cone)->starttime;
  (*cone)->endtime;
*/
    
  return success;
}

dd_ConePtr ConeDataLoad(dd_PolyhedraPtr poly)
{
  dd_ConePtr cone=NULL;
  dd_colrange d,j;
  dd_rowrange m,i;

  m=poly->m;
  d=poly->d;
  if (!(poly->homogeneous) && poly->representation==Inequality) m=poly->m+1;

  InitializeConeData(m, d, &cone);
  cone->representation=poly->representation;

/* Points to the original polyhedra data, and reversely */
  cone->parent=poly;
  poly->child=cone;

  for (i=1; i<=poly->m; i++)
    for (j=1; j<=cone->d; j++)
      dd_set(cone->A[i-1][j-1],poly->A[i-1][j-1]);  
  
  if (poly->representation==Inequality && !(poly->homogeneous)){
    dd_set(cone->A[m-1][0],dd_one);
    for (j=2; j<=d; j++) dd_set(cone->A[m-1][j-1],dd_purezero);
  }

  return cone;
}

void SetLinearity(dd_MatrixPtr M, char *line)
{
  int i=0;
  dd_rowrange eqsize,var;
  char *next;
  const char ct[]=", ";  /* allows separators "," and " ". */

  next=strtok(line,ct);
  eqsize=atol(next); 
  while (i < eqsize && (next=strtok(NULL,ct))!=NULL) {
     var=atol(next);
     set_addelem(M->linset,var); i++;
  }
  if (i!=eqsize) {
    printf("* Warning: there are inconsistencies in linearity setting.\n");
  }
  return;
}

dd_MatrixPtr dd_PolyFile2Matrix (FILE *f, dd_ErrorType *Error)
{
  dd_MatrixPtr M=NULL;
  dd_rowrange m_input,i;
  dd_colrange d_input,j;
  dd_RepresentationType rep=Inequality;
  mytype value;
  boolean found=FALSE, newformat=FALSE, successful=FALSE, linearity=FALSE;
  char command[dd_linelenmax], comsave[dd_linelenmax], numbtype[dd_wordlenmax];
  dd_NumberType NT;
#if !defined(GMPRATIONAL)
  double rvalue;
#endif

  dd_init(value);
  (*Error)=NoError;
  while (!found)
  {
    if (fscanf(f,"%s",command)==EOF) {
      (*Error)=ImproperInputFormat;
      goto _L99;
    }
    else {
      if (strncmp(command, "V-representation", 16)==0) {
        rep=Generator; newformat=TRUE;
      }
      if (strncmp(command, "H-representation", 16)==0){
        rep=Inequality; newformat=TRUE;
      }
      if (strncmp(command, "partial_enum", 12)==0 || 
          strncmp(command, "equality", 8)==0  ||
          strncmp(command, "linearity", 9)==0 ) {
        linearity=TRUE;
        fgets(comsave,dd_linelenmax,f);
      }
      if (strncmp(command, "begin", 5)==0) found=TRUE;
    }
  }
  fscanf(f, "%ld %ld %s", &m_input, &d_input, numbtype);
  printf("size = %ld x %ld\nNumber Type = %s\n", m_input, d_input, numbtype);
  NT=GetNumberType(numbtype);
  if (NT==Unknown) {
      (*Error)=ImproperInputFormat;
      goto _L99;
    } 
  M=dd_CreateMatrix(m_input, d_input);
  M->representation=rep;
  M->numbtype=NT;

  for (i = 1; i <= m_input; i++) {
    for (j = 1; j <= d_input; j++) {
      if (NT==Real) {
#if defined GMPRATIONAL
        *Error=NoRealNumberSupport;
        goto _L99;
#else
        fscanf(f, "%lf", &rvalue);
        dd_set_d(value, rvalue);
#endif
      } else {
        fread_rational_value (f, value);
      }
      dd_set(M->matrix[i-1][j - 1],value);
      if (debug) {printf("a(%3ld,%5ld) = ",i,j); dd_WriteNumber(stdout,value);}
    }  /*of j*/
  }  /*of i*/
  if (fscanf(f,"%s",command)==EOF) {
   	 (*Error)=ImproperInputFormat;
  	 goto _L99;
  }
  else if (strncmp(command, "end", 3)!=0) {
     if (debug) printf("'end' missing or illegal extra data: %s\n",command);
     (*Error)=ImproperInputFormat;
     goto _L99;
  }
  
  successful=TRUE;
  if (linearity) {
    SetLinearity(M,comsave);
  }
  while (!feof(f)) {
    fscanf(f,"%s", command);
    ProcessCommandLine(f, M, command);
    fgets(command,dd_linelenmax,f); /* skip the CR/LF */
  } 

_L99: ;
  dd_clear(value);
  /* if (f!=NULL) fclose(f); */
  return M;
}


dd_PolyhedraPtr dd_Matrix2Poly(dd_MatrixPtr M, dd_ErrorType *err)
{
  dd_rowrange i;
  dd_colrange j;
  dd_PolyhedraPtr poly;

  *err=NoError;
  poly=CreatePolyhedraData(M->rowsize, M->colsize);
  poly->representation=M->representation;
  poly->homogeneous=TRUE;

  for (i = 1; i <= M->rowsize; i++) {
    if (set_member(i, M->linset)) {
      poly->EqualityIndex[i]=1;
    }
    for (j = 1; j <= M->colsize; j++) {
      dd_set(poly->A[i-1][j-1], M->matrix[i-1][j-1]);
      if (j==1 && dd_Nonzero(M->matrix[i-1][j-1])) poly->homogeneous = FALSE;
    }  /*of j*/
  }  /*of i*/
  return poly; 
}

void MatrixIntegerFilter(dd_MatrixPtr M)
{   /* setting an almost integer to the integer. */
  dd_rowrange i;
  dd_colrange j;
  mytype x;

  dd_init(x);
  for (i=0; i< M->rowsize; i++)
    for (j=0; j< M->colsize; j++){
       SnapToInteger(x, M->matrix[i][j]);
       dd_set(M->matrix[i][j],x);
    }
  dd_clear(x);
}

void CopyRay(mytype *a, dd_colrange d_origsize, dd_RayPtr RR, 
  dd_RepresentationType rep, dd_colindex reducedcol)
{
  long j,j1;
  mytype b;

  dd_init(b);
  for (j = 1; j <= d_origsize; j++){
    j1=reducedcol[j];
    if (j1>0){
      dd_set(a[j-1],RR->Ray[j1-1]); 
        /* the original column j is mapped to j1, and thus
           copy the corresponding component */
    } else {
      dd_set(a[j-1],dd_purezero);  
        /* original column is redundant and removed for computation */
    }
  }

  dd_set(b,a[0]);
  if (rep==Generator && dd_Nonzero(b)){
    dd_set(a[0],dd_one);
    for (j = 2; j <= d_origsize; j++)
       dd_div(a[j-1],a[j-1],b);    /* normalization for generators */
  }
  dd_clear(b);
}

void WriteRay(FILE *f, dd_colrange d_origsize, dd_RayPtr RR, dd_RepresentationType rep, dd_colindex reducedcol)
{
  dd_colrange j;
  static dd_colrange d_last=0;
  static dd_Arow a;

  if (d_last< d_origsize){
    if (d_last>0) free(a);
    dd_InitializeArow(d_origsize+1, &a);
    d_last=d_origsize+1;
  }

  CopyRay(a, d_origsize, RR, rep, reducedcol);
  for (j = 0; j < d_origsize; j++) dd_WriteNumber(f, a[j]);
  fprintf(f, "\n");
}

void WriteArow(FILE *f, dd_colrange d, dd_Arow a)
{
  dd_colrange j;

  for (j = 0; j < d; j++) dd_WriteNumber(f, a[j]);
  fprintf(f, "\n");
}



void dd_WriteAmatrix(FILE *f, dd_Amatrix A, long rowmax, long colmax)
{
  long i,j;

  if (A==NULL){
    fprintf(f, "WriteAmatrix: The requested matrix is empty\n");
    goto _L99;
  }
  fprintf(f, "begin\n");
#if defined GMPRATIONAL
  fprintf(f, " %ld %ld rational\n",rowmax, colmax);
#else
  fprintf(f, " %ld %ld real\n",rowmax, colmax);
#endif
  for (i=1; i <= rowmax; i++) {
    for (j=1; j <= colmax; j++) {
      dd_WriteNumber(f, A[i-1][j-1]);
    }
    fprintf(f,"\n");
  }
  fprintf(f, "end\n");
_L99:;
}

void dd_WriteBmatrix(FILE *f, dd_colrange d_size, dd_Bmatrix B)
{
  dd_colrange j1, j2;

  if (B==NULL){
    fprintf(f, "WriteBmatrix: The requested matrix is empty\n");
    goto _L99;
  }
  for (j1 = 0; j1 < d_size; j1++) {
    for (j2 = 0; j2 < d_size; j2++) {
      dd_WriteNumber(f, B[j1][j2]);
    }  /*of j2*/
    putc('\n', f);
  }  /*of j1*/
  putc('\n', f);
_L99:;
}

void dd_WriteSetFamily(FILE *f, dd_SetFamilyPtr F)
{
  dd_bigrange i;

  if (F==NULL){
    fprintf(f, "WriteSetFamily: The requested family is empty\n");
    goto _L99;
  }
  fprintf(f,"begin\n");
  fprintf(f,"  %ld    %ld\n", F->famsize, F->setsize);
  for (i=0; i<F->famsize; i++) {
    fprintf(f, " %ld %ld : ", i+1, set_card(F->set[i]));
    set_fwrite(f, F->set[i]);
  }
  fprintf(f,"end\n");
_L99:;
}

void dd_WriteSetFamilyCompressed(FILE *f, dd_SetFamilyPtr F)
{
  dd_bigrange i,card;

  if (F==NULL){
    fprintf(f, "WriteSetFamily: The requested family is empty\n");
    goto _L99;
  }
  fprintf(f,"begin\n");
  fprintf(f,"  %ld    %ld\n", F->famsize, F->setsize);
  for (i=0; i<F->famsize; i++) {
    card=set_card(F->set[i]);
    if (F->setsize - card >= card){
      fprintf(f, " %ld %ld : ", i+1, card);
      set_fwrite(f, F->set[i]);
    } else {
      fprintf(f, " %ld %ld : ", i+1, -card);
      set_fwrite_compl(f, F->set[i]);
    }
  }
  fprintf(f,"end\n");
_L99:;
}

void dd_WriteMatrix(FILE *f, dd_MatrixPtr M)
{
  dd_rowrange i, linsize;

  if (M==NULL){
    fprintf(f, "WriteMatrix: The requested matrix is empty\n");
    goto _L99;
  }
  switch (M->representation) {
    case Inequality:
      fprintf(f, "H-representation\n"); break;
    case Generator:
      fprintf(f, "V-representation\n"); break;
    case Unspecified:
      break;
  }
  linsize=set_card(M->linset);
  if (linsize>0) {
    fprintf(f, "linearity %ld ", linsize);
    for (i=1; i<=M->rowsize; i++) 
      if (set_member(i, M->linset)) fprintf(f, " %ld", i);
    fprintf(f, "\n");
  }
  dd_WriteAmatrix(f, M->matrix, M->rowsize, M->colsize);
_L99:;
}

void SnapToInteger(mytype y, mytype x)
{
 /* this is broken.  It does nothing. */
   dd_set(y,x);
}


void dd_WriteReal(FILE *f, mytype x)
{
  long ix1,ix2,ix;
  double ax;

  ax=dd_get_d(x);  
  ix1= fabs(ax) * 10000. + 0.5;
  ix2= (fabs(ax) + 0.5);
  ix2= ix2*10000;
  if ( ix1 == ix2) {
    if (dd_Positive(x)) {
      ix = ax + 0.5;
    } else {
      ix = -ax + 0.5;
      ix = -ix;
    }
    fprintf(f, " %2ld", ix);
  } else
    fprintf(f, " % .9E",ax);
}

void dd_WriteNumber(FILE *f, mytype x)
{
#if defined GMPRATIONAL
  mpz_t zn,zd;

  mpz_init(zn); mpz_init(zd);
  mpq_canonicalize(x);
  mpq_get_num(zn,x);
  mpq_get_den(zd,x);
  fprintf(f," ");
  if (mpz_sgn(zn)==0){
    fprintf(f,"0");
  } else if (mpz_cmp_ui(zd,1U)==0){
    mpz_out_str(f,10,zn);
  } else {
    mpz_out_str(f,10,zn);fprintf(f,"/");mpz_out_str(f,10,zd);
  }
  mpz_clear(zn); mpz_clear(zd);
#else
  dd_WriteReal(f, x);
#endif
}


void dd_WriteIncidence(FILE *f, dd_PolyhedraPtr poly)
{
  dd_SetFamilyPtr I;

  switch (poly->representation) {
    case Inequality:
       fprintf(f, "ecd_file: Incidence of generators and inequalities\n");
      break;
    case Generator:
       fprintf(f, "icd_file: Incidence of inequalities and generators\n");
      break;

    default:
      break;
  }
  I=dd_CopyIncidence(poly);
  dd_WriteSetFamilyCompressed(f,I);
  dd_FreeSetFamily(I);
}

void dd_WriteAdjacency(FILE *f, dd_PolyhedraPtr poly)
{
  dd_SetFamilyPtr A;

  switch (poly->representation) {
    case Inequality:
       fprintf(f, "ead_file: Adjacency of generators\n");
      break;
    case Generator:
       fprintf(f, "iad_file: Adjacency of inequalities\n");
      break;

    default:
      break;
  }
  A=dd_CopyAdjacency(poly);
  dd_WriteSetFamilyCompressed(f,A);
  dd_FreeSetFamily(A);
}


void ComputeAicd(dd_ConePtr cone)
{
  dd_RayPtr RayPtr1;
  dd_rowrange i,k, elem;
  boolean redundant;
  long pos1, scard;

  if (cone->AicdGenerated==FALSE){
    cone->Aicd=(set_type*)calloc(cone->m_alloc, sizeof(set_type));
    for(i=1; i<=cone->m_alloc; i++) set_initialize(&(cone->Aicd[i-1]),cone->RayCount);
    set_initialize(&(cone->Ared), cone->m); 
    set_initialize(&(cone->Adom), cone->m); 
  }
  if (cone->RayCount==0){
    goto _L99;
  }
  if (cone->AicdGenerated==FALSE){
    cone->LastRay->Next=NULL;
    for (RayPtr1=cone->FirstRay, pos1=1;RayPtr1 != NULL; RayPtr1 = RayPtr1->Next, pos1++){
      scard=set_card(RayPtr1->ZeroSet);
      elem = 0; i = 0;
      while (elem < scard && i < cone->m){
        i++;
        if (set_member(i,RayPtr1->ZeroSet)) {
          set_addelem(cone->Aicd[i-1],pos1);
          elem++;
       }
      }
    }
  }
  for (i=1; i<=cone->m; i++){
    if (set_card(cone->Aicd[i-1])==cone->RayCount){
      set_addelem(cone->Adom, i);
    }  
  }
  for (i=cone->m; i>=1; i--){
    if (set_card(cone->Aicd[i-1])==0){
      redundant=TRUE;
      set_addelem(cone->Ared, i);
    }else {
      redundant=FALSE;
      for (k=1; k<=cone->m; k++) {
        if (k!=i && !set_member(k, cone->Ared)  && !set_member(k, cone->Adom) && 
            set_subset(cone->Aicd[i-1], cone->Aicd[k-1])){
          if (!redundant){
            redundant=TRUE;
          }
          set_addelem(cone->Ared, i);
        }
      }
    }
  }
_L99:;
  cone->AicdGenerated=TRUE;
}

boolean InputAdjacentQ(dd_ConePtr cone, 
  dd_rowrange i1, dd_rowrange i2)
/* Before calling this function, RedundantSet must be 
   a set of row indices whose removal results in a minimal
   nonredundant system to represent the input polyhedron,
   DominantSet must be the set of row indices which are
   active at every extreme points/rays.
*/
{
  boolean adj=TRUE;
  dd_rowrange i;
  static set_type common;
  static long lastRayCount=0;

  if (cone->AicdGenerated==FALSE) ComputeAicd(cone);
  if (lastRayCount!=cone->RayCount){
    if (lastRayCount >0) set_free(common);
    set_initialize(&common, cone->RayCount);
    lastRayCount=cone->RayCount;
  }
  if (set_member(i1, cone->Ared) || set_member(i2, cone->Ared)){
    adj=FALSE;
    goto _L99;
  }
  if (set_member(i1, cone->Adom) || set_member(i2, cone->Adom)){
  // dominant inequality is considered adjacencent to all others.
    adj=TRUE;
    goto _L99;
  }
  set_int(common, cone->Aicd[i1-1], cone->Aicd[i2-1]);
  i=0;
  while (i<cone->m && adj==TRUE){ 
    i++; 
    if (i!=i1 && i!=i2 && !set_member(i, cone->Ared) &&
        !set_member(i, cone->Adom) && set_subset(common,cone->Aicd[i-1])){
      adj=FALSE;
    }
  }
_L99:;
  return adj;
} 

void dd_WriteInputIncidence(FILE *f, dd_PolyhedraPtr poly)
{
  dd_SetFamilyPtr I;

  if (poly->child->AicdGenerated==FALSE) ComputeAicd(poly->child);
  switch (poly->representation) {
  case Inequality:
    fprintf(f,"icd_file: Incidence of inequalities and generators\n");
    break;

  case Generator:
   fprintf(f,"ecd_file: Incidence of generators and inequalities\n");
    break;

  default:
    break;
  }

  I=dd_CopyInputIncidence(poly);
  dd_WriteSetFamilyCompressed(f,I);
  dd_FreeSetFamily(I);
}

void dd_WriteInputAdjacency(FILE *f, dd_PolyhedraPtr poly)
{
  dd_SetFamilyPtr A;

  if (poly->child->AicdGenerated==FALSE){
    ComputeAicd(poly->child);
  }
  switch (poly->representation) {
  case Inequality:
    fprintf(f, "iad_file: Adjacency of inequalities\n");
    break;

  case Generator:
    fprintf(f, "ead_file: Adjacency of generators\n");
    break;

  default:
    break;
  }
  A=dd_CopyInputAdjacency(poly);
  dd_WriteSetFamilyCompressed(f,A);
  dd_FreeSetFamily(A);
}


void dd_WriteProgramDescription(FILE *f)
{
  fprintf(f, "* cdd: a double description code:%s\n", DDVERSION);
  fprintf(f, "* compiled for %s arithmetic.\n", ARITHMETIC);
  fprintf(f,"* %s\n",COPYRIGHT);
}

void dd_WriteTimes(FILE *f, time_t starttime, time_t endtime)
{
  long ptime,ptime_sec,ptime_minu, ptime_hour;
  
/* 
   ptime=difftime(endtime,starttime); 
   This function is ANSI standard, but not available sometime 
*/
  ptime=(endtime - starttime);     
 /* This is to replace the line above, but it may not give 
    correct time in seconds */ 
  ptime_hour=ptime/3600;
  ptime_minu=(ptime-ptime_hour*3600)/60;
  ptime_sec=ptime%60;
  fprintf(f,"*Computation starts     at %s",asctime(localtime(&starttime)));
  fprintf(f,"*            terminates at %s",asctime(localtime(&endtime)));
  fprintf(f,"*Total processor time = %ld seconds\n",ptime);
  fprintf(f,"*                     = %ld h %ld m %ld s\n",ptime_hour, ptime_minu, ptime_sec);
}

void dd_WriteDDTimes(FILE *f, dd_PolyhedraPtr poly)
{
  dd_WriteTimes(f,poly->child->starttime,poly->child->endtime);
}

void dd_WriteLPTimes(FILE *f, dd_LPPtr lp)
{
  dd_WriteTimes(f,lp->starttime,lp->endtime);
}

void dd_WriteRunningMode(FILE *f, dd_ConePtr cone)
{

  switch (cone->HalfspaceOrder) {

    case MinIndex:
      fprintf(f, "minindex\n");
      break;

    case MaxIndex:
      fprintf(f, "maxindex\n");
      break;

    case MinCutoff:
      fprintf(f, "mincutoff\n");
      break;

    case MaxCutoff:
      fprintf(f, "maxcutoff\n");
    break;

    case MixCutoff:
      fprintf(f, "mixcutoff\n");
      break;

    case LexMin:
      fprintf(f, "lexmin\n");
      break;

    case LexMax:
      fprintf(f, "lexmax\n");
      break;

    case RandomRow:
      fprintf(f, "random  %d\n",cone->rseed);
      break;

    default: break;
  }
}


void WriteCompletionStatus(FILE *f, dd_ConePtr cone)
{
  if (cone->Iteration<cone->m && cone->CompStatus==AllFound) {
    fprintf(f,"*Computation completed at Iteration %4ld.\n", cone->Iteration);
  } 
  if (cone->CompStatus == RegionEmpty) {
    fprintf(f,"*Computation completed at Iteration %4ld because the region found empty.\n",cone->Iteration);
  }   
}

void dd_WritePolyFile(FILE *f, dd_PolyhedraPtr poly)
{
  dd_WriteAmatrix(f,poly->A,poly->m,poly->d);
}


void dd_WriteErrorMessages(FILE *f, dd_ErrorType Error)
{
  switch (Error) {

  case DimensionTooLarge:
    fprintf(f, "*Input Error: Input matrix is too large:\n");
    fprintf(f, "*Please increase MMAX and/or NMAX in the source code and recompile.\n");
    break;

  case IFileNotFound:
    fprintf(f, "*Input Error: Specified input file does not exist.\n");
    break;

  case OFileNotOpen:
    fprintf(f, "*Output Error: Specified output file cannot be opened.\n");
    break;

  case ImproperInputFormat:
    fprintf(f,"*Input Error: Input format is not correct.\n");
    fprintf(f,"*Format:\n");
    fprintf(f," begin\n");
    fprintf(f,"   m   n  NumberType(real, rational or integer)\n");
    fprintf(f,"   b  -A\n");
    fprintf(f," end\n");
    break;

  case EmptyVrepresentation:
    fprintf(f, "*Input Error: V-representation is empty:\n");
    fprintf(f, "*cddlib does not accept this trivial case for which output can be any inconsistent system.\n");
    break;

  case NoLPObjective:
    fprintf(f, "*LP Error: No LP objective (max or min) is set.\n");
    break;

  case NoRealNumberSupport:
    fprintf(f, "*LP Error: The binary (with GMP Rational) does not support Real number input.\n");
    fprintf(f, "         : Use a binary compiled with -DCDOUBLE option.\n");
    break;

  case NoError:
    fprintf(f,"*No Error found.\n");
    break;
  }
}

dd_SetFamilyPtr dd_CopyIncidence(dd_PolyhedraPtr poly)
{
  dd_RayPtr RayPtr;
  dd_SetFamilyPtr F=NULL;
  dd_bigrange i=0;

  if (poly->child==NULL || poly->child->CompStatus!=AllFound) goto _L99;
  F=dd_CreateSetFamily(poly->child->FeasibleRayCount, poly->child->m);
  RayPtr = poly->child->FirstRay;
  while (RayPtr != NULL) {
    if (RayPtr->feasible) {
      set_copy(F->set[i], RayPtr->ZeroSet);
      i++;  /* 086 */
    }
    RayPtr = RayPtr->Next;
  }
_L99:;
  return F;
}

dd_SetFamilyPtr dd_CopyInputIncidence(dd_PolyhedraPtr poly)
{
  dd_rowrange i;
  dd_SetFamilyPtr F=NULL;

  if (poly->child==NULL || poly->child->CompStatus!=AllFound) goto _L99;
  if (poly->child->AicdGenerated==FALSE) ComputeAicd(poly->child);
  F=dd_CreateSetFamily(poly->child->m, poly->child->FeasibleRayCount);
  for(i=0; i< poly->child->m; i++){
    set_copy(F->set[i], poly->child->Aicd[i]);
  }
_L99:;
  return F;
}

dd_SetFamilyPtr dd_CopyAdjacency(dd_PolyhedraPtr poly)
{
  dd_RayPtr RayPtr1,RayPtr2;
  dd_SetFamilyPtr F=NULL;
  long pos1, pos2;
  boolean adj;

  if (poly->child==NULL || poly->child->CompStatus!=AllFound) goto _L99;
  F=dd_CreateSetFamily(poly->child->FeasibleRayCount, 
      poly->child->FeasibleRayCount);
  poly->child->LastRay->Next=NULL;
  for (RayPtr1=poly->child->FirstRay, pos1=1;RayPtr1 != NULL; 
				RayPtr1 = RayPtr1->Next, pos1++){
    for (RayPtr2=poly->child->FirstRay, pos2=1; RayPtr2 != NULL; 
					RayPtr2 = RayPtr2->Next, pos2++){
      if (RayPtr1!=RayPtr2){
        CheckAdjacency(poly->child, &RayPtr1, &RayPtr2, &adj);
        if (adj){
          set_addelem(F->set[pos1-1], pos2);
        }
      }
    }
  }
_L99:;
  return F;
}

dd_SetFamilyPtr dd_CopyInputAdjacency(dd_PolyhedraPtr poly)
{
  dd_rowrange i,j;
  dd_SetFamilyPtr F=NULL;

  if (poly->child==NULL || poly->child->CompStatus!=AllFound) goto _L99;
  if (poly->child->AicdGenerated==FALSE) ComputeAicd(poly->child);
  F=dd_CreateSetFamily(poly->child->m, poly->child->m);
  for (i=1; i<=poly->child->m; i++){
    for (j=1; j<=poly->child->m; j++){
      if (i!=j && InputAdjacentQ(poly->child, i, j)) {
        set_addelem(F->set[i-1],j);
      }
    }
  }
_L99:;
  return F;
}

dd_MatrixPtr dd_CopyGenerators(dd_PolyhedraPtr poly)
{
  dd_RayPtr RayPtr;
  dd_MatrixPtr M=NULL;
  dd_rowrange i=0,total;
  dd_colrange j,j1;
  mytype b;

  dd_init(b);
  total=poly->child->LinearityDim + poly->child->FeasibleRayCount;
  if (poly->child->newcol[1]==0) total=total-1;
  if (poly->child==NULL || poly->child->CompStatus!=AllFound) goto _L99;
  if (poly->representation==Inequality){
    M=dd_CreateMatrix(total, poly->d);

    RayPtr = poly->child->FirstRay;
    while (RayPtr != NULL) {
      if (RayPtr->feasible) {
        CopyRay(M->matrix[i], poly->d, RayPtr, Generator, poly->child->newcol);
        i++;  /* 086 */
      }
      RayPtr = RayPtr->Next;
    }
    for (j=2; j<=poly->d; j++){
      if (poly->child->newcol[j]==0){
        dd_set(b,poly->child->Bsave[0][j-1]);
        if (dd_Nonzero(b)){
          dd_set(M->matrix[i][0],dd_one);  /* Normalize */
          for (j1=1; j1<poly->d; j1++) 
            dd_div(M->matrix[i][j1],(poly->child->Bsave[j1][j-1]),b);
        } else {
          for (j1=0; j1<poly->d; j1++)
            dd_set(M->matrix[i][j1],poly->child->Bsave[j1][j-1]);
        }
        set_addelem(M->linset, i+1);
        i++;
      }     
    }
  } else {
    M=dd_CreateMatrix(poly->m, poly->d);
    CopyAmatrix(M->matrix, poly->A, poly->m, poly->d);
    for (i=1; i<=poly->m; i++) 
      if (poly->EqualityIndex[i]==1) set_addelem(M->linset,i);
  }
  MatrixIntegerFilter(M);
  M->representation=Generator;
_L99:;
  dd_clear(b);
  return M;
}

dd_MatrixPtr dd_CopyInequalities(dd_PolyhedraPtr poly)
{
  dd_RayPtr RayPtr;
  dd_MatrixPtr M=NULL;
  dd_rowrange i=0,total;
  dd_colrange j,j1;

  total=poly->child->LinearityDim + poly->child->FeasibleRayCount;
  if (poly->child->newcol[1]==0) total=total-1;
  if (poly->child==NULL || poly->child->CompStatus!=AllFound) goto _L99;
  if (poly->representation==Generator){
    M=dd_CreateMatrix(total, poly->d);

    RayPtr = poly->child->FirstRay;
    while (RayPtr != NULL) {
      if (RayPtr->feasible) {
        CopyRay(M->matrix[i], poly->d, RayPtr, Inequality, poly->child->newcol);
        i++;  /* 086 */
      }
      RayPtr = RayPtr->Next;
    }
    for (j=2; j<=poly->d; j++){
      if (poly->child->newcol[j]==0){
        for (j1=0; j1<poly->d; j1++) 
          dd_set(M->matrix[i][j1],poly->child->Bsave[j1][j-1]);
        set_addelem(M->linset, i+1);
        i++;
      }     
    }
  } else {
    M=dd_CreateMatrix(poly->m, poly->d);
    CopyAmatrix(M->matrix, poly->A, poly->m, poly->d);
    for (i=1; i<=poly->m; i++) 
      if (poly->EqualityIndex[i]==1) set_addelem(M->linset,i);
  }
  MatrixIntegerFilter(M);
  M->representation=Inequality;
_L99:;
  return M;
}

/****************************************************************************************/
/*  rational number (a/b) read is taken from Vinci by Benno Bueeler and Andreas Enge    */
/****************************************************************************************/
void sread_rational_value (char *s, mytype value)
   /* reads a rational value from the specified string "s" and assigns it to "value"    */
   
{
   char     *numerator_s=0, *denominator_s=0, *position;
   int      sign = 1;
   double   numerator, denominator;
#if defined GMPRATIONAL
   long num;
   unsigned long den;
#else
   double rvalue;
#endif
  
   /* determine the sign of the number */
   numerator_s = s;
   if (s [0] == '-')
   {  sign = -1;
      numerator_s++;
   }
   else if (s [0] == '+')
      numerator_s++;
      
   /* look for a sign '/' and eventually split the number in numerator and denominator */
   position = strchr (numerator_s, '/');
   if (position != NULL)
   {  *position = '\0'; /* terminates the numerator */
      denominator_s = position + 1;
   };

   /* determine the floating point values of numerator and denominator */
   numerator=atol (numerator_s);
  
   if (position != NULL)
   { 
     denominator=atol (denominator_s);  
   }
   else denominator = 1;

/* 
   printf("\nrational_read: numerator %f\n",numerator);
   printf("rational_read: denominator %f\n",denominator);
   printf("rational_read: sign %d\n",sign); 
*/

#if defined GMPRATIONAL
   num=(long)sign * numerator;
   den=(unsigned long) denominator;
   mpq_set_si(value, num, den);
#elif defined GMPFLOAT
   rvalue=sign * numerator/ (signed long) denominator;
   mpf_set_d(value, rvalue);
#else
   rvalue=sign * numerator/ (signed long) denominator;
   ddd_set_d(value, rvalue);
#endif
   if (debug) {
     printf("rational_read: "); dd_WriteNumber(stdout,value); printf("\n");
   }
}
   

void fread_rational_value (FILE *f, mytype value)
   /* reads a rational value from the specified file "f" and assigns it to "value"      */
   
{
   char     number_s [255];
   mytype rational_value;
   
   dd_init(rational_value);
   fscanf (f, "%s ", number_s);
   sread_rational_value (number_s, rational_value);
   dd_set(value,rational_value);
   dd_clear(rational_value);
}
   
/****************************************************************************************/


/* end of cddio.c */

