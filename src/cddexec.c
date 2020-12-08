/* cddexec.c: executing cdd functions.  It mimics the interface of the program
   cdd_both_reps by Volker Braun <vbraun@stp.dias.ie> with much faster executions.

   The input is taken from stdin and can be either a
   H or V representation.

   Written by Komei Fukuda <fukuda@math.ethz.ch>, by using some part of cdd_both_reps.

*/

/*  This program is free software; you can redistribute it and/or modify
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

#include "setoper.h"
#include "cdd.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>


void compute_adjacency(dd_MatrixPtr Rep, dd_ErrorType* err_ptr)
{
  dd_SetFamilyPtr AdjacencyGraph;
  if (*err_ptr != dd_NoError) return;

  switch (Rep->representation) {
  case dd_Inequality:
    fprintf(stdout, "Facet adjacency\n");
    break;
  case dd_Generator:
    fprintf(stdout, "Vertex adjacency\n");
    break;
  case dd_Unspecified:
    fprintf(stderr, "unknown representation type!\n");
  default:
    fprintf(stderr, "This should be unreachable!\n");
    exit(2);
  }

  /* Output adjacency of vertices/rays/lines */
  if (Rep->rowsize > 0) {  /* workaround for bug with empty polyhedron */
    /* compute adjacent vertices/rays/lines */
    AdjacencyGraph = dd_Matrix2Adjacency(Rep, err_ptr);
    if (*err_ptr == dd_NoError) {
      dd_WriteSetFamily(stdout,AdjacencyGraph);
      dd_FreeSetFamily(AdjacencyGraph);
    }
  } else {
    printf("begin\n");
    printf("  0    0\n");
    printf("end\n");
  }

  printf("\n");
}


void compute_the_second_rep(dd_MatrixPtr M,
               dd_ErrorType* err_ptr)
{
  dd_PolyhedraPtr poly;
    dd_MatrixPtr A;
   /* compute the second representation */
  poly = dd_DDMatrix2Poly(M, err_ptr);
  if (*err_ptr!=dd_NoError) goto _L99;

  switch (poly->representation) {
    case dd_Inequality:
        fprintf(stdout, "The second representation:\n");
        A=dd_CopyGenerators(poly);
        dd_WriteMatrix(stdout,A);
        dd_FreeMatrix(A);
        break;

    case dd_Generator:
        fprintf(stdout, "The second representation:\n");
        A=dd_CopyInequalities(poly);
        dd_WriteMatrix(stdout,A);
        dd_FreeMatrix(A);
        break;

    default:
        break;
  }

  _L99:;
 }

void compute_the_second_repall(dd_MatrixPtr M,
                            dd_ErrorType* err_ptr)
{
    dd_PolyhedraPtr poly;
    dd_MatrixPtr A,B;
    /* compute the second representation */
    poly = dd_DDMatrix2Poly(M, err_ptr);
    if (*err_ptr!=dd_NoError) goto _L99;

    switch (poly->representation) {
        case dd_Inequality:
            A=dd_CopyGenerators(poly);
            B=dd_CopyInequalities(poly);

            fprintf(stdout, "The second representation:\n");
            dd_WriteMatrix(stdout,A);
            fprintf(stdout, "\nVertex incidence\n");
            dd_WriteIncidence(stdout,poly);
            fprintf(stdout, "\nVertex adjacency\n");
            dd_WriteAdjacency(stdout,poly);

            fprintf(stdout, "\nThe first (input) representation\n");
            dd_WriteMatrix(stdout,B);
            fprintf(stdout, "\nFacet incidence\n");
            dd_WriteInputIncidence(stdout,poly);
            fprintf(stdout, "\nFacet adjacency\n");
            dd_WriteInputAdjacency(stdout,poly);

            dd_FreeMatrix(A);
            dd_FreeMatrix(B);
            break;

        case dd_Generator:
            A=dd_CopyInequalities(poly);
            B=dd_CopyGenerators(poly);

            fprintf(stdout, "The second representation:\n");
            dd_WriteMatrix(stdout,A);
            fprintf(stdout, "\nFacet incidence\n");
            dd_WriteIncidence(stdout,poly);
            fprintf(stdout, "\nFacet adjacency\n");
            dd_WriteAdjacency(stdout,poly);

            fprintf(stdout, "\nThe first (input) representation\n");
            dd_WriteMatrix(stdout,B);
            fprintf(stdout, "\nVertex incidence\n");
            dd_WriteInputIncidence(stdout,poly);
            fprintf(stdout, "\nVertex adjacency\n");
            dd_WriteInputAdjacency(stdout,poly);

            dd_FreeMatrix(A);
            dd_FreeMatrix(B);
            break;

        default:
            break;
    }

_L99:;
}

void compute_redundancy(dd_MatrixPtr M,
                            dd_ErrorType* err_ptr)
{
    dd_rowrange i,m;
    dd_rowindex newpos;
    dd_rowset impl_linset,redset;

    m=M->rowsize;

    fprintf(stdout, "Canonicalize the matrix.\n");

    dd_MatrixCanonicalize(&M, &impl_linset, &redset, &newpos, err_ptr);

    if (*err_ptr!=dd_NoError) goto _L99;

    fprintf(stdout, "Implicit linearity rows are: "); set_fwrite(stdout, impl_linset);

    fprintf(stdout, "\nRedundant rows are: "); set_fwrite(stdout, redset);
    fprintf(stdout, "\n");

    fprintf(stdout, "Nonredundant representation:\n");
    fprintf(stdout, "The new row positions are as follows (orig:new).\nEach redundant row has the new number 0.\nEach deleted duplicated row has a number nagative of the row that\nrepresents its equivalence class.\n");

    for (i=1; i<=m; i++){
        fprintf(stdout, " %ld:%ld",i, newpos[i]);
    }
    fprintf(stdout, "\n");
    dd_WriteMatrix(stdout, M);

    set_free(redset);
    set_free(impl_linset);
    free(newpos);
    dd_FreeMatrix(M);

_L99:;
}


void usage(char *name)
{
  fprintf(stderr, "No known option specified, I don't know what to do!\n"
     "Usage:\n"
     "%s --option\n"
     "where --option is precisely one of the following:\n\n"
     "  --rep: Compute the second (H- or V-) representation.\n"
     "         The computed representation is minimal (without redundancy).\n"
     "\n"
     "  --repall: Compute the second (H- or V-) representation.\n"
     "    It outputs both the input and output representations,\n"
     "    as well as their incidence and adjacency relations.\n"
     "\n"
     "  --redcheck: Compute a minimal (non-redundant) representation.\n"
     "              This is sometimes called the redundancy removal.\n"
     "\n"
     "  --adjacency: Compute adjacency information only.\n"
     "    The input is assumed to be a minimal representation, as, for example, computed\n"
     "    by --redcheck. Warning, you will not get the correct answer if the input\n"
     "    representation is not minimal! The output is the vertex or facet graph,\n"
     "    depending on the input.\n"
     "\n"
     "The input data is a H- or V-representation in cdd's ine/ext format and\n"
     "is in each case read from stdin.\n",
     name);
}


enum command_line_arguments {REP, REPALL, ADJACENCY, REDCHECK};


int parse_arguments(char* arg, enum command_line_arguments* option)
{
  if (strcmp(arg,"--repall")==0) {
    *option = REPALL;
    return 0;
  }
  if (strcmp(arg,"--rep")==0) {
    *option = REP;
    return 0;
  }
  if (strcmp(arg,"--adjacency")==0) {
    *option = ADJACENCY;
    return 0;
  }
  if (strcmp(arg,"--redcheck")==0) {
    *option = REDCHECK;
    return 0;
  }

  fprintf(stderr, "Unknown option: %s\n", arg);
  return 1;
}


int main(int argc, char *argv[])
{
  dd_ErrorType err=dd_NoError;
  dd_MatrixPtr M=NULL;
  enum command_line_arguments option;

  if (argc!=2 || parse_arguments(argv[1],&option)) {
    usage(argv[0]);
    return 1;
  }

  dd_set_global_constants();

  /* Read data from stdin */
  M = dd_PolyFile2Matrix(stdin, &err);
  if (err != dd_NoError) {
    fprintf(stderr, "I was unable to parse the input data!\n");
    dd_WriteErrorMessages(stdout,err);
    dd_free_global_constants();
    return 1;
  }

  switch (option) {
  case REPALL:
    compute_the_second_repall(M,&err);
    break;
  case REP:
    compute_the_second_rep(M,&err);
    break;
  case ADJACENCY:
    compute_adjacency(M,&err);
    break;
  case REDCHECK:
    compute_redundancy(M,&err);
    break;
  default:
    fprintf(stderr, "unreachable option %d\n", option);
    exit(3); /* unreachable */
  }

  /* cleanup */
  if (option!=REDCHECK) dd_FreeMatrix(M); /* compute_redundancy modifies M and frees M. */
  if (err != dd_NoError) {
    dd_WriteErrorMessages(stdout,err);
  }

  dd_free_global_constants();
  return 0;
}
