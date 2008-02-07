/* Header file for setoper.c  */

/* setoper.c: 
 * A set operation library 
 * created by Komei Fukuda, Nov.14, 1993
 * last modified on June 1, 2000
 */

#ifndef  __SETOPER_H
#define  __SETOPER_H
#endif  /* __SETOPER_H */

#include <stdio.h>
#include <stdlib.h>

typedef unsigned long *set_type;   /* set type definition */

typedef unsigned char set_card_lut_t;

#if defined(__cplusplus)
extern "C" {
#endif

unsigned long set_blocks(long len);
void set_initialize(set_type *setp,long len);
void set_free(set_type set);
void set_emptyset(set_type set);
void set_copy(set_type setcopy,set_type set);
void set_addelem(set_type set, long elem);
void set_delelem(set_type set, long elem);
void set_int(set_type set,set_type set1,set_type set2);
void set_uni(set_type set,set_type set1,set_type set2);
void set_diff(set_type set,set_type set1,set_type set2);
void set_compl(set_type set,set_type set1);
int set_subset(set_type set1,set_type set2);
int set_member(long elem, set_type set);
long set_card(set_type set);
long set_groundsize(set_type set); /* output the size of the ground set */
void set_write(set_type set);
void set_fwrite(FILE *f,set_type set);
void set_fwrite_compl(FILE *f,set_type set); /* write the complement */
void set_binwrite(set_type set);
void set_fbinwrite(FILE *f,set_type set);

#if defined(__cplusplus)
}
#endif

/* End of File: setoper.h */

