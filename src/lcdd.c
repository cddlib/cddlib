/* lcdd.c: Main test program to call the cdd library cddlib
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   and  David Avis, avis@mutt.cs.mcgill.ca
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

/*  This program behaves like lrs by David Avis.
    Usage: lcdd filein fileout
           lcdd filein          output to stdout
           lcdd                 input stdin, output stdout

    This allows things like
    lcdd file | lcdd    (should give a minimal rep of the input file on stdout)
    lcdd file | lrs
    lcdd < filein
*/

#include "setoper.h"
#include "cdd.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>


int main(int argc, char *argv[])
{
  dd_ErrorType err;

  dd_set_global_constants();  /* First, this must be called. */
  dd_log=dd_TRUE; /* Output log */

  if (argc > 2)  
    dd_DDFile2File(argv[1],argv[2],&err);

  else if (argc > 1)
    dd_DDFile2File(argv[1],"**stdout",&err);
  else
    dd_DDFile2File("**stdin","**stdout",&err);
  return 0;
}

/* end of lcdd.c */
