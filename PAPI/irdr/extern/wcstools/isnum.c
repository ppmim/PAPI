/* File isnum.c
 * April 3, 2006
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 * Send bug reports to dmink@cfa.harvard.edu
 *
 * Return 1 if argument is an integer, 2 if it is floating point, else 0
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libwcs/fitshead.h"

const char *RevMsg = "ISNUM WCSTools 3.6.4, 3 May 2006, Doug Mink (dmink@cfa.harvard.edu)";

int
main (ac, av)
int ac;
char **av;
{
    char *fn;
    char *str;

    /* Check for version or help command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help")) {
	fprintf (stderr,"%s\n",RevMsg);
	fprintf (stderr,"Usage: Return 1 if argument is an integer, ");
	fprintf (stderr,"2 if it is floating point, else 0\n");
	exit (1);
	}
    else if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	fprintf (stderr,"%s\n",RevMsg);
	exit (1);
	}

    /* check to see if this is a number */
    else
	printf ("%d\n", isnum (str));

    exit (0);
}
/* Nov  7 2001	New program
 *
 * Apr 11 2005	Print version
 *
 * Apr  3 2006	Declare main to be int
 */
