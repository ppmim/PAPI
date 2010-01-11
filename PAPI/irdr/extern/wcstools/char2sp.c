/* File char2sp.c
 * April 3, 2006
 * By Doug Mink Harvard-Smithsonian Center for Astrophysics)
 * Send bug reports to dmink@cfa.harvard.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>
#include "libwcs/fitsfile.h"

#define MAXKWD 50
static int maxnkwd = MAXKWD;

static void usage();
static void SetValues ();

static int verbose = 0;		/* verbose/debugging flag */
static int newimage0 = 0;
static int keyset = 0;
static int histset = 0;
static int krename = 0;
static char prefix[2];
static int version = 0;		/* If 1, print only program name and version */
static int first = 1;
static char spchar = '_';

static char *RevMsg = "CHAR2SP WCSTools 3.6.4, 3 May 2006, Doug Mink (dmink@cfa.harvard.edu)";

int
main (ac, av)
int ac;
char **av;
{
    char *str;
    char **kwd;
    int nkwd = 0;
    int lstr;
    int readlist = 0;
    int ifile;
    int ikwd, i, nc;

    nkwd = 0;
    kwd = (char **)calloc (maxnkwd, sizeof(char *));
    for (ikwd = 0; ikwd < maxnkwd; ikwd++) {
	kwd[ikwd] = NULL;
	}

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage();
	}

    /* crack arguments */
    for (av++; --ac > 0; av++) {
	if (*(str = *av) == '-') {
	    char c;
	    while (c = *++str)
	    switch (c) {

		case 's':	/* Replace this character with spaces in string arguments */
		    if (ac > 1) {
			spchar= *++av[0];
			ac--;
			}
		    break;

		case 'v':	/* more verbosity */
		    verbose++;
		    break;

		default:
		    usage();
		    break;
		}
	    }

	/* String component */
	else {
	    if (nkwd >= maxnkwd) {
		maxnkwd = maxnkwd * 2;
		kwd = (char **) realloc ((void *)kwd, maxnkwd);
		}
	    kwd[nkwd] = *av;
	    nkwd++;
	    }
	}

    if (nkwd <= 0 )
	usage ();
    else if (nkwd <= 0) {
	fprintf (stderr, "CHAR2SP: no string components specified\n");
	exit (1);
	}

    for (ikwd = 0; ikwd < nkwd; ikwd++) {
	lstr = strlen (kwd[ikwd]);
	for (i = 0; i < lstr; i++) {
	    if (kwd[ikwd][i] == spchar)
		kwd[ikwd][i] = ' ';
	    }
	printf ("%s", kwd[ikwd]);
	if (ikwd < nkwd - 1)
	    printf (" ");
	else
	    printf ("\n");
	}

    free (kwd);

    return (0);
}

static void
usage ()
{
    fprintf (stderr,"%s\n",RevMsg);
    if (version)
	exit (-1);
    fprintf (stderr,"Replaces a specified character in a string with spaces (def=_)\n");
    fprintf(stderr,"Usage: [-v][-s char] string [string2][...][string n]\n");
    fprintf(stderr,"  -s [char]: Replace this character with spaces in output\n");
    fprintf(stderr,"  -v: Verbose\n");
    exit (1);
}

/* Jan  8 2002	New program
 *
 * Apr  3 2006	Declare main to be int
 */
