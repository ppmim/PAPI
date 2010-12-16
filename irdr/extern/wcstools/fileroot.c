/* File fileroot.c
 * March 1, 2001
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 * Send bug reports to dmink@cfa.harvard.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>

static int verbose = 0;         /* verbose/debugging flag */
static int replace = 0;         /* character replacement flag */
static char c1, c2;
static void usage();
static void PrintRoot();

int
main (ac, av)
int ac;
char **av;
{
    char *fn;
    char *str;
    char *ext;
    int i, lroot;

    /* crack arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
        char c;
        while (c = *++str)
        switch (c) {

        case 'v':       /* more verbosity */
            verbose++;
            break;

        case 'r':       /* replace next character with one after it */
	    if (ac < 3) 
		usage();
	    av++;
	    c1 = *av[0];
	    ac--;
	    av++;
	    c2 = *av[0];
	    ac--;
            replace++;
            break;

        default:
            usage();
            break;
        }
    }

    /* There are ac remaining file names starting at av[0] */
    if (ac == 0)
        usage ();

    while (ac-- > 0) {
	fn = *av++;
	if (verbose)
    	    printf ("%s -> ", fn);
	ext = strrchr (fn, ',');
	if (ext != NULL)
	    *ext = 0;
	else {
	    ext = strrchr (fn, '.');
	    if (ext != NULL)
		*ext = 0;
	    }
	if (replace) {
	    lroot= strlen (fn);
	    for (i = 0; i < lroot; i++)
		if (fn[i] == c1) fn[i] = c2;
	    }
	printf ("%s\n", fn);
	}

    return (0);
}

static void
usage ()
{
    fprintf (stderr,"FILEROOT: Drop file name extension\n");
    fprintf(stderr,"Usage:  fileroot file1 file2 file3 ...\n");
    fprintf(stderr,"        fileroot -r c1 c2 file1 file2 file3 ...\n");
    fprintf(stderr,"        -r replaces c1 with c2 in file name\n");
    exit (1);
}
/* May  3 2000	New program
 * Sep 12 2000	Truncate at comma as well as period
 *
 * Mar  1 2000	Add character replacement
 */
