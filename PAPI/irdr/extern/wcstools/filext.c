/* File filext.c
 * April 29, 2002
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

int
main (ac, av)
int ac;
char **av;
{
    char *fn;
    char *str;
    char *ext, *ext2;
    int i, lroot;

    /* crack arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
        char c;
        while (c = *++str)
        switch (c) {

        case 'v':       /* more verbosity */
            verbose++;
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
	ext = strrchr (fn, '[');
	if (ext != NULL) {
	    ext = ext + 1;
	    ext2 = strrchr (fn, ']');
	    if (ext2 != NULL)
		*ext2 = (char) 0;
	    }
	else {
	    ext = strrchr (fn, ',');
	    if (ext != NULL)
		ext = ext + 1;
	    else {
		ext = strrchr (fn, '.');
		if (ext != NULL)
		    ext = ext + 1;
		}
	    }
	printf ("%s\n", ext);
	}

    return (0);
}

static void
usage ()
{
    fprintf (stderr,"FILEXT: Print file name extension\n");
    fprintf(stderr,"Usage:  filext file1 file2 file3 ...\n");
    exit (1);
}
/* Apr 29 2002	New program
 */
