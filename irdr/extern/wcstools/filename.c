/* File fname.c
 * August 3, 1998
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 * Send bug reports to dmink@cfa.harvard.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>

static int verbose = 0;         /* verbose/debugging flag */
static void usage();

int
main (ac, av)
int ac;
char **av;
{
    char *fn;
    char *str;
    char *name;

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
	name = strrchr (fn, '/');
	if (name == NULL)
	    name = fn;
	else
	    name = name + 1;

	printf ("%s\n", name);
	}

    return (0);
}

static void
usage ()
{
    fprintf (stderr,"FILENAME: Drop directory from pathname\n");
    fprintf(stderr,"Usage:  filename path1 path2 path3 ...\n");
    exit (1);
}
