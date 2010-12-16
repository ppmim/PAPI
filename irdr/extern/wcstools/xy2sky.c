/* File xy2sky.c
 * February 23, 2006
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 * Send bug reports to dmink@cfa.harvard.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>
#include "libwcs/wcs.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcscat.h"

static void PrintUsage();
extern struct WorldCoor *GetWCSFITS();	/* Read WCS from FITS or IRAF header */
static void PrintHead();

static int verbose = 0;		/* verbose/debugging flag */
static int append = 0;		/* append input line flag */
static int tabtable = 0;	/* tab table output flag */
static int identifier = 0;	/* 1st column=id flag */
static char coorsys[16];
static int linmode = -1;
static int face = 1;
static int ncm = 0;
static int printhead = 0;
static char printonly = 'n';
static int version = 0;		/* If 1, print only program name and version */

const char *RevMsg = "XY2SKY WCSTools 3.6.4, 3 May 2006, Doug Mink (dmink@cfa.harvard.edu)";

int
main (ac, av)
int ac;
char **av;
{
    char *str;
    char wcstring[64];
    char lstr = 64;
    int ndecset = 0;
    int degout = 0;
    int ndec = 3;
    int nlines, iline, entx, enty, entmag;
    int i;
    double x, y, mag;
    FILE *fd;
    char *ln, *listname;
    char linebuff[1024];
    char *line;
    char *fn;
    struct WorldCoor *wcs;
    char xstr[32], ystr[32], mstr[32];
    char keyword[16];
    char temp[64];
    int ncx = 0;
    char *cofile;
    char *cobuff;
    char *ctab;
    char *space, *cstr, *dstr;
    int nterms;
    double mag0, magx;
    double coeff[5];
    struct Tokens tokens;
    struct TabTable *tabxy;

    *coorsys = 0;
    cofile = NULL;

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	PrintUsage (str);
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	PrintUsage (str);
	}

    /* crack arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
	char c;
	while (c = *++str)
    	switch (c) {

    	case 'v':	/* more verbosity */
    	    verbose++;
    	    break;

	case 'a':	/* Append input line to sky position */
	    append++;
    	    break;

	case 'b':	/* Output B1950 coordinates */
	    strcpy (coorsys,"B1950");
    	    break;

	case 'c':	/* Magnitude conversion coefficients */
	    if (ac < 2)
		PrintUsage (str);
	    cofile = *++av;
	    ac--;
    	    break;

	case 'd':	/* Output degrees instead of hh:mm:ss dd:mm:ss */
            degout++;
            break;

	case 'e':	/* Output galactic coordinates */
	    strcpy (coorsys,"ecliptic");
            degout++;
    	    break;

	case 'f':	/* Face to use for projections with 3rd dimension */
	    if (ac < 2)
		PrintUsage (str);
	    face = atoi (*++av);
	    (void) wcszin (face);
	    ac--;
	    break;

	case 'g':	/* Output galactic coordinates */
	    strcpy (coorsys,"galactic");
            degout++;
    	    break;

	case 'h':	/* Print column heading */
	    printhead++;
	    break;

	case 'i':	/* 1st column = star id */
	    identifier++;
	    break;

	case 'j':	/* Output J2000 coordinates */
	    strcpy (coorsys,"J2000");
    	    break;

	case 'l':	/* Mode for output of linear coordinates */
	    if (ac < 2)
		PrintUsage (str);
	    linmode = atoi (*++av);
	    ac--;
	    break;

	case 'm':	/* column for magnitude */
	    ncm = atoi (*++av);
	    ac--;
    	    break;

	case 'n':	/* Number of decimal places in output sec or deg */
	    if (ac < 2)
		PrintUsage (str);
	    ndec = atoi (*++av);
	    ndecset++;
	    ac--;
	    break;

	case 'o':   /* Output only the following part of the coordinates */
	    if (ac < 2)
		PrintUsage (str);
	    av++;
	    printonly = **av;
	    ac--;
	    break;


	case 'q':	/* Equinox for output */
	    if (ac < 2)
		PrintUsage (str);
	    strcpy (coorsys, *++av);
	    ac--;
	    break;

	case 't':	/* Output tab table */
	    tabtable++;
    	    break;

	case 'x':	/* column for X; Y is next column */
	    ncx = atoi (*++av);
	    ac--;
    	    break;

    	case 'z':	/* Use AIPS classic WCS */
    	    setdefwcs (WCS_ALT);
    	    break;

    	default:
    	    PrintUsage (str);
    	    break;
    	}
    }

    /* There are ac remaining file names starting at av[0] */
    if (ac == 0)
	PrintUsage ();

    fn = *av++;
    if (isfits(fn) || isiraf(fn) || istiff(fn) || isjpeg(fn) || isgif(fn)) {
	wcs = GetWCSFITS (fn, verbose);
	if (nowcs (wcs)) {
	    printf ("%s: No WCS for file, cannot compute image size\n", fn);
	    wcsfree (wcs);
	    exit(1);
	    }
	if (linmode > -1)
	    setwcslin (wcs, linmode);
	if (wcs->sysout == WCS_B1950)
	    wcsoutinit (wcs, "B1950");
	if (wcs->sysout == WCS_J2000)
	    wcsoutinit (wcs, "J2000");
	if (*coorsys)
	    wcsoutinit (wcs, coorsys);

	/* Get magnitude conversion polynomial coefficients */
	if (cofile != NULL) {
	    cobuff = getfilebuff (cofile);
	    mag0 = 0.0;
	    (void) agetr8 (cobuff, "mag0", &mag0);
	    for (i = 0; i < 5; i++) {
		coeff[i] = 0.0;
		sprintf (keyword, "mcoeff%d", i);
		(void) agetr8 (cobuff, keyword, &coeff[i]);
		if (coeff[i] != 0.0) nterms = i + 1;
		}
	    if (ncm == 0)
		ncm = 3;
	    }
	if (degout) {
	    wcs->degout = degout;
            if (!ndecset) {
        	ndec = 5;
		ndecset++;
		}
	    }
	if (tabtable)
	    wcs->tabsys = 1;
	if (ndecset)
	    wcs->ndec = ndec;
	}

    else {
	fprintf (stderr, "*** XY2SKY error: Image filename must be first argument ***\n");
	exit (1);
	}

    /* Loop through arguments */
    while (ac-- > 1) {

	/* Process file of image coordinates */
	listname = *av;
	if (listname[0] == '@') {
	    ln = listname;
	    while (*ln++)
		*(ln-1) = *ln;
	    if (strcmp (listname,"STDIN")==0 || strcmp (listname,"stdin")==0) {
		fd = stdin;
		nlines = 10000;
		tabxy = NULL;
		if (printhead || verbose || tabtable)
		    PrintHead (fn, wcs, NULL, listname);
		}
	    else if (istab (listname)) {
		tabxy = tabopen (listname);
		nlines = tabxy->nlines;
		tabtable = 1;
		wcs->tabsys = 1;
		if (append) {
		    tabtable++;
		    PrintHead (fn, wcs, tabxy, listname);
		    }
		else
		    PrintHead (fn, wcs, NULL, listname);

		/* Find columns for X and Y */
		if (!(entx = tabcol (tabxy, "X")))
		    entx = tabcol (tabxy, "x");
		if (!(enty = tabcol (tabxy, "Y")))
		    enty = tabcol (tabxy, "y");

		/* Find column for magnitude */
		if (!(entmag = tabcol (tabxy, "MAG"))) {
		    if (!(entmag = tabcol (tabxy, "mag"))) {
			if (!(entmag = tabcol (tabxy, "magv"))) {
			    if (!(entmag = tabcol (tabxy, "magj")))
				entmag = tabcol (tabxy, "magr");
			    }
			}
		    }
		}
	    else {
		if (printhead || verbose || tabtable)
		    PrintHead (fn, wcs, NULL, listname);
		tabxy = NULL;
		nlines = getfilelines (listname);
		fd = fopen (listname, "r");
		if (fd == NULL) {
		    fprintf (stderr, "Cannot read file %s\n", listname);
		    nlines = 0;
		    }
		}
	    for (iline = 0; iline < nlines; iline++) {
		if (tabxy != NULL) {

		    /* Read line for next position */
		    if ((line = gettabline (tabxy, iline+1)) == NULL) {
			fprintf (stderr,"Cannot read star %d\n", iline);
			break;
			}

		    /* Extract x, y, and magnitude */
		    (void) setoken (&tokens, line, "tab");
		    x = tabgetr8 (&tokens, entx);
		    y = tabgetr8 (&tokens, enty);
		    if (ncm)
			mag = tabgetr8 (&tokens, entmag);
		    }
		else {
		    if (!fgets (linebuff, 1023, fd))
			break;
		    line = linebuff;
		    if (line[0] == '#')
			continue;
		    if (ncm || ncx)
			setoken (&tokens, line, "");
		    if (ncx) {
			getoken (&tokens, ncx, xstr, 31);
			getoken (&tokens, ncx+1, ystr, 31);
			}
		    else if (identifier)
			sscanf (line,"%s %s %s", temp, xstr, ystr);
		    else
			sscanf (line,"%s %s", xstr, ystr);
		    x = atof (xstr);
		    y = atof (ystr);
		    if (ncm) {
			wcs->printsys = 0;
			getoken (&tokens, ncm, mstr, 31);
			mag = atof (mstr);
			}
		    }
		if (pix2wcst (wcs, x, y, wcstring, lstr)) {
		    /* Remove coordinate system if tab table output
		    if (tabtable) {
			ctab = strrchr (wcstring, (char) 9);
			*ctab = (char) 0;
			} */
		    if (wcs->sysout == WCS_ECLIPTIC) {
			sprintf(temp,"%.5f",wcs->epoch);
			strcat (wcstring, " ");
			strcat (wcstring, temp);
			}
		    if (ncm) {
			magx = polcomp (mag, mag0, nterms, coeff);
			if (tabtable)
			    sprintf(temp,"	%6.2f", magx);
			else
			    sprintf(temp," %6.2f", magx);
			strcat (wcstring, temp);
			}
		    if (append) {
			if (tabtable)
			    printf ("%s	%s", wcstring, line);
			else
			    printf ("%s %s", wcstring, line);
			}
		    else {
			if (tabtable)
			    printf ("%s	", wcstring);
			else
			    printf ("%s ", wcstring);

			if (verbose && !tabtable)
			    printf (" <- ");
			if (wcs->nxpix > 9999 || wcs->nypix > 9999) {
			    if (tabtable)
				printf ("%9.3f	%9.3f	",x, y);
			    else
				printf ("%9.3f %9.3f ",x, y);
			    }
			else if (wcs->nxpix > 999 || wcs->nypix > 999) {
			    if (tabtable)
				printf ("%8.3f	%8.3f	",x, y);
			    else
				printf ("%8.3f %8.3f ",x, y);
			    }
			else {
			    if (tabtable)
				printf ("%7.3f	%7.3f	",x, y);
			    else
				printf ("%7.3f %9.3f ",x, y);
			    }
			if (wcs->naxis > 2) {
			    if (tabtable)
				printf ("%2d	", face);
			    else
				printf ("%2d  ", face);
			    }
			}
		    printf ("\n");
		    }
		}
	    av++;
	    }

	/* Process image coordinates from the command line */
	else if (ac > 1) {
	    if (printhead || verbose) {
		PrintHead (fn, wcs, NULL, NULL);
		printhead = 0;
		}
	    x = atof (*av);
	    ac--;
	    y = atof (*++av);
	    if (pix2wcst (wcs, x, y, wcstring, lstr)) {
		if (wcs->sysout == WCS_ECLIPTIC) {
		    sprintf(temp,"%.5f",wcs->epoch);
		    strcat (wcstring, " ");
		    strcat (wcstring, temp);
		    }
		if (ncm) {
		    magx = polcomp (mag, mag0, nterms, coeff);
		    if (tabtable)
			sprintf(temp,"	%6.2f", magx);
		    else
			sprintf(temp," %6.2f", magx);
		    strcat (wcstring, temp);
		    }
		if (printonly == 'r') {
		    space = strchr (wcstring, ' ');
		    *space = (char) 0;
		    printf ("%s ", wcstring);
		    }
		else if (printonly == 'd') {
		    dstr = strchr (wcstring, ' ') + 1;
		    space = strchr (dstr, ' ');
		    *space = (char) 0;
		    printf ("%s ", dstr);
		    }
		else if (printonly == 's') {
		    dstr = strchr (wcstring, ' ') + 1;
		    cstr = strchr (dstr, ' ') + 1;
		    printf ("%s ", cstr);
		    }
		else {
		    printf ("%s ", wcstring);
		    if (verbose && !tabtable)
			printf (" <- ");
		    if (wcs->nxpix > 9999 || wcs->nypix > 9999) {
			if (tabtable)
			    printf ("%9.3f	%9.3f	",x, y);
			else
			    printf ("%9.3f %9.3f ",x, y);
			}
		    else if (wcs->nxpix > 999 || wcs->nypix > 999) {
			if (tabtable)
			    printf ("%8.3f	%8.3f	",x, y);
			else
			    printf ("%8.3f %8.3f ",x, y);
			}
		    else {
			if (tabtable)
			    printf ("%7.3f	%7.3f	",x, y);
			else
			    printf ("%7.3f %9.3f ",x, y);
			}
		    if (wcs->naxis > 2) {
			if (tabtable)
			    printf ("%2d	", face);
			else
			    printf ("%2d  ", face);
			}
		    }
		printf ("\n");
		}
	    av++;
	    }
	}

    wcsfree (wcs);
    return (0);
}

static void
PrintUsage (command)
char	*command;

{
    fprintf (stderr,"%s\n",RevMsg);
    if (version)
	exit (-1);
    if (command != NULL) {
	if (command[0] == '*')
	    fprintf (stderr, "%s\n", command);
	else
	    fprintf (stderr, "* Missing argument for command: %c\n", command[0]);
	exit (1);
	}
    fprintf (stderr,"Compute RA Dec from X Y using WCS in FITS and IRAF image files\n");
    fprintf (stderr,"Usage: [-abdjgv] [-n ndec] file.fits x1 y1 ... xn yn\n");
    fprintf (stderr,"  or : [-abdjgv] [-n ndec] file.fits @listfile\n");
    fprintf (stderr,"  -a: append input line after output position\n");
    fprintf (stderr,"  -b: B1950 (FK4) output\n");
    fprintf (stderr,"  -c: file with coefficients for magnitude conversion\n");
    fprintf (stderr,"  -d: RA and Dec output in degrees\n");
    fprintf (stderr,"  -e: ecliptic longitude and latitude output\n");
    fprintf (stderr,"  -f: Number of face to use for 3-d projection\n");
    fprintf (stderr,"  -g: galactic longitude and latitude output\n");
    fprintf (stderr,"  -i: first column is star id; 2nd, 3rd are x,y position\n");
    fprintf (stderr,"  -j: J2000 (FK5) output\n");
    fprintf (stderr,"  -l: mode for output of LINEAR WCS coordinates\n");
    fprintf (stderr,"  -m: column for magnitude (defaults to 3 if -c used)\n");
    fprintf (stderr,"  -n: number of decimal places in output RA seconds\n");
    fprintf (stderr,"  -o r|d|s: print only ra, dec, or coordinate system\n");
    fprintf (stderr,"  -q: output equinox if not 2000 or 1950\n");
    fprintf (stderr,"  -t: tab table output\n");
    fprintf (stderr,"  -v: verbose\n");
    fprintf (stderr,"  -x: Column for image X coordinate; Y follows\n");
    fprintf (stderr,"  -z: use AIPS classic projections instead of WCSLIB\n");
    exit (1);
}

static void
PrintHead (fn, wcs, tabxy, listfile)

char	*fn;		/* Name of file containing list of x,y coordinates */
struct WorldCoor *wcs;	/* World coordinate system structure */
struct TabTable *tabxy;	/* tab table structure */
char *listfile;		/* Name of file with list of input coordinates */

{
    char newline = (char) 10;
    char *eol;
    char *ctab, *tabhead;
    int lhead, i;

    if (tabtable) {
	printf ("image	%s\n", fn);
	if (listfile != NULL)
	    printf ("listfile	%s\n", listfile);
	printf ("radecsys	%s\n",wcs->radecout);
	printf ("epoch	%.4f\n",wcs->epoch);
	printf ("program	%s\n",RevMsg);
	if (wcs->sysout == WCS_B1950 || wcs->sysout == WCS_J2000)
	    printf ("ra         	dec         	");
	else if (wcs->sysout == WCS_GALACTIC)
	    printf ("glon     	glat     	");
	else if (wcs->sysout == WCS_ECLIPTIC)
	    printf ("elon     	elat     	");
	if (ncm)
	    printf ("mag   	");
	if (tabxy != NULL) {
	    eol = strchr (tabxy->tabhead, newline);
	    *eol = (char) 0;
	    if (wcs->sysout == WCS_B1950 || wcs->sysout == WCS_J2000) {
		lhead = strlen (tabxy->tabhead);
		tabhead = (char *) calloc (lhead+4, sizeof (char));
		ctab = tabhead;
		for (i = 0; i < lhead; i++) {
		    *ctab++ = tabxy->tabhead[i];
		    if (*(ctab-2) == 'r' && *(ctab-1) == 'a')
			*ctab++ = '0';
		    if (*(ctab-3) == 'd' && *(ctab-2) == 'e' && *(ctab-1) == 'c')
			*ctab++ = '0';
		    }
		printf ("%s", tabhead);
		}
	    else
		printf ("%s", tabxy->tabhead);
	    *eol = newline;
	    }
	else {
	    printf ("x       	y       	");
	    if (wcs->naxis > 2)
		printf ("z    	");
	    }
	printf ("\n");
	if (wcs->degout)
	    printf ("---------	---------	");
	else
	    printf ("------------	------------	");
	if (ncm)
	    printf ("------	");
	if (tabxy != NULL) {
	    eol = strchr (tabxy->tabdash, newline);
	    *eol = (char) 0;
	    printf ("%s", tabxy->tabdash);
	    *eol = newline;
	    }
	else {
	    printf ("--------	--------");
	    if (wcs->naxis > 2)
		printf ("	-----");
	    }
	printf ("\n");
	}
    else {
        fprintf (stderr,"%s\n",RevMsg);
	if (listfile == NULL)
	    fprintf (stderr,
		"Print sky coordinates from %s image coordinates\n", fn);
	else
	    fprintf (stderr,
		"Print sky coordinates from %s image coordinates in %s\n",
		fn, listfile);
	if (wcs->sysout == WCS_ECLIPTIC || wcs->sysout == WCS_GALACTIC)
	    printf ("Longitude  Latitude   Sys    ");
	else
	    printf ("    RA           Dec       Sys  ");
	if (wcs->sysout == WCS_ECLIPTIC)
	    printf("  Epoch    ",wcs->epoch);
	if (verbose) printf ("    ");
	printf ("    X        Y\n");
	}
    return;
}

/*
 * Feb 23 1996	New program
 * Apr 24 1996	Version 1.1: Add B1950, J2000, or galactic coordinate output options
 * Jun 10 1996	Change name of subroutine which reads WCS
 * Aug 28 1996	Remove unused variables after lint
 * Nov  1 1996	Add options to set number of decimal places and output degrees
 *
 * Dec 15 1997	Print message if no WCS; read IRAF 2.11 header format
 * Dec 15 1997	Drop -> if output sky coordinates are in degrees
 * Dec 31 1997	Allow entire input line to be appended to sky position
 *
 * Jan  7 1998	Apply WFPC and WFPC2 pixel corrections if requested
 * Jan  7 1998	Add tab table output using -t
 * Jan 26 1998	Implement Mark Calabretta's WCSLIB
 * Jan 29 1998	Add -z for AIPS classic WCS projections
 * Feb 18 1998	Version 2.0: Full Calabretta implementation
 * Mar 12 1998	Version 2.1: IRAF TNX projection added
 * Mar 27 1998	Version 2.2: Polynomial plate fit added
 * Apr 24 1998	Increase size of WCS string from 40 to 64
 * Apr 28 1998	Change coordinate system flag to WCS_*
 * Apr 28 1998	Add output mode for linear coordinates
 * Apr 28 1998	Add ecliptic coordinate system output
 * May 13 1998	Allow arbitrary equinox for output coordinates
 * Jun 25 1998	Set WCS subroutine choice with SETDEFWCS()
 * Jul  7 1998	Change setlinmode() to setwcslin()
 * Jul  7 1998	Add -f for face to use in 3-d projection
 * Jul  7 1998	Add 3rd dimension in output
 * Nov 30 1998	Add version and help commands for consistency
 *
 * Mar 29 1999	Add -i option for X,Y after ID in input file (J.-B. Marquette)
 * Apr 29 1999	Drop pix2wcst declaration; it is in wcs.h
 * Jun  2 1999	Make tab output completely-defined catalog
 * Oct 15 1999	Free wcs using wcsfree()
 * Oct 22 1999	Link includes to libwcs
 *
 * Jan 28 2000	Call setdefwcs() with WCS_ALT instead of 1
 *
 * Jul 19 2001	Add -x to specify column for X coordinate; Y follows immediately
 * Jul 19 2001	Add -m to specify column for magnitude
 * Jul 19 2001	Add -c to specify column for magnitude coefficients
 * Jul 23 2001	Add code to calibrate magnitudes using polynomial from immatch
 * Jul 25 2001	Ignore lines with # in first column
 * Sep 12 2001	Fix output to match column headings
 * Oct 16 2001	Increase maximum input line length from 200 to 1024
 * Oct 16 2001	Add option to prepend coordinates to tab-separated table
 * Oct 19 2001	Change names of old ra and dec columns if prepending ra and dec
 * Dec 13 2001	Add -h for headings and add to verbose output
 *
 * Apr  8 2002	Free wcs structure if no WCS is found in file header
 *
 * Jun 19 2002	Add verbose argument to GetWCSFITS()
 *
 * Jan  7 2003	Fix bug which dropped declination for tab output from file
 * Jan  7 2003	Fix bug which failed to ignore #-commented-out input file lines
 * Apr  7 2003	Add -o to output only RA, Dec, or system
 * Oct 14 2003	Change naxes to naxes in wcs structure
 *
 * Feb 23 2006	Allow appended headers in TIFF, JPEG, and GIF files
 */
