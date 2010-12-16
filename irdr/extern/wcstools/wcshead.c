/* File wcshead.c
 * July 19, 2004
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
#include "libwcs/wcs.h"

static void usage();
static void ListWCS();

static int verbose = 0;		/* verbose/debugging flag */
static int tabout = 0;		/* tab table output flag */
static int ndec = 3;		/* maximum number of decimal places for output*/
static int nchar = 24;		/* maximum number of characters for filename */
static int hms = 0;		/* 1 for output in hh:mm:ss dd:mm:ss */
static int nf = 0;		/* Number of files */
static int version = 0;		/* If 1, print only program name and version */
static int wave = 0;		/* If 1, print first dimension limits */
static int restwave = 0;	/* If 1, print first dimension limits */
static int printhead = 1;	/* 1 until header has been printed */

static char *RevMsg = "WCSHEAD WCSTools 3.6.4, 3 May 2006, Doug Mink SAO";

int
main (ac, av)
int ac;
char **av;
{
    char *str;
    int readlist = 0;
    int lfile;
    char *lastchar;
    char filename[256];
    FILE *flist;
    char *listfile;

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage();
	}

    /* crack arguments */
    for (av++; --ac > 0 && (*(str = *av)=='-' || *str == '@'); av++) {
	char c;
	if (*str == '@')
	    str = str - 1;
	while (c = *++str)
	switch (c) {

	case 'h':	/* hh:mm:ss output for crval, cdelt in arcsec/pix */
	    hms++;
	    break;

	case 'n':	/* hh:mm:ss output */
	    tabout++;
	    break;

	case 'r':	/* Print first dimension as rest wavelength, first and last values */
	    restwave++;
	    break;

	case 't':	/* tab table output */
	    tabout++;
	    break;

	case 'v':	/* more verbosity */
	    verbose++;
	    break;

	case 'w':	/* Print only first dimension, first and last values */
	    wave++;
	    break;

    	case 'z':	/* Use AIPS classic WCS */
    	    setdefwcs (WCS_ALT);
    	    break;

	case '@':	/* List of files to be read */
	    readlist++;
	    listfile = ++str;
	    str = str + strlen (str) - 1;
	    av++;
	    ac--;
	    break;

	default:
	    usage();
	    break;
	}
    }

    /* Read filenames of images from listfile */
    if (readlist) {

	/* Find maximimum filename length */
	if ((flist = fopen (listfile, "r")) == NULL) {
	    fprintf (stderr,"WCSHEAD: List file %s cannot be read\n",
		     listfile);
	    usage();
	    }
	nchar = 8;
	while (fgets (filename, 256, flist) != NULL) {
	    lfile = strlen (filename) - 1;
	    if (lfile > nchar) nchar = lfile;
	    }
	fclose (flist);

	/* Process files */
	if ((flist = fopen (listfile, "r")) == NULL) {
	    fprintf (stderr,"WCSHEAD: List file %s cannot be read\n",
		     listfile);
	    usage();
	    }
	while (fgets (filename, 256, flist) != NULL) {
	    lastchar = filename + strlen (filename) - 1;
	    *lastchar = 0;
	    ListWCS (filename);
	    if (verbose)
		printf ("\n");
	    }
	fclose (flist);
	}

    /* If no arguments left, print usage */
    if (ac == 0)
	usage();

    if (verbose)
	fprintf (stderr,"%s\n",RevMsg);

    nf = 0;
    while (ac-- > 0) {
	char *fn = *av++;
	nf++;
	ListWCS (fn);
	}

    return (0);
}

static void
usage ()
{
    fprintf (stderr,"%s\n",RevMsg);
    if (version)
	exit (-1);
    fprintf (stderr,"Print WCS part of FITS or IRAF image header\n");
    fprintf (stderr,"usage: wcshead [-htv] file.fit ...\n");
    fprintf (stderr,"  -h: Print CRVALs as hh:mm:ss dd:mm:ss\n");
    fprintf (stderr,"  -r: Print first dimension as rest wavelength limiting values\n");
    fprintf (stderr,"  -t: Print tab table output\n");
    fprintf (stderr,"  -v: Verbose\n");
    fprintf (stderr,"  -w: Print only first dimension limiting values\n");
    fprintf (stderr,"  -z: Use AIPS classic WCS subroutines\n");
    exit (1);
}

static void
ListWCS (filename)

char	*filename;	/* FITS or IRAF image file name */
{
    double w1, w2, dx, dy;
    int i, nxpix, nypix;
    char str[256], temp[256], *header;
    char *GetFITShead();
    char rastr[32], decstr[32], fform[8];
    struct WorldCoor *wcs, *GetWCSFITS();
    double wlast;

    wcs = GetWCSFITS (filename, verbose);
    if (nowcs (wcs)) {
	wcsfree (wcs);
	wcs = NULL;
	header = GetFITShead (filename, verbose);
	}

    if (wcs && wcs->ctype[0][0] == (char) 0) {
	wcsfree (wcs);
	wcs = NULL;
	header = GetFITShead (filename, verbose);
	}
    if (tabout && printhead) {
	strcpy (str, "filename");
	for (i = 1; i < nchar - 7; i++) strcat (str, " ");
	if (wave || restwave) {
	    strcat (str, "	naxis1	ctype1	first    ");
	    strcat (str, "	last     	delt   \n");
	    for (i = 0; i < nchar; i++) strcat (str, "-");
	    strcat (str, "	------	------");
	    strcat (str, "	---------	---------	-------\n");
	    }
	else {
	    strcat (str, "	naxis1	naxis2");
	    strcat (str, "	ctype1  	ctype2  ");
	    strcat (str, "	crval1 	crval2 	system  ");
	    strcat (str, "	crpix1	crpix2");
	    strcat (str, "	cdelt1	cdelt2");
	    strcat (str, "	crota2\n");
	    for (i = 0; i < nchar; i++) strcat (str, "-");
	    strcat (str, "	------	------");
	    strcat (str, "	--------	--------");
	    strcat (str, "	-------	-------");
	    strcat (str, "	--------	-------");
	    strcat (str, "	-------	-------");
	    strcat (str, "	-------	-------\n");
	    }
	printf ("%s", str);
	printhead = 0;
	}

    sprintf (fform,"%%%d.%ds",nchar, nchar);
    if (tabout) {
	sprintf (str, fform, filename);
	if (!wcs) {
	    if (wave) {
		strcat (str, "	___	___	_________	_________	_______");
		}
	    else if (restwave) {
		strcat (str, "	___	___	_________	_________	_______");
		}
	    else {
		hgeti4 (header, "NAXIS1", &nxpix);
		hgeti4 (header, "NAXIS2", &nypix);
		sprintf (temp, "	%d	%d", nxpix, nypix);
		strcat (str, temp);
		strcat (str, "	________	________	_______	_______	________	_______	_______	_______	_______	_______");
		}
	    }
	else if (wave) {
	    wlast = wcs->xref + ((wcs->nxpix - 1.0) * wcs->cdelt[0]);
	    sprintf (temp, "	%.0f	%s	%9.4f	%9.4f	%7.4f",
		     wcs->nxpix, wcs->ctype[0], wcs->xref, wlast, wcs->xinc);
	    strcat (str, temp);
	    }
	else if (restwave) {
	    w1 = wcs->xref / (1.0 + wcs->zvel);
	    wlast = wcs->xref + ((wcs->nxpix - 1.0) * wcs->cdelt[0]);
	    w2 = wlast / (1.0 + wcs->zvel);
	    sprintf (temp, "	%.0f	%s	%9.4f	%9.4f	%7.4f",
		     wcs->nxpix, wcs->ctype[0], w1, w2, wcs->xinc);
	    strcat (str, temp);
	    }
	else {
	    sprintf (temp, "	%.0f	%.0f", wcs->nxpix, wcs->nypix);
	    strcat (str, temp);
	    if (strlen (wcs->ctype[0]) < 8)
		sprintf (temp, "	%8.8s	%8.8s", wcs->ctype[0], wcs->ctype[1]);
	    else
		sprintf (temp, "	%s	%s", wcs->ctype[0], wcs->ctype[1]);
	    strcat (str, temp);
	    if (hms) {
		if (wcs->coorflip) {
			ra2str (rastr, 32, wcs->yref, ndec);
			dec2str (decstr, 32, wcs->xref, ndec-1);
			}
		else {
			ra2str (rastr, 32, wcs->xref, ndec);
			dec2str (decstr, 32, wcs->yref, ndec-1);
			}
		sprintf (temp, " %s %s %s", rastr, decstr, wcs->radecsys);
		}
	    else
		sprintf (temp, "	%7.2f	%7.2f	%8s",
			 wcs->xref, wcs->yref, wcs->radecsys);
	    strcat (str, temp);
	    sprintf (temp, "	%7.2f	%7.2f", wcs->xrefpix, wcs->yrefpix);
	    strcat (str, temp);
	    dx = 3600.0 * wcs->xinc;
	    dy = 3600.0 * wcs->yinc;
	    if (dx >= 10000.0 || dx <= -1000.0)
		sprintf (temp, "	%7.1f	%7.1f", dx, dy);
	    else if (dx >= 1000.0 || wcs->xinc <= -100.0)
		sprintf (temp, "	%7.2f	%7.2f", dx, dy);
	    else if (dx >= 100.0 || dx <= -10.0)
		sprintf (temp, "	%7.3f	%7.3f", dx, dy);
	    else
		sprintf (temp, "	%7.4f	%7.4f", dx, dy);
	    strcat (str, temp);
	    sprintf (temp, "	%7.4f", wcs->rot);
	    strcat (str, temp);
	    }
	}
    else {
	sprintf (str, fform, filename);
	if (!wcs) {
	    if (wave) {
		strcpy (temp, " ___ ___ _________ _________ _______");
		}
	    else if (restwave) {
		strcpy (temp, " ___ ___ _________ _________ _______");
		}
	    else {
		hgeti4 (header, "NAXIS1", &nxpix);
		hgeti4 (header, "NAXIS2", &nypix);
		sprintf (temp, " %4d %4d", nxpix, nypix);
		strcat (str, temp);
		strcpy (temp, " ________ ________ _______ _______ ________ _______ _______ _______ _______ _______");
		}
	    strcat (str, temp);
	    }
	else if (wave) {
	    sprintf (temp, " %.0f %s", wcs->nxpix, wcs->ctype[0]);
	    strcat (str, temp);
	    wlast = wcs->xref + ((wcs->nxpix - 1.0) * wcs->cdelt[0]);
	    sprintf (temp, " %9.4f %9.4f", wcs->xref, wlast);
	    strcat (str, temp);
	    sprintf (temp, " %7.4f", wcs->xinc);
	    strcat (str, temp);
	    }
	else if (restwave) {
	    w1 = wcs->xref / (1.0 + wcs->zvel);
	    wlast = wcs->xref + ((wcs->nxpix - 1.0) * wcs->cdelt[0]);
	    w2 = wlast / (1.0 + wcs->zvel);
	    sprintf (temp, " %.0f %s %9.4f %9.4f %7.4f",
		     wcs->nxpix, wcs->ctype[0], w1, w2, wcs->xinc);
	    strcat (str, temp);
	    }
	else {
	    sprintf (temp, " %4.0f %4.0f", wcs->nxpix, wcs->nypix);
	    strcat (str, temp);
	    if (strlen (wcs->ctype[0]) < 8)
		sprintf (temp, " %8.8s %8.8s", wcs->ctype[0], wcs->ctype[1]);
	    else
		sprintf (temp, " %s %s", wcs->ctype[0], wcs->ctype[1]);
	    strcat (str, temp);
	    if (hms) {
		if (wcs->coorflip) {
		    ra2str (rastr, 32, wcs->yref, ndec);
		    dec2str (decstr, 32, wcs->xref, ndec-1);
		    }
		else {
		    ra2str (rastr, 32, wcs->xref, ndec);
		    dec2str (decstr, 32, wcs->yref, ndec-1);
		    }
		sprintf (temp, " %s %s %8s", rastr, decstr, wcs->radecsys);
		}
	    else
		sprintf (temp, " %7.2f %7.2f %8s",
			 wcs->xref, wcs->yref, wcs->radecsys);
	    strcat (str, temp);
	    sprintf (temp, " %7.2f %7.2f", wcs->xrefpix, wcs->yrefpix);
	    strcat (str, temp);
	    dx = 3600.0 * wcs->xinc;
	    dy = 3600.0 * wcs->yinc;
	    if (dx >= 10000.0 || dx <= -1000.0)
		sprintf (temp, " %7.1f %7.1f", dx, dy);
	    else if (dx >= 1000.0 || wcs->xinc <= -100.0)
		sprintf (temp, " %7.2f %7.2f", dx, dy);
	    else if (dx >= 100.0 || dx <= -10.0)
		sprintf (temp, " %7.3f %7.3f", dx, dy);
	    else
		sprintf (temp, " %7.4f %7.4f", dx, dy);
	    strcat (str, temp);
	    sprintf (temp, " %7.4f", wcs->rot);
	    strcat (str, temp);
	    }
	}

    printf ("%s\n", str);

    wcsfree (wcs);
    wcs = NULL;

    return;
}
/* Feb 18 1998	New program
 * May 27 1998	Include fitsio.h instead of fitshead.h
 * Jun 24 1998	Add string lengths to ra2str() and dec2str() calls
 * Jul 10 1998	Add option to use AIPS classic WCS subroutines
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 * Nov 30 1998	Add version and help commands for consistency
 *
 * Apr  7 1999	Print lines all at once instead of one variable at a time
 * Jun  3 1999	Change PrintWCS to ListWCS to avoid name conflict
 * Oct 15 1999	Free wcs using wcsfree()
 * Oct 22 1999	Drop unused variables after lint
 * Nov 30 1999	Fix declaration of ListWCS()
 *
 * Jan 28 2000	Call setdefwcs() with WCS_ALT instead of 1
 * Jun 21 2000	Add -w option to print limits for 1-d WCS
 * Aug  4 2000	Add -w option to printed option list
 *
 * Apr  8 2002	Free wcs structure if no WCS is found in file header
 * May  9 2002	Add option to print rest wavelength limits
 * May 13 2002	Set wcs pointer to NULL after freeing data structure
 * Jun 19 2002	Add verbose argument to GetWCSFITS()
 *
 * Jul 19 2004	Print header if flag is set, not if first file
 * Jul 19 2004	Print underscores to fill lines if no WCS
 */
