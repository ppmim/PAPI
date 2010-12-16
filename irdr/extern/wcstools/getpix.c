/* File getpix.c
 * July 29, 2005
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics)
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

static void usage();
static void PrintPix();
static void procpix();

const char *RevMsg = "GETPIX WCSTools 3.6.4, 3 May 2006, Doug Mink SAO";

static int verbose = 0;		/* verbose/debugging flag */
static int version = 0;		/* If 1, print only program name and version */
static int nline = 10;		/* Number of pixels printer per line */
static char *pform = NULL;	/* Format in which to print pixels */
static int pixlabel = 0;	/* If 1, label pixels in output */
static int gtcheck = 0;		/* If 1, list pixels greater than gtval */
static int ltcheck = 0;		/* If 1, list pixels less than ltval */
static int nopunct=0;		/* If 1, print output with no punctuation */
static int printrange = 0;	/* If 1, print range of values, not values */
static int printmean = 0;	/* If 1, print mean of values, not values*/
static double gtval = 0.0;
static double ltval = 0.0;
static double ra0 = -99.0;	/* Initial center RA in degrees */
static double dec0 = -99.0;	/* Initial center Dec in degrees */
static double rad0 = 0.0;	/* Search box radius */
static double dra0 = 0.0;	/* Search box width */
static double ddec0 = 0.0;	/* Search box height */
static double eqcoor = 2000.0;  /* Equinox of search center */
static int syscoor = 0;         /* Input search coordinate system */

int
main (ac, av)

int ac;
char **av;
{
    char *str;
    char *fn;
    char *rrange;       /* Row range string */
    char *crange;       /* Column range string */
    char *rstr, *dstr, *cstr, *ccom;
    int systemp;
    int i;
    int npix = 0;
    int maxnpix = 10;
    int ix, iy;
    int *xpix, *ypix;
    int nmaxpix, *xp1, *yp1;

    xpix = calloc (maxnpix, sizeof (int));
    ypix = calloc (maxnpix, sizeof (int));

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage();
	}

    crange = NULL;
    rrange = NULL;
    fn = NULL;

    /* crack arguments */
    for (av++; --ac > 0; av++) {
	str = *av;

	/* output format */
	if (str[0] == '%') {
	    pform = *av;
	    }

	/* other command */
	else if (str[0] == '-') {
	    char c;
	    while (c = *++str)
	    switch (c) {

		case 'v':	/* more verbosity */
		    verbose++;
		    break;

		case 'd':	/* Print range of values */
		    printrange=1;
		    break;

		case 'g':	/* Keep pixels greater than this */
		    if (ac < 2)
			usage();
		    gtval = atof (*++av);
		    gtcheck++;
		    ac--;
		    break;

		case 'l':	/* Keep pixels less than this */
		    if (ac < 2)
			usage();
		    ltval = atof (*++av);
		    ltcheck++;
		    ac--;
		    break;

		case 'm':	/* Print mean, sigma of values */
		    printmean = 1;
		    break;

		case 'n':	/* Number of pixels per line */
		    if (ac < 2)
			usage();
		    nline = atoi (*++av);
		    ac--;
		    break;

		case 'p':	/* label pixels */
		    pixlabel++;
		    break;

		case 'r':	/* Box radius in arcseconds */
		    if (ac < 2)
    			usage ();
		    av++;
		    if ((dstr = strchr (*av, ',')) != NULL) {
			*dstr = (char) 0;
			dstr++;
			}
		    if (strchr (*av,':'))
			rad0 = 3600.0 * str2dec (*av);
		    else
			rad0 = atof (*av);
		    if (dstr != NULL) {
			dra0 = rad0;
			rad0 = 0.0;
			if (strchr (dstr, ':'))
			    ddec0 = 3600.0 * str2dec (dstr);
			else
			    ddec0 = atof (dstr);
			if (ddec0 <= 0.0)
			    ddec0 = dra0;
			/* rad0 = sqrt (dra0*dra0 + ddec0*ddec0); */
			}
    		    ac--;
    		    break;

		case 's':	/* Print x y value without punctuation */
		    nopunct++;
		    break;

		default:
		    usage();
		    break;
		}
	    }

	/* Set search RA, Dec, and equinox if colon in argument */
	else if (strsrch (*av,":") != NULL) {
	    if (ac < 2)
		usage ();
	    else {
		strcpy (rstr, *av);
		ac--;
		strcpy (dstr, *++av);
		ra0 = str2ra (rstr);
		dec0 = str2dec (dstr);
		ac--;
		if (ac < 1) {
		    syscoor = WCS_J2000;
		    eqcoor = 2000.0;
		    }
		else if ((syscoor = wcscsys (*(av+1))) >= 0)
		    eqcoor = wcsceq (*++av);
		else {
		    syscoor = WCS_J2000;
		    eqcoor = 2000.0;
		    }
		}
	    }

	/* Center coordinates in degrees */
	else if (isnum (str) == 2 && ac > 1 && isnum (*(av+1)) == 2) {
	    rstr = *av++;
	    ac--;
	    dstr = *av;
	    ac--;
	    ra0 = atof (rstr);
	    dec0 = atof (dstr);
	    if (ac > 0 && (systemp = wcscsys (*av)) > 0) {
		cstr = *av++;
		syscoor = systemp;
		eqcoor = wcsceq (cstr);
		}
	    else {
		cstr = (char *) malloc (8);
		strcpy (cstr, "J2000");
		syscoor = WCS_J2000;
		eqcoor = 2000.0;
		}
	    }

	/* Coordinate pairs for pixels to print */
        else if (isnum (str) && ac > 1 && isnum (*(av+1))) {
	    if (npix+1 > maxnpix) {
		nmaxpix = 2 * maxnpix;
		xp1 = calloc (nmaxpix, sizeof (int));
		yp1 = calloc (nmaxpix, sizeof (int));
		for (i = 0; i < nmaxpix; i++) {
		    xp1[i] = xpix[i];
		    yp1[i] = ypix[i];
		    }
		free (xpix);
		free (ypix);
		xpix = xp1;
		ypix = yp1;
		}
	    ix = atoi (*av);
	    iy = atoi (*(av+1));
	    if (ix == 0 || iy == 0) {
		crange = *av++;
		rrange = *av;
		}
	    else {
	        xpix[npix] = ix;
	        ypix[npix] = iy;
	        npix++;
		}
	    ac--;
	    }

	/* Range of pixels to print (only one allowed) */
        else if (isrange (str) && ac > 1 && isrange (*(av+1))) {
	    crange = *av++;
	    ac--;
	    rrange = *av;
	    }

	/* file name */
	else
	    fn = str;

	}

    if (fn && ((crange && rrange) || npix > 0))
        PrintPix (fn, crange, rrange, npix, xpix, ypix);

    return (0);
}

static void
usage ()
{
    fprintf (stderr,"%s\n",RevMsg);
    if (version)
	exit (-1);
    fprintf (stderr,"Print FITS or IRAF pixel values\n");
    fprintf(stderr,"Usage: getpix [-vp][-n num][-g val][-l val][format] file.fit x_range y_range\n");
    fprintf(stderr,"  or   getpix [-vp][-n num][-g val][-l val][format] file.fit x1 y1 x2 y2 ... xn yn\n");
    fprintf(stderr,"  format: C-style (%%f, %%d, ...) format for pixel values\n");
    fprintf(stderr,"  -d: Print range of pixel values in specified image region\n");
    fprintf(stderr,"  -f name: Write specified region to a FITS file\n");
    fprintf(stderr,"  -g num: keep pixels with values greater than this\n");
    fprintf(stderr,"  -l num: keep pixels with values less than this\n");
    fprintf(stderr,"  -m: Print mean of pixel values in specified image region\n");
    fprintf(stderr,"  -n num: number of pixel values printed per line\n");
    fprintf(stderr,"  -p: label pixels\n");
    fprintf(stderr,"  -r num: radius (<0=box) to extract in degrees/arcsec\n");
    fprintf(stderr,"  -s: print x y value with no punctuation\n");
    fprintf(stderr,"  -v: verbose\n");
    fprintf(stderr,"   %%: C format for each pixel value\n");
    exit (1);
}


static void
PrintPix (name, crange, rrange, npix, xpix, ypix)

char *name;
char *crange;		/* Column range string */
char *rrange;		/* Row range string */
int npix;		/* Number of coordinate pairs */
int *xpix, *ypix;	/* Vectors of x,y coordinate pairs */

{
    char *header;	/* FITS image header */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    char *irafheader;	/* IRAF image header */
    char *image;	/* FITS or IRAF image */
    double bzero;	/* Zero point for pixel scaling */
    double bscale;	/* Scale factor for pixel scaling */
    int iraffile;
    double dpix, dsum, dmean, dmin, dmax, dnpix;
    char *c;
    int *yi;
    int bitpix,xdim,ydim, ipix, i, nx, ny, ix, iy, x, y;
    char pixname[255];
    char nform[8];
    struct Range *xrange;    /* X range structure */
    struct Range *yrange;    /* Y range structure */

    /* Open IRAF image if .imh extension is present */
    if (isiraf (name)) {
	iraffile = 1;
	if ((irafheader = irafrhead (name, &lhead)) != NULL) {
	    header = iraf2fits (name, irafheader, lhead, &nbhead);
	    free (irafheader);
	    if (header == NULL) {
		fprintf (stderr, "Cannot translate IRAF header %s/n",name);
		return;
		}
	    if ((image = irafrimage (header)) == NULL) {
		hgetm (header,"PIXFIL", 255, pixname);
		fprintf (stderr, "Cannot read IRAF pixel file %s\n", pixname);
		free (irafheader);
		free (header);
		return;
		}
	    }
	else {
	    fprintf (stderr, "Cannot read IRAF file %s\n", name);
	    return;
	    }
	}

    /* Open FITS file if .imh extension is not present */
    else {
	iraffile = 0;
	if ((header = fitsrhead (name, &lhead, &nbhead)) != NULL) {
	    if ((image = fitsrimage (name, nbhead, header)) == NULL) {
		fprintf (stderr, "Cannot read FITS image %s\n", name);
		free (header);
		return;
		}
	    }
	else {
	    fprintf (stderr, "Cannot read FITS file %s\n", name);
	    return;
	    }
	}
    if (verbose) {
	fprintf (stderr,"%s\n",RevMsg);
	if (npix > 0)
	    fprintf (stderr,"Print pixels in ");
	else if (!strcmp (crange, "0"))
	    fprintf (stderr,"Print rows %s in ", rrange);
	else if (!strcmp (rrange, "0"))
	    fprintf (stderr,"Print columns %s in ", crange);
	else
	    fprintf (stderr,"Print rows %s, columns %s in ", rrange, crange);
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s\n", name);
	else
	    fprintf (stderr,"FITS image file %s\n", name);
	if (ltcheck & gtcheck)
	    fprintf (stderr, "%f < pixel values < %f\n", gtval, ltval);
	else if (ltcheck)
	    fprintf (stderr, "pixel values < %f\n", ltval);
	else if (gtcheck)
	    fprintf (stderr, "pixel values > %f\n", gtval);
	}

/* Get value of specified pixel */

    /* Get size of image and scaling factors */
    hgeti4 (header,"BITPIX",&bitpix);
    xdim = 1;
    hgeti4 (header,"NAXIS1",&xdim);
    ydim = 1;
    hgeti4 (header,"NAXIS2",&ydim);
    bzero = 0.0;
    hgetr8 (header,"BZERO",&bzero);
    bscale = 1.0;
    hgetr8 (header,"BSCALE",&bscale);

    /* Set initial values */
    dsum = 0.0;
    dnpix = 0.0;
    dmin = 0.0;
    dmax = 0.0;

    /* Set format if not already set */
    if (pform == NULL) {
	pform = (char *) calloc (8,1);
	if (bitpix > 0)
	    strcpy (pform, "%d");
	else
	    strcpy (pform, "%.2f");
	}

/* Print values at specified coordinates in an image  */
    if (npix > 0) {

	/* Loop through rows starting with the last one */
	for (i = 0; i < npix; i++) {
            dpix = getpix1(image,bitpix,xdim,ydim,bzero,bscale,xpix[i],ypix[i]);
	    if (gtcheck || ltcheck) {
		if (gtcheck && dpix > gtval ||
		    ltcheck && dpix < ltval) {
		    procpix (&dsum, &dnpix, &dmin, &dmax, dpix);
		    if (nopunct)
			printf ("%d %d %f\n", xpix[i], ypix[i], dpix);
		    else
			printf ("[%d,%d] = %f\n", xpix[i], ypix[i], dpix);
		    }
		continue;
		}
	    else
		procpix (&dsum, &dnpix, &dmin, &dmax, dpix);
	    if (printrange || printmean)
		continue;
	    if (bitpix > 0) {
		if ((c = strchr (pform,'f')) != NULL)
		    *c = 'd';
		if (dpix > 0)
	 	    ipix = (int) (dpix + 0.5);
		else if (dpix < 0)
		     ipix = (int) (dpix - 0.5);
		else
		    ipix = 0;
		}
	    else {
		if ((c = strchr (pform,'d')) != NULL)
		    *c = 'f';
		}
	    if (verbose) {
		printf ("%s[%d,%d] = ",name,xpix[i],ypix[i]);
		if (bitpix > 0)
		    printf (pform, ipix);
		else
		    printf (pform, dpix);
		printf ("\n");
		}
	    else {
		if (bitpix > 0)
		    printf (pform, ipix);
		else
		    printf (pform, dpix);
		if ((i+1) % nline == 0)
		    printf ("\n");
		else
		    printf (" ");
		}
	    }
	if (!verbose && !ltcheck && !gtcheck)
	    printf ("\n");
	free (xpix);
	free (ypix);
	}

/* Print entire image */
    else if (!strcmp (rrange, "0") && !strcmp (crange, "0")) {
	if (printmean || printrange) {
	    nx = xdim;
	    ny = ydim;
	    for (y = 0; y < ny; y++) {
		for (x = 0; x < nx; x++) {
        	    dpix = getpix (image,bitpix,xdim,ydim,bzero,bscale,x,y);
		    if (!gtcheck && !ltcheck)
			procpix (&dsum, &dnpix, &dmin, &dmax, dpix);
		    else if (gtcheck && ltcheck) {
			if (dpix > gtval && dpix < ltval)
			    procpix (&dsum, &dnpix, &dmin, &dmax, dpix);
			}
		    else if (gtcheck && dpix > gtval ||
			ltcheck && dpix < ltval)
			procpix (&dsum, &dnpix, &dmin, &dmax, dpix);
		    }
		}
	    }
	else if (gtcheck || ltcheck) {
	    nx = xdim;
	    ny = ydim;
	    for (y = 0; y < ny; y++) {
		for (x = 0; x < nx; x++) {
        	    dpix = getpix (image,bitpix,xdim,ydim,bzero,bscale,x,y);
		    if (gtcheck && ltcheck) {
			if (dpix > gtval && dpix < ltval) {
			    if (nopunct)
				printf ("%d %d %f\n", x+1, y+1, dpix);
			    else
				printf ("[%d,%d] = %f\n", x+1, y+1, dpix);
			    }
			}
		    else if (gtcheck && dpix > gtval ||
			ltcheck && dpix < ltval) {
			if (nopunct)
			    printf ("%d %d %f\n", x+1, y+1, dpix);
			else
			    printf ("[%d,%d] = %f\n", x+1, y+1, dpix);
			}
		    }
		}
	    }
	else
	    printf ("GETPIX will not print this %d x %d image; use ranges\n",
		xdim, ydim);
	}

/* Print entire columns */
    else if (!strcmp (rrange, "0")) {
	xrange = RangeInit (crange, xdim);
	nx = rgetn (xrange);
	ny = ydim;
	for (ix = 0; ix < nx; ix++) {
	    rstart (xrange);
	    x = rgeti4 (xrange) - 1;
	    for (y = 0; y < ny; y++) {
        	dpix = getpix (image,bitpix,xdim,ydim,bzero,bscale,x,y);
		if (gtcheck || ltcheck) {
		    if (gtcheck && dpix > gtval ||
			ltcheck && dpix < ltval) {
			procpix (&dsum, &dnpix, &dmin, &dmax, dpix);
			if (nopunct)
			    printf ("%d %d %f\n", x+1, y+1, dpix);
			else
			    printf ("[%d,%d] = %f\n", x+1, y+1, dpix);
			}
		    continue;
		    }
		else
		    procpix (&dsum, &dnpix, &dmin, &dmax, dpix);
	        if (bitpix > 0) {
		    if (dpix > 0)
	 		ipix = (int) (dpix + 0.5);
		    else if (dpix < 0)
		 	ipix = (int) (dpix - 0.5);
		    else
			ipix = 0;
		    }
		if (verbose) {
		    printf ("%s[%d,%d] = ",name,x,y);
		    if (bitpix > 0)
			printf (pform, ipix);
		    else
			printf (pform, dpix);
		    printf ("\n");
		    }
		else {
		    if (bitpix > 0)
			printf (pform, ipix);
		    else
			printf (pform, dpix);
		    if ((y+1) % nline == 0)
			printf ("\n");
		    else
			printf (" ");
		    }
		}
	    if (y % nline != 0 && !gtcheck && !ltcheck)
		printf ("\n");
	    if (nx > 1 && !gtcheck && !ltcheck)
		printf ("\n");
	    }
	free (xrange);
	}

/* Print entire rows */
    else if (!strcmp (crange, "0")) {
	yrange = RangeInit (rrange, xdim);
	ny = rgetn (yrange);
	nx = xdim;
	for (iy = 0; iy < ny; iy++) {
	    y = rgeti4 (yrange) - 1;
	    for (x = 0; x < nx; x++) {
        	dpix = getpix (image,bitpix,xdim,ydim,bzero,bscale,x,y);
		if (gtcheck || ltcheck) {
		    if (gtcheck && dpix > gtval ||
			ltcheck && dpix < ltval) {
			procpix (&dsum, &dnpix, &dmin, &dmax, dpix);
			if (nopunct)
			    printf ("%d %d %f\n", x+1, y+1, dpix);
			else
			    printf ("[%d,%d] = %f\n", x+1, y+1, dpix);
			}
		    continue;
		    }
		else
		    procpix (&dsum, &dnpix, &dmin, &dmax, dpix);
	        if (bitpix > 0) {
		    if (dpix > 0)
	 		ipix = (int) (dpix + 0.5);
		    else if (dpix < 0)
		 	ipix = (int) (dpix - 0.5);
		    else
			ipix = 0;
		    }
		if (verbose) {
		    printf ("%s[%d,%d] = ",name,x,y);
		    if (bitpix > 0)
			printf (pform, ipix);
		    else
			printf (pform, dpix);
		    printf ("\n");
		    }
		else {
		    if (bitpix > 0)
			printf (pform, ipix);
		    else
			printf (pform, dpix);
		    if ((x+1) % nline == 0)
			printf ("\n");
		    else
			printf (" ");
		    }
		}
	    if (x % nline != 0 && !gtcheck && !ltcheck)
		printf ("\n");
	    if (ny > 1 && !gtcheck && !ltcheck)
		printf ("\n");
	    }
	free (yrange);
	}

/* Print a region of a two-dimensional image */
    else {
	xrange = RangeInit (crange, xdim);
	nx = rgetn (xrange);

	/* Make list of y coordinates */
	yrange = RangeInit (rrange, ydim);
	ny = rgetn (yrange);
	yi = (int *) calloc (ny, sizeof (int));
	for (i = 0; i < ny; i++) {
	    yi[i] = rgeti4 (yrange) - 1;
	    }

	/* Label horizontal pixels */
	if (!verbose && pixlabel) {
	    printf ("     ");
	    rstart (xrange);
	    strcpy (nform, pform);
	    if ((c = strchr (nform,'.')) != NULL) {
		*c = 'd';
		c[1] = (char) 0;
		}
	    else if ((c = strchr (nform,'f')) != NULL) {
		*c = 'd';
		}
	    for (ix = 0; ix < nx; ix++) {
		x = rgeti4 (xrange);
		printf (" ");
		printf (nform, x);
		}
	    printf ("\n");
	    }
	if (verbose)
	    iy = -1;
	else
	    iy = ny;

	/* Loop through rows starting with the last one */
	for (i = 0; i < ny; i++) {
	    if (verbose)
		iy++;
	    else
		iy--;
	    rstart (xrange);
	    if (!verbose && pixlabel)
		printf ("%4d: ",yi[iy]+1);

	    /* Loop through columns */
	    for (ix = 0; ix < nx; ix++) {
		x = rgeti4 (xrange) - 1;
        	dpix = getpix (image,bitpix,xdim,ydim,bzero,bscale,x,yi[iy]);
		if (gtcheck || ltcheck) {
		    if (gtcheck && dpix > gtval ||
			ltcheck && dpix < ltval) {
			procpix (&dsum, &dnpix, &dmin, &dmax, dpix);
			if (nopunct)
			    printf ("%d %d %f\n", x+1, yi[iy]+1, dpix);
			else
			    printf ("[%d,%d] = %f\n", x+1, yi[iy]+1, dpix);
			}
		    continue;
		    }
		else
		    procpix (&dsum, &dnpix, &dmin, &dmax, dpix);
	        if (bitpix > 0) {
		    if ((c = strchr (pform,'f')) != NULL)
			*c = 'd';
		    if (dpix > 0)
	 		ipix = (int) (dpix + 0.5);
		    else if (dpix < 0)
		 	ipix = (int) (dpix - 0.5);
		    else
			ipix = 0;
		    }
		else {
		    if ((c = strchr (pform,'d')) != NULL)
			*c = 'f';
		    }
		if (verbose) {
		    printf ("%s[%d,%d] = ",name,x+1,yi[i]+1);
		    if (bitpix > 0)
			printf (pform, ipix);
		    else
			printf (pform, dpix);
		    printf ("\n");
		    }
		else {
		    if (bitpix > 0)
			printf (pform, ipix);
		    else
			printf (pform, dpix);
		    if ((ix+1) % nline == 0)
			printf ("\n");
		    else
			printf (" ");
		    }
		}
	    if (!verbose && !ltcheck && !gtcheck)
		printf ("\n");
	    }
	free (xrange);
	free (yrange);
	}

    if (printmean) {
	dmean = dsum / dnpix;
	printf ("Mean= %.4f ", dmean);
	}
    if (printrange)
	printf ("Range = %.4f - %.4f ", dmin, dmax);
    if (printmean || printrange)
	printf ("for %d pixels\n", dnpix);

    free (header);
    free (image);
    return;
}


static void
procpix (dsum, dmin, dmax, dpix)

double	*dsum;	/* Sum of pixel values */
double	*dmin;	/* Minimum pixel value */
double	*dmax;	/* Maximum pixel value */
double	dpix;	/* Current pixel value */
{
	*dsum = *dsum + dpix;
	if (*dmin == 0.0 && *dmax == 0.0) {
	    *dmin = dpix;
	    *dmax = dpix;
	    }
	else if (dpix < *dmin)
	    *dmin = dpix;
	else if (dpix > *dmax)
	    *dmax = dpix;
}
/* Dec  6 1996	New program
 *
 * Feb 21 1997  Check pointers against NULL explicitly for Linux
 * Dec 15 1997	Add capability of reading and writing IRAF 2.11 images
 *
 * May 27 1998	Include fitsio.h instead of fitshead.h
 * Jul 24 1998	Make irafheader char instead of int
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 * Oct 14 1998	Use isiraf() to determine file type
 * Nov 30 1998	Add version and help commands for consistency
 *
 * Feb 12 1999	Initialize dxisn to 1 so it works for 1-D images
 * Apr 29 1999	Add BZERO and BSCALE
 * Jun 29 1999	Fix typo in BSCALE setting
 * Jul  2 1999	Use ranges instead of individual pixels
 * Oct 15 1999	Fix format statement
 * Oct 22 1999	Drop unused variables after lint
 * Dec  9 1999	Add -g -l limits
 * Dec 13 1999	Fix bug so that -g and -l limits can be ANDed
 *
 * Mar 23 2000	Use hgetm() to get the IRAF pixel file name, not hgets()
 *
 * Jan 30 2001	Fix format specification in help message
 *
 * Jun  3 2002	Add -s option to print x y value with no punctuation
 * Oct 30 2002	Add code to count lines when printing a region
 *
 * Feb 20 2003	Add option to enter multiple pixel (x,y) as well as ranges
 * Mar 26 2003	Fix pixel counter bug in individual pixel printing
 * Sep 17 2003	Fix bug which broke use of 0 as substitute for 1-naxisn range
 *
 * Apr 26 2004	Fix handling of 0 0 for entire image
 * Aug 30 2004	Fix declarations
 * Sep 21 2004	Fix bug which used x instead of ix for number of elements printed
 *
 * Jul 29 2005	Add mean and range computation
 */
