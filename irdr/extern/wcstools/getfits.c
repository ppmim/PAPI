/* File getfits.c
 * September 30, 2005
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
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"

#define MAXRANGE 20
#define MAXFILES 2000
static int maxnfile = MAXFILES;

#define MAXKWD 500
static int maxnkwd = MAXKWD;

/* Structure for dealing with ranges */
struct Range {
    double first;	/* Current minimum value */
    double last;	/* Current maximum value */
    double step;	/* Current step in value */
    double value;	/* Current value */
    double ranges[MAXRANGE*3];  /* nranges sets of first, last, step */
    int nvalues;	/* Total number of values in all ranges */
    int nranges;	/* Number of ranges */
    int irange;		/* Index of current range */
};

/* Subroutines for dealing with ranges */
static struct Range *RangeInit(); /* Initialize range structure from string */
static int isrange();	/* Return 1 if string is a range of numbers, else 0 */
static int rgetn();	/* Return number of values in all ranges */
static int rgeti4();	/* Return next number in range as integer */
static double rgetr8();	/* Return next number in range as double */
static void rstart();	/* Restart range */

static void usage();
static void nextname();	/* Find next available name (namea, nameb, ...) */
static int ExtractFITS();

static int verbose = 0;		/* verbose/debugging flag */
static char *RevMsg = "GETFITS WCSTools 3.6.4, 3 May 2006, Doug Mink (dmink@cfa.harvard.edu)";
static int version = 0;		/* If 1, print only program name and version */
static char outname[128];	/* Name for output image */
static char outdir[256];	/* Output directory pathname */
static int first = 1;
static int nlog = 100;
static int xcpix = 0;
static int ycpix = 0;
static int xdpix = 0;
static int ydpix = 0;
static char *rrange;	/* Row range string */
static char *crange;	/* Column range string */
static double ra0 = -99.0;	/* Initial center RA in degrees */
static double dec0 = -99.0;	/* Initial center Dec in degrees */
static double rad0 = 0.0;	/* Search box radius */
static double dra0 = 0.0;	/* Search box width */
static double ddec0 = 0.0;	/* Search box height */
static double epoch0 = 0.0;	/* Epoch for coordinates */
static int syscoor = 0;		/* Input search coordinate system */
static double eqcoor = 0.0;	/* Input search coordinate system */

int
main (ac, av)
int ac;
char **av;
{			
    char *str;
    char *listfile;
    char **fn;
    char filename[256];
    char temp[80];
    int ifile, nfile, nbytes;
    FILE *flist;
    char *nextarg;
    char rastr[32], decstr[32];
    int nkwd = 0;		/* Number of keywords to delete */
    int nkwd1 = 0;		/* Number of keywords in delete file */
    char **kwd;			/* List of keywords to be deleted */
    char **kwdnew;
    int ikwd;
    FILE *fdk;
    char *klistfile;

    crange = NULL;
    rrange = NULL;
    outname[0] = 0;
    outdir[0] = 0;
    nfile = 0;
    fn = (char **)calloc (maxnfile, sizeof(char *));
    klistfile = NULL;

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage("");
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage("");
	}

    /* crack arguments */
    for (av++; --ac > 0; av++) {

	/* Set center RA, Dec, and equinox if colon in argument */
	if (strsrch (*av,":") != NULL) {
	    if (ac < 2)
		usage ("Right ascension given but no declination");
	    else {
		strcpy (rastr, *av);
		strcpy (decstr, *++av);
		ra0 = str2ra (rastr);
		dec0 = str2dec (decstr);
		ac--;
		if (ac < 1) {
		    syscoor = WCS_J2000;
		    eqcoor = 2000.0;
		    }
		else if ((syscoor = wcscsys (*(av+1))) >= 0) {
		    eqcoor = wcsceq (*++av);
		    ac--;
		    }
		else {
		    syscoor = WCS_J2000;
		    eqcoor = 2000.0;
		    }
		xcpix = -1;
		ycpix = -1;
		}
	    }

	/* Read command */
	else if (*(str = *av) == '-') {
	    char c;
	    while (c = *++str) {
		switch (c) {

		case 'v':	/* more verbosity */
		    verbose++;
		    break;

		case 'c':	/* center coordinates for extraction */
		    if (ac < 3)
			usage ("-c needs two arguments");
		    xcpix = atoi (*++av);
		    ac--;
		    ycpix = atoi (*++av);
		    ac--;
		    break;

		case 'd':	/* output directory */
		    if (ac < 2)
			usage ("-d needs a directory");
		    strcpy (outdir, *++av);
		    ac--;
		    break;

		case 'k':	/* Read keywords to delete from this file */
		    if (ac < 2)
			usage ("-k needs a file of keywords");
		    klistfile = *av++;
		    nkwd1 = getfilelines (klistfile);
		    if (nkwd1 + nkwd > maxnkwd) {
			maxnkwd = maxnkwd + nkwd1;
			kwdnew = (char **)calloc (maxnkwd, sizeof(char *));
			if (nkwd > 0) {
			    for (ikwd = 0; ikwd < nkwd; ikwd++)
				kwdnew[ikwd] = kwd[ikwd];
			    free (kwd);
			    }
			kwd = kwdnew;
			}
		    if ((fdk = fopen (klistfile, "r")) == NULL) {
			fprintf (stderr,"GETFITS: File %s cannot be read\n",
				 klistfile);
			}
		    else {
			for (ikwd = nkwd; ikwd < nkwd+nkwd1; ikwd++) {
			    kwd[ikwd] = (char *) calloc (32, 1);
			    first_token (fdk, 31, kwd[ikwd]);
			    }
			fclose (fdk);
			}
		    ac--;
		    break;

		case 'i':	/* logging interval */
		    if (ac < 2)
			usage ("-i needs an argument");
		    nlog = atoi (*++av);
		    ac--;
		    break;

		case 'o':	/* output file name */
		    if (ac < 2)
			usage ("-o needs an argument");
		    strcpy (outname, *++av);
		    ac--;
		    break;

		case 's':	/* output to stdout */
		    strcpy (outname, "stdout");
		    break;

		case 'x':	/* width and height for extraction */
		    if (ac < 2)
			usage ("-x needs at least one argument");
		    xdpix = atoi (*++av);
		    ac--;
		    if (ac > 1) {
			nextarg = *(av+1);
			if (isnum (nextarg)) {
			    ydpix = atoi (nextarg);
			    av++;
			    ac--;
			    }
			else
			    ydpix = xdpix;
			}
		    else
			ydpix = xdpix;
		    break;

	        default:
		    sprintf (temp, "Illegal argument '%c'", c);
		    usage(temp);
		    break;
		}
		}
    	    }

        /* center or center and size of section to extract */
        else if (isnum (*av)) {
	    if (ac > 2 && isnum (*(av+1)) && (syscoor = wcscsys (*(av+2))) > 0) {
		ra0 = str2ra (*av++);
		ac--;
		dec0 = str2dec (*av++);
		ac--;
		eqcoor = wcsceq (*av);
		xcpix = -1;
		ycpix = -1;
		}
	    else {
		if (!xcpix)
		    xcpix = atoi (str);
		else if (!ycpix)
		    ycpix = atoi (str);
		else if (!xdpix)
		    xdpix = atoi (str);
		else if (!ydpix)
		    ydpix = atoi (str);
		}
	    }

        /* range of pixels to extract */
        else if (isrange (*av)) {
	    if (crange == NULL)
		crange = str;
	    else
		rrange = str;
	    }

 	else if (*av[0] == '@') {
	    listfile = *av + 1;
	    }

        /* Image file */
        else if (isfits (*av)) {
            if (nfile >= maxnfile) {
                maxnfile = maxnfile * 2;
                nbytes = maxnfile * sizeof (char *);
                fn = (char **) realloc ((void *)fn, nbytes);
                }
            fn[nfile] = *av;
            nfile++;
            }
	}

    /* If center is set, but not size, extract 500x500 image */
    if (!crange) {
	if (xcpix && ycpix && xdpix == 0)
	    xdpix = 500;
	if (xcpix && ycpix && ydpix == 0)
	    ydpix = xdpix;
	}

    /* If one side is set, set the other */
    if (xdpix && !ydpix)
	ydpix = xdpix;

    /* now there are ac remaining file names starting at av[0] */
    if (listfile && isimlist (listfile)) {
	nfile = getfilelines (listfile);
	if ((flist = fopen (listfile, "r")) == NULL) {
	    fprintf (stderr,"I2F: Image list file %s cannot be read\n",
		     listfile);
	    usage ();
	    }
	for (ifile = 0; ifile < nfile; ifile++) {
	    first_token (flist, 254, filename);
	    ExtractFITS (filename,kwd,nkwd);
	    }
	fclose (flist);
	}

    /* Process files from command line */
    else if (fn) {
	for (ifile = 0; ifile < nfile; ifile++) {
	    (void) ExtractFITS (fn[ifile],kwd,nkwd);
	    }
	}

    return (0);
}

static void
usage (errmsg)

char *errmsg;	/* Error message */
{
    fprintf (stderr,"%s\n",RevMsg);
    if (version)
	exit (-1);
    if (*errmsg)
	fprintf (stderr, "*** %s ***\n", errmsg);
    fprintf (stderr,"Extract FITS files from FITS image files\n");
    fprintf(stderr,"Usage: getfits -sv [-i num] [-o name] [-d dir] file.fits [xrange yrange] [x y dx [dy]] ...\n");
    fprintf(stderr,"  or : getfits [-sv1][-i num][-o name] [-d path] @fitslist [xrange yrange] [x y dx [dy]]\n");
    fprintf(stderr,"  xrange: Columns to extract in format x1-x2\n");
    fprintf(stderr,"  yrange: Rows to extract in format y1-y2\n");
    fprintf(stderr,"  x y: Center pixel (column row) of region to extract\n");
    fprintf(stderr,"  hh:mm:ss dd:mm:ss sys: Center pixel in sky coordintes\n");
    fprintf(stderr,"  dx dy: Width and height in pixels of region to extract\n");
    fprintf(stderr,"         (Height is same as width if omitted)\n");
    fprintf(stderr,"  -d dir: write FITS file(s) to this directory\n");
    fprintf(stderr,"  -i num: log rows as they are copied at this interval\n");
    fprintf(stderr,"  -k file: file of keyword names to delete from output header\n");
    fprintf(stderr,"  -o name: output name for one file\n");
    fprintf(stderr,"  -s: write output to standard output\n");
    fprintf(stderr,"  -x dx dy: dimensions of image section to be extracted\n");
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}

static int
ExtractFITS (name, kwd, nkwd)

char	*name;
char	**kwd;
int	nkwd;
{
    char *header;	/* FITS header */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    char pixname[256];	/* Pixel file name */
    char history[128];	/* for HISTORY line */
    char *filename;	/* Pointer to start of file name */
    char inname[256];	/* Name of TIFF file */
    char fitsname[256];	/* Name of FITS file */
    char fitspath[256];	/* Pathname of FITS file  */
    char *hdrfile1;
    char *fitsfile;
    char *image;        /* FITS image */
    int nbytes;
    double cra;         /* Center right ascension in degrees (returned) */
    double cdec;        /* Center declination in degrees (returned) */
    double dra;         /* Right ascension half-width in degrees (returned) */
    double ddec;        /* Declination half-width in degrees (returned) */
    double secpix;      /* Arcseconds per pixel (returned) */
    double crpix1, crpix2;
    double cxpix, cypix;
    double ra, dec;
    int offscl;
    int wp;             /* Image width in pixels (returned) */
    int hp;             /* Image height in pixels (returned) */
    int sysout=0;       /* Coordinate system to return (0=image, returned) */
    double eqout=0.0;   /* Equinox to return (0=image, returned) */
    int i, j;
    int ifrow1, ifrow2;	/* First and last rows to extract */
    int ifcol1, ifcol2;	/* First and last columns to extracT */
    int nbimage;        /* Number of bytes in image */
    int nbfile;         /* Number of bytes in file */
    int nbread;         /* Number of bytes actually read from file */
    char *ext;		/* Pointer to start of extension */
    int ikwd;
    char *endchar;
    char *ltime;
    char *newimage;
    char *kw, *kwl;
    int nc, ier, bitpix;
    int tpix, bytepix, nprow;
    int nbrow, nblock, nbleft, sampix, nbout;
    int nrows, ncols;
    char rastr[32], decstr[32], cstr[32];
    char temp[80];
    int vpix;
    int xdim = 1;
    int ydim = 1;
    int drow = 0;
    int dcol = 0;
    struct Range *xrange;    /* Column range structure */
    struct Range *yrange;    /* Row range structure */
    struct WorldCoor *wcs, *GetWCSFITS();

    /* Check to see if this is a FITS file */
    if (!isfits (name)) {
	fprintf (stderr, "File %s is not a FITS file\n", name);
	return (-1);
	}

    /* Read FITS header */
    if ((header = fitsrhead (name, &lhead, &nbhead)) == NULL) {
	fprintf (stderr, "Cannot read FITS file %s\n", name);
	return (-1);
	}

    /* If requested, delete keywords one at a time */
    for (ikwd = 0; ikwd < nkwd; ikwd++) {

	/* Make keyword all upper case */
	kwl = kwd[ikwd] + strlen (kwd[ikwd]);
	for (kw = kwd[ikwd]; kw < kwl; kw++) {
	    if (*kw > 96 && *kw < 123)
		*kw = *kw - 32;
	    }

	/* Delete keyword */
	if (hdel (header, kwd[ikwd]) && verbose)
	    printf ("%s: %s deleted\n", name, kwd[ikwd]);
	}


    if (verbose && first) {
	fprintf (stderr,"%s\n",RevMsg);
	fprintf (stderr, "Extract from FITS image file %s\n", name);
	}

    if (ra0 > -99.0 && dec0 > -99.0) {
	wcs = GetWCSFITS (name, header, verbose);
	if (iswcs (wcs)) {
	    ra = ra0;
	    dec = dec0;
	    wcscon (syscoor, wcs->syswcs, eqcoor, eqout, &ra, &dec, wcs->epoch);
	    wcs2pix (wcs, ra, dec, &cxpix, &cypix, &offscl);
	    xcpix = (int) (cxpix + 0.5);
	    ycpix = (int) (cypix + 0.5);
	    wcsfree (wcs);
	    }
	else {
	    fprintf (stderr, "No WCS in FITS image file %s\n", name);
	    wcsfree (wcs);
	    return (-1);
	    }
	}

    ncols = 0;
    (void) hgeti4 (header, "NAXIS1", &ncols);
    nrows = 1;
    (void) hgeti4 (header, "NAXIS2", &nrows);
    bitpix = 16;
    (void) hgeti4 (header, "BITPIX", &bitpix);
    if (bitpix < 0)
	bytepix = -bitpix / 8;
    else
	bytepix = bitpix / 8;

    /* Set up limiting rows to read */
    if (rrange) {
	yrange = RangeInit (rrange, ydim);
	ifrow1 = yrange->ranges[0];
	ifrow2 = yrange->ranges[1];
	}
    else if (xcpix && ycpix) {
	int ydpix2 = ydpix / 2;
	ifrow1 = ycpix - ydpix2 + 1;
	ifrow2 = ifrow1 + ydpix - 1;
	}
    else{
	ifrow1 = 1;
	ifrow2 = nrows;
	}
    if (ifrow1 < 1)
	ifrow1 = 1;
    if (ifrow2 < 1)
	ifrow2 = 1;
    if (ifrow1 > nrows)
	ifrow1 = nrows;
    if (ifrow2 > nrows)
	ifrow2 = nrows;
    hp = ifrow2 - ifrow1 + 1;

    /* Set up first and number of columns to write */
    if (crange) {
	xrange = RangeInit (crange, xdim);
	ifcol1 = xrange->ranges[0];
	ifcol2 = xrange->ranges[1];
	}
    else if (xdpix) {
	int xdpix2 = xdpix / 2;
	ifcol1 = xcpix - xdpix2 + 1;
	ifcol2 = ifcol1 + xdpix - 1;
	}
    else {
	ifcol1 = 1;
	ifcol2 = ncols;
	}
    if (ifcol1 < 1)
	ifcol1 = 1;
    if (ifcol2 < 1)
	ifcol2 = 1;
    if (ifcol1 > ncols)
	ifcol1 = ncols;
    if (ifcol2 > ncols)
	ifcol2 = ncols;
    wp = ifcol2 - ifcol1 + 1;

    /* Extract image */
    newimage = fitsrsect (name, header, nbhead, ifcol1, ifrow1, wp, hp, nlog);

    hputi4 (header, "NAXIS", 2);
    hputi4 (header, "NAXIS1", wp);
    hputi4 (header, "NAXIS2", hp);
    hdel (header, "NAXIS3");
    hdel (header, "NAXIS4");
    hdel (header, "NAXIS5");
    nbout = wp * bytepix;
    nbimage = nbout * hp;
    nblock = nbimage / 2880;
    nbleft = (nblock + 1) * 2880 - nbimage;
    if (nbleft > 2880)
	nbleft = 0;

    /* Reset image WCS if present */
    if (ifcol1 > 1 || ifrow1 > 1) {
	if (hgetr8 (header, "CRPIX1", &crpix1)) {
	    crpix1 = crpix1 - (double) (ifcol1 - 1);
	    hputr8 (header, "CRPIX1", crpix1);
	    }
	else {
	    hputs (header, "CTYPE1", "PIXEL");
	    hputi4 (header, "CRPIX1", 1);
	    hputi4 (header, "CRVAL1", ifcol1);
	    hputi4 (header, "CDELT1", 1);
	    }
	if (hgetr8 (header, "CRPIX1P", &crpix1)) {
	    crpix1 = crpix1 - (double) (ifcol1 - 1);
	    hputr8 (header, "CRPIX1P", crpix1);
	    }
	else {
	    hputs (header, "WCSNAMEP", "PLATE");
	    hputs (header, "CTYPE1P", "PIXEL");
	    hputi4 (header, "CRPIX1P", 1);
	    hputi4 (header, "CRVAL1P", ifcol1);
	    hputi4 (header, "CDELT1P", 1);
	    }
	if (hgetr8 (header, "CRPIX2", &crpix2)) {
	    crpix2 = crpix2 - (double) (ifrow1 - 1);
	    hputr8 (header, "CRPIX2", crpix2);
	    }
	else {
	    hputs (header, "CTYPE2", "PIXEL");
	    hputi4 (header, "CRPIX2", 1);
	    hputi4 (header, "CRVAL2", ifrow1);
	    hputi4 (header, "CDELT2", 1);
	    }
	if (hgetr8 (header, "CRPIX2P", &crpix2)) {
	    crpix2 = crpix2 - (double) (ifrow1 - 1);
	    hputr8 (header, "CRPIX2P", crpix2);
	    }
	else {
	    hputs (header, "CTYPE2P", "PIXEL");
	    hputi4 (header, "CRPIX2P", 1);
	    hputi4 (header, "CRVAL2P", ifrow1);
	    hputi4 (header, "CDELT2P", 1);
	    }
	}

    /* Add HISTORY notice of this conversion */
    if (ra0 > -99.0 && dec0 > -99.0) {
	strcpy (history, RevMsg);
	endchar = strchr (history, ',');
	*endchar = (char) 0;
	strcat (history, " ");
	ltime = lt2fd ();
	strcat (history, ltime);
	endchar = strrchr (history,':');
	*endchar = (char) 0;
	ra2str (rastr, 31, ra0, 3);
	dec2str (decstr, 31, dec0, 2);
	wcscstr (cstr, syscoor, eqcoor, 0.0);
	sprintf (temp, " %d x %d centered on %s %s %s",
		 xdpix, ydpix, rastr, decstr, cstr);
	strcat (history, temp);
	hputc (header, "HISTORY", history);
	}
    if (xcpix && ycpix) {
	strcpy (history, RevMsg);
	endchar = strchr (history, ',');
	*endchar = (char) 0;
	strcat (history, " ");
	ltime = lt2fd ();
	strcat (history, ltime);
	endchar = strrchr (history,':');
	*endchar = (char) 0;
	sprintf (temp, " %d x %d centered on [%d,%d]",
		 xdpix, ydpix, xcpix, ycpix);
	strcat (history, temp);
	hputc (header, "HISTORY", history);
	}
    else if (crange && rrange) {
	strcpy (history, RevMsg);
	endchar = strchr (history, ',');
	*endchar = (char) 0;
	strcat (history, " ");
	ltime = lt2fd ();
	strcat (history, ltime);
	endchar = strrchr (history,':');
	*endchar = (char) 0;
	sprintf (temp, " Rows %d - %d, Columns %d - %d",
		 ifrow1, ifrow2, ifcol1, ifcol2);
	strcat (history, temp);
	hputc (header, "HISTORY", history);
	}

    /* If assigned output name, use it */
    if (outname[0] != (char) 0)
	strcpy (fitsname, outname);

    /* Otherwise, create a new name */
    else
	nextname (name, fitsname);

    /* Write FITS file to a specified directory */
    if (outdir[0] > 0) {
	strcpy (fitspath, outdir);
	strcat (fitspath, "/");
	fitsfile = strrchr (fitsname,'/');
	if (fitsfile == NULL)
	    strcat (fitspath, fitsname);
	else
	    strcat (fitspath, fitsfile+1);
	}
    else
	strcpy (fitspath, fitsname);

    if (fitswimage (fitspath, header, newimage) > 0) {
	if (verbose)
	    fprintf (stderr, "%s: written successfully.\n", fitspath);
	else
	    printf ("%s\n", fitspath);
	}
    else if (verbose)
	fprintf (stderr, "NEWFITS: File %s not written.\n", fitspath);

    if (verbose)
	fprintf (stderr, "\n");
    free (newimage);
    free (header);
    return (0);
}


/* RANGEINIT -- Initialize range structure from string */

static struct Range *
RangeInit (string, ndef)

char	*string;	/* String containing numbers separated by , and - */
int	ndef;		/* Maximum allowable range value */

{
    struct Range *range;
    int ip, irange;
    char *slast;
    double first, last, step;

    if (!isrange (string) && !isnum (string))
	return (NULL);
    ip = 0;
    range = (struct Range *)calloc (1, sizeof (struct Range));
    range->irange = -1;
    range->nvalues = 0;
    range->nranges = 0;

    for (irange = 0; irange < MAXRANGE; irange++) {

	/* Default to entire list */
	first = 1.0;
	last = ndef;
	step = 1.0;

	/* Skip delimiters to start of range */
	while (string[ip] == ' ' || string[ip] == '	' ||
	       string[ip] == ',')
	    ip++;

	/* Get first limit
	 * Must be a number, '-', 'x', or EOS.  If not return ERR */
	if (string[ip] == (char)0) {	/* end of list */
	    if (irange == 0) {

		/* Null string defaults */
		range->ranges[0] = first;
		if (first < 1)
		    range->ranges[1] = first;
		else
		    range->ranges[1] = last;
		range->ranges[2] = step;
		range->nvalues = range->nvalues + 1 +
			  ((range->ranges[1]-range->ranges[0])/step);
		range->nranges++;
		return (range);
		}
	    else
		return (range);
	    }
	else if (string[ip] > (char)47 && string[ip] < 58) {
	    first = strtod (string+ip, &slast);
	    ip = slast - string;
	    }
	else if (strchr ("-:x", string[ip]) == NULL) {
	    free (range);
	    return (NULL);
	    }

	/* Skip delimiters */
	while (string[ip] == ' ' || string[ip] == '	' ||
	       string[ip] == ',')
	    ip++;

	/* Get last limit
	* Must be '-', or 'x' otherwise last = first */
	if (string[ip] == '-' || string[ip] == ':') {
	    ip++;
	    while (string[ip] == ' ' || string[ip] == '	' ||
	   	   string[ip] == ',')
		ip++;
	    if (string[ip] == (char)0)
		last = first + ndef;
	    else if (string[ip] > (char)47 && string[ip] < 58) {
		last = strtod (string+ip, &slast);
		ip = slast - string;
		}
	    else if (string[ip] != 'x')
		last = first + ndef;
	    }
	else if (string[ip] != 'x')
	    last = first;

	/* Skip delimiters */
	while (string[ip] == ' ' || string[ip] == '	' ||
	       string[ip] == ',')
	    ip++;

	/* Get step
	 * Must be 'x' or assume default step. */
	if (string[ip] == 'x') {
	    ip++;
	    while (string[ip] == ' ' || string[ip] == '	' ||
	   	   string[ip] == ',')
		ip++;
	    if (string[ip] == (char)0)
		step = 1.0;
	    else if (string[ip] > (char)47 && string[ip] < 58) {
		step = strtod (string+ip, &slast);
		ip = slast - string;
		}
	    else if (string[ip] != '-' && string[ip] != ':')
		step = 1.0;
            }

	/* Output the range triple */
	range->ranges[irange*3] = first;
	range->ranges[irange*3 + 1] = last;
	range->ranges[irange*3 + 2] = step;
	range->nvalues = range->nvalues + ((last-first+(0.1*step)) / step + 1);
	range->nranges++;
	}

    return (range);
}


/* ISRANGE -- Return 1 if string is a range, else 0 */

static int
isrange (string)

char *string;		/* String which might be a range of numbers */

{
    int i, lstr;

    /* If range separators present, check to make sure string is range */
    if (strchr (string+1, '-') || strchr (string+1, ',')) {
	lstr = strlen (string);
	for (i = 0; i < lstr; i++) {
	    if (strchr ("0123456789-,.x", (int)string[i]) == NULL)
		return (0);
	    }
	return (1);
	}
    else
	return (0);
}


/* RSTART -- Restart at beginning of range */

static void
rstart (range)

struct Range *range;	/* Range structure */

{
    range->irange = -1;
    return;
}


/* RGETN -- Return number of values from range structure */

static int
rgetn (range)

struct Range *range;	/* Range structure */

{
    return (range->nvalues);
}


/*  RGETR8 -- Return next number from range structure as 8-byte f.p. number */

static double
rgetr8 (range)

struct Range *range;	/* Range structure */

{
    int i;

    if (range == NULL)
	return (0.0);
    else if (range->irange < 0) {
	range->irange = 0;
	range->first = range->ranges[0];
	range->last = range->ranges[1];
	range->step = range->ranges[2];
	range->value = range->first;
	}
    else {
	range->value = range->value + range->step;
	if (range->value > range->last) {
	    range->irange++;
	    if (range->irange < range->nranges) {
		i = range->irange * 3;
		range->first = range->ranges[i];
		range->last = range->ranges[i+1];
		range->step = range->ranges[i+2];
		range->value = range->first;
		}
	    else
		range->value = 0.0;
	    }
	}
    return (range->value);
}


/*  RGETI4 -- Return next number from range structure as 4-byte integer */

static int
rgeti4 (range)

struct Range *range;	/* Range structure */

{
    double value;

    value = rgetr8 (range);
    return ((int) (value + 0.000000001));
}

static void
nextname (name, newname)

char *name;
char *newname;
{
    char *ext, *sufchar;
    int lext, lname;

    ext = strrchr (name, '.');
    if (ext)
	lname = ext - name;
    else
	lname = strlen (name);
    strncpy (newname, name, lname);
    sufchar = newname + lname;
    *sufchar = 'a';
    *(sufchar+1) = (char) 0;
    strcat (newname, ext);
    while (!access (newname, F_OK)) {
	if (*sufchar < 'z')
	    (*sufchar)++;
	else {
	    *sufchar = 'a';
	    sufchar++;
	    *sufchar = 'a';
	    *(sufchar+1) = (char) 0;
	    strcat (newname, ext);
	    }
	}
    return;
}

/* Oct 22 2002	New program based on t2f
 * Dec  6 2002	Initialize bytepix, which wasn't
 * Dec 16 2002	Add -k option to delete FITS keywords when copying
 *
 * Jan 30 2003	Fix typo in variable name 
 * May  2 2003	Fix bug if no keywords are deleted
 *
 * Apr 16 2004	Delete NAXISn for n > 2 in output image
 * Sep 15 2004	Fix bug dealing with center specified as sky coordinates
 * Sep 17 2004	Add option to set extraction center as decimal degrees
 * Dec  6 2004	Don't print gratuitous newline at end of process
 *
 * Sep 30 2005	Convert input center coordinates to image system
 */
