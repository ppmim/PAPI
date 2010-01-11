/* File sethead.c
 * March 7, 2005
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

#define MAXKWD 100
#define MAXFILES 1000
static int maxnkwd = MAXKWD;
static int maxnfile = MAXFILES;
#define MAXNEW 1024

#define KEY_ADD 1
#define KEY_SUB 2
#define KEY_MUL 3
#define KEY_DIV 4

static void usage();
static void SetValues ();

static int verbose = 0;		/* verbose/debugging flag */
static int newimage0 = 0;
static int keyset = 0;
static int histset = 0;
static int krename = 0;
static char prefix[2];
static int version = 0;		/* If 1, print only program name and version */
static int logfile = 0;
static int first_file = 1;
static char spchar = (char) 0;	/* Character to replace with spaces */
static int nproc = 0;

static char *RevMsg = "SETHEAD WCSTools 3.6.4, 3 May 2006, Doug Mink (dmink@cfa.harvard.edu)";

int
main (ac, av)
int ac;
char **av;
{
    char *str;
    char **kwd, **kwdnew;
    char **comment, **comnew;
    int nkwd = 0;
    int nkwd1 = 0;
    char **fn;
    int nfile = 0;
    int readlist = 0;
    int ifile;
    int icom, lcom;
    char filename[1024];
    char *keybuff, *kw1, *kw2, *kwdi, *kwdc;
    FILE *flist;
    char *listfile;
    char *ilistfile;
    char *klistfile;
    int ikwd, i, nc;
    char *dq, *sq, *sl;
    char lf = (char) 10;
    char cr = (char) 13;
    char dquote = (char) 34;
    char squote = (char) 39;

    ilistfile = NULL;
    klistfile = NULL;
    keybuff = NULL;
    nkwd = 0;
    nfile = 0;
    fn = (char **)calloc (maxnfile, sizeof(char *));
    for (ifile = 0; ifile < maxnfile; ifile++)
	fn[ifile] = NULL;
    kwd = (char **)calloc (maxnkwd, sizeof(char *));
    comment = (char **)calloc (maxnkwd, sizeof(char *));
    for (ikwd = 0; ikwd < maxnkwd; ikwd++) {
	kwd[ikwd] = NULL;
	comment[ikwd] = NULL;
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

		case 'h':	/* Set HISTORY */
		    histset++;
		    break;
	
		case 'k':	/* Set SETHEAD keyword */
		    keyset++;
		    break;
	
		case 'l':	/* Log files changed */
		    logfile++;
		    break;

		case 'm':	/* Maximum number of keywords to be changed */
		    if (ac > 1) {
			maxnkwd = atoi (++av[0]);
			ac--;
			free (kwd);
			kwd = (char **) calloc (maxnkwd, sizeof (void *));
			for (ikwd = 0; ikwd < maxnkwd; ikwd++)
			    kwd[ikwd] = NULL;
			kwd = kwdnew;
			free (comment);
			comment = (char **) calloc (maxnkwd, sizeof (void *));
			for (ikwd = 0; ikwd < maxnkwd; ikwd++)
			    comment[ikwd] = NULL;
			}
		    break;

		case 'n':	/* Write new file */
		    newimage0++;
		    break;

		case 'r':	/* Rename keywords with replaced values */
		    krename++;
		    if (ac > 1) {
			strncpy (prefix, *(++av), 1);
			ac--;
			}
		    else
			prefix[0] = 'X';
		    prefix[1] = (char) 0;
		    break;

		case 's':	/* Replace this character with spaces in string arguments */
		    if (ac > 1) {
			strncpy (&spchar, *(++av), 1);
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

	/* File containing a list of keywords or files */
	else if (*av[0] == '@') {
	    readlist++;
	    listfile = *av + 1;
	    if (nfile == 0 && isimlist (listfile)) {
		ilistfile = listfile;
		nfile = getfilelines (ilistfile);
		}
	    else {
		klistfile = listfile;
		nkwd1 = getfilelines (klistfile);
		if (nkwd1 > 0) {
		    if (nkwd+nkwd1 >= maxnkwd) {
			maxnkwd = nkwd + nkwd1 + 32;
			kwdnew = (char **) calloc (maxnkwd, sizeof (void *));
			for (ikwd = 0; ikwd < nkwd; ikwd++)
			    kwdnew[ikwd] = kwd[ikwd];
			free (kwd);
			kwd = kwdnew;
			comnew = (char **) calloc (maxnkwd, sizeof (void *));
			for (ikwd = 0; ikwd < nkwd; ikwd++)
			    comnew[ikwd] = comment[ikwd];
			free (comment);
			comment = comnew;
			}
		    keybuff = getfilebuff (klistfile);
		    if (keybuff != NULL) {
			kw1 = keybuff;

			/* One keyword per line of buffer */
			for (ikwd = 0; ikwd < nkwd1; ikwd++) {
			    kwd[nkwd] = kw1;
			    kwdi = kw1;

			    /* Replace LF and/or CR with NULL */
			    if (ikwd < nkwd1 - 1) {
				kw2 = strchr (kwdi, lf);
				if (kw2 != NULL) {
				    kw1 = kw2 + 1;
				    *kw2 = (char) 0;
				    }
				kw2 = strchr (kwdi, cr);
				if (kw2 != NULL) {
				    kw1 = kw2 + 1;
				    *kw2 = (char) 0;
				    if (kw1 == (char) 0)
					kw1 = kw1 + 1;
				    }
				}

			    /* Replace final LF and/or CR with NULL */
			    else if (ikwd < nkwd1) {
				kw2 = strchr (kwdi, cr);
				if (kw2 != NULL)
				    *kw2 = (char) 0;
				kw2 = strchr (kwdi, lf);
				if (kw2 != NULL)
				    *kw2 = (char) 0;
				}

			    /* Check for presence of a comment past quotes */
			    sq = strchr (kwdi, squote);
			    dq = strchr (kwdi, dquote);
			    kwdc = kwdi;
			    if (sq != NULL && (dq > sq || dq == NULL)) {
				kwdc = strchr (sq+1, squote);
				if (kwdc != NULL)
				    kwdc = kwdc + 1;
				}
			    else if (dq != NULL && (sq > dq || sq == NULL)) {
				kwdc = strchr (dq+1, dquote);
				if (kwdc != NULL)
				    kwdc = kwdc + 1;
				}
			    if (kwdc != NULL) {
				sl = strsrch (kwdc, " / ");
				if (sl != NULL) {
				    *sl = (char) 0;
				    comment[nkwd] = sl + 3;
				    if (spchar)
					stc2s (spchar, comment[nkwd]);
				    }
				}
			    nkwd++;
			    }
			}
		    }
		}
	    }

	/* Comment for preceding keyword */
	else if (*av[0] == '/' && strlen (*av) == 1) {
	    av++;
	    ac--;
	    if (nkwd > 0) {
		if (spchar)
		    stc2s (spchar, *av);
		comment[nkwd-1] = *av;
		}
	    }

	/* Image file */
	else if (isfits (*av) || isiraf (*av)) {
	    if (nfile >= maxnfile) {
		maxnfile = maxnfile * 2;
		fn = (char **) realloc ((void *)fn, maxnfile);
		}
	    fn[nfile] = *av;
	    nfile++;
	    }

	/* Keyword */
	else {
	    if (nkwd >= maxnkwd) {
		maxnkwd = maxnkwd * 2;
		kwdnew = (char **) calloc (maxnkwd, sizeof (void *));
		for (ikwd = 0; ikwd < nkwd; ikwd++)
		    kwdnew[ikwd] = kwd[ikwd];
		free (kwd);
		kwd = kwdnew;
		comnew = (char **) calloc (maxnkwd, sizeof (void *));
		for (ikwd = 0; ikwd < nkwd; ikwd++)
		    comnew[ikwd] = comment[ikwd];
		free (comment);
		comment = comnew;
		}
	    kwd[nkwd] = *av;
	    nkwd++;
	    }
	}

    if (nkwd <= 0 && nfile <= 0 )
	usage ();
    else if (nkwd <= 0) {
	fprintf (stderr, "SETHEAD: no keywords specified\n");
	exit (1);
	}
    else if (nfile <= 0 ) {
	fprintf (stderr, "SETHEAD: no files specified\n");
	exit (1);
	}

    /* Open file containing a list of images, if there is one */
    if (ilistfile != NULL) {
	if ((flist = fopen (ilistfile, "r")) == NULL) {
	    fprintf (stderr,"SETHEAD: Image list file %s cannot be read\n",
		     ilistfile);
	    usage ();
	    }
	}

    /* Read through headers of images */
    for (ifile = 0; ifile < nfile; ifile++) {
	if (ilistfile != NULL) {
	    first_token (flist, 254, filename);
	    SetValues (filename, nkwd, kwd, comment);
	    }
	else
	    SetValues (fn[ifile], nkwd, kwd, comment);
	}
    if (ilistfile != NULL)
	fclose (flist);

    if (keybuff != NULL)
	free (keybuff);

    return (0);
}

static void
usage ()
{
    fprintf (stderr,"%s\n",RevMsg);
    if (version)
	exit (-1);
    fprintf (stderr,"Set FITS or IRAF header keyword values\n");
    fprintf(stderr,"Usage: [-hknv][-r [char]][-s char] file1.fits [... filen.fits] kw1=val1 [ / comment] [ ... kwn=valuen]\n");
    fprintf(stderr,"  or : [-hknv][-r [char]][-s char] file1.fits [... filen.fits] @keywordfile]\n");
    fprintf(stderr,"  or : [-hknv][-r [char]][-s char] @listfile kw1=val1 [ ... kwn=valuen]\n");
    fprintf(stderr,"  or : [-hknv][-r [char]][-s char] @listfile @keywordfile\n");
    fprintf(stderr,"  -h: Write HISTORY line\n");
    fprintf(stderr,"  -k: Write SETHEAD keyword\n");
    fprintf(stderr,"  -l: Log files as processed (on one line)\n");
    fprintf(stderr,"  -m num: Change max number of keywords changed to num\n");
    fprintf(stderr,"  -n: Write a new file (add e before the extension)\n");
    fprintf(stderr,"  -r [char]: Rename reset keywords with char or X prefixed\n");
    fprintf(stderr,"  -s [char]: Replace this character with space in string values\n");
    fprintf(stderr,"  -v: Verbose\n");
    fprintf(stderr,"  / comment: Add this comment to previous keyword\n");
    exit (1);
}


static void
SetValues (filename, nkwd, kwd, comment)

char	*filename;	/* Name of FITS or IRAF image file */
int	nkwd;		/* Number of keywords for which to set values */
char	*kwd[];		/* Names and values of those keywords */
char	*comment[];	/* Comments for those keywords (none if NULL) */

{
    char *header;	/* FITS image header */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    char *irafheader;	/* IRAF image header */
    int iraffile;	/* 1 if IRAF image, 0 if FITS image */
    int newimage;	/* 1 to write new image file, else 0 */
    int i, lext, lroot, ndec, isra;
    char *image;
    char newname[MAXNEW];
    char *newval;
    char string[80];
    char *ext, *fname, *imext, *imext1;
    char *kw, *kwv, *kwl, *kwv0, *knl;
    char *v, *vq0, *vq1;
    char echar;
    int ikwd, lkwd, lkwv, lhist;
    int fdr, fdw, ipos, nbr, nbw;
    int squote = 39;
    int dquote = 34;
    int naxis = 0;
    int nbold, nbnew;
    int imageread = 0;
    char cval[24];
    char history[72];
    char *endchar;
    char *ltime;
    char newkey[10];
    char value[80];
    char keyroot[8];
    char ctemp;
    int lval, ii, lv, lnl;
    int ival0, ival1;
    double dval0, dval1;
    int bitpix = 0;
    char newline = 10;
    int keyop = 0;
    char *opkey;
    char ops[8];

    strcpy (ops, "=+-*/");
    newimage = newimage0;

    /* Open IRAF image if .imh extension is present */
    if (isiraf (filename)) {
	iraffile = 1;
	if ((irafheader = irafrhead (filename, &lhead)) != NULL) {
	    if ((header = iraf2fits (filename, irafheader, lhead, &nbhead)) == NULL) {
		fprintf (stderr, "Cannot translate IRAF header %s/n",filename);
		free (irafheader);
		irafheader = NULL;
		return;
		}
	    }
	else {
	    fprintf (stderr, "Cannot read IRAF file %s\n", filename);
	    return;
	    }
	}

    /* Open FITS file if .imh extension is not present */
    else {
	iraffile = 0;
	setfitsinherit (0);
	if ((header = fitsrhead (filename, &lhead, &nbhead)) != NULL) {
	    hgeti4 (header,"NAXIS",&naxis);
	    hgeti4 (header,"BITPIX",&bitpix);
	    if (naxis > 0 && bitpix != 0) {
		if ((image = fitsrfull (filename, nbhead, header)) == NULL) {
		    if (verbose)
			fprintf (stderr, "No FITS image in %s\n", filename);
		    imageread = 0;
		    }
		else
		    imageread = 1;
		}
	    }
	else {
	    fprintf (stderr, "Cannot read FITS file %s\n", filename);
	    return;
	    }
	}
    if (verbose && first_file) {
	fprintf (stderr,"%s\n",RevMsg);
	fprintf (stderr,"Set Header Parameter Values in ");
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s\n", filename);
	else
	    fprintf (stderr,"FITS image file %s\n", filename);
	first_file = 0;
	}

    if (nkwd < 1)
	return;

    nbold = fitsheadsize (header);

    /* Set keywords one at a time */
    for (ikwd = 0; ikwd < nkwd; ikwd++) {
        strcpy (cval,"                    ");
	lkwd = strlen (kwd[ikwd]);
	if (!strncmp (kwd[ikwd], "COMMENT", 7) || !strncmp (kwd[ikwd], "HISTORY", 7) ||
	    !strncmp (kwd[ikwd], "comment", 7) || !strncmp (kwd[ikwd], "history", 7)) {
	    if (lkwd > 8) {
		kwv0 = kwd[ikwd] + 7;
		if (spchar)
		    stc2s (spchar, kwv0);
		}
	    else
		kwv0 = NULL;
	    }
	else {
	    kwv0 = strchr (kwd[ikwd], '=');
	    keyop = 0;
	    if (kwv0 == NULL) {
		kwv0 = strchr (kwd[ikwd], '+');
		if (kwv0 != NULL)
		    keyop = KEY_ADD;
		else {
		    kwv0 = strchr (kwd[ikwd], '/');
		    if (kwv0 != NULL)
			keyop = KEY_DIV;
		    else {
			kwv0 = strchr (kwd[ikwd], '*');
			if (kwv0 != NULL)
			    keyop = KEY_MUL;
			else {
			    kwv0 = strchr (kwd[ikwd], '-');
			    if (kwv0 != NULL)
				keyop = KEY_SUB;
			    }
			}
		    }
		}
	    }
	if (kwv0 == NULL) {

	    /* Get current length of header buffer */
	    lhead = gethlength (header);

	    /* Make keyword all upper case */
	    kwl = kwd[ikwd] + lkwd;
	    for (kw = kwd[ikwd]; kw < kwl; kw++) {
		if (*kw > 96 && *kw < 123)
		    *kw = *kw - 32;
		}
	    if (!strcmp (kwd[ikwd], "COMMENT") || !strcmp (kwd[ikwd], "HISTORY"))
		hputc (header, kwd[ikwd], "");
	    else if (hgets (header, kwd[ikwd], 80, value)) {
		if (isdate (value)) {
		    newval = fd2fd (value);
		    hputs (header, kwd[ikwd], newval);
		    }
		}
	    if (comment[ikwd]) {
		hputcom (header, kwd[ikwd], comment[ikwd]);
		if (verbose) {
		    hgets (header, kwd[ikwd], 80, value);
		    printf ("%s = %s / %s\n", kwd[ikwd], value, comment[ikwd]);
		    }
		}
	    }
	else {
	    *kwv0 = 0;
	    lkwd = kwv0 - kwd[ikwd];
	    kwv = kwv0 + 1;
	    lkwv = strlen (kwv);

	    /* If end of line before end of string, terminate value there */
	    knl = strchr (kwv,newline);
	    if (knl != NULL) {
		lnl = knl - kwv;
		if (lnl < lkwv) lkwv = lnl;
		}

	    /* Get current length of header buffer */
	    lhead = gethlength (header);

	    /* Make keyword all upper case */
	    kwl = kwd[ikwd] + lkwd;
	    for (kw = kwd[ikwd]; kw < kwl; kw++) {
		if (*kw > 96 && *kw < 123)
		    *kw = *kw - 32;
		}

	    /* If keyword is already in header, rename it if requested */
	    if (krename && ksearch (header, kwd[ikwd])) {
		strcpy (newkey, prefix);
		strcat (newkey, kwd[ikwd]);
		if (strlen (newkey) > 8)
		    newkey[8] = (char) 0;
		hchange (header, kwd[ikwd], newkey);
		opkey = newkey;
		}
	    else
		opkey = kwd[ikwd];

	    /* Add, subtract, multiply, or divide keyword value by constant */
	    if (keyop) {
		if (!hgets (header, opkey, 72,string)) {
		    if (verbose)
			printf ("* %s %c %s keyword not in header.\n",
			    opkey, ops[keyop], kwv);
		    continue;
		    }
		if (!isnum (string) && !strchr (string,':')) {
		    if (verbose) 
			printf ("* %s = %s in header is not a number.\n",
			    opkey, ops[keyop], string);
		    continue;
		    }
		if (!isnum (kwv) && !strchr (kwv, ':')) {
		    if (verbose) 
			printf ("* %s %c %s not a number.\n",
			    opkey, ops[keyop], kwv);
		    continue;
		    }

		/* Make sexagesimal output */
		if (strchr (string, ':')) {
		    if (strsrch (opkey, "RA") || strsrch (opkey, "HA")) {
			isra = 1;
			dval0 = str2ra (string);
			}
		    else {
			isra = 0;
			dval0 = str2dec (string);
			}
		    if (strchr (kwv, ':')) {
			if (isra)
			    dval1 = str2ra (kwv);
			else
			    dval1 = str2dec (kwv);
			}
		    else
			dval1 = atof (kwv);
		    
		    if (keyop == KEY_ADD)
			dval0 = dval0 + dval1;
		    else if (keyop == KEY_SUB)
			dval0 = dval0 - dval1;
		    else if (keyop == KEY_MUL)
			dval0 = dval0 * dval1;
		    else if (keyop == KEY_DIV && dval1 != 0)
			dval0 = dval0 / dval1;
		    ndec = numdec (string);
		    if (isra)
			ra2str (string, 32, dval0, ndec);
		    else
			dec2str (string, 32, dval0, ndec);
		    hputs (header, opkey, string);
		    }

		/* Make integer output */
		else if (isnum (string) == 1) {
		    if (kwv[0] == '-')
			ival1 = (int) (atof (kwv) - 0.5);
		    else
			ival1 = (int) (atof (kwv) + 0.5);
		    ival0 = atoi (string);
		    if (keyop == KEY_ADD)
			ival0 = ival0 + ival1;
		    else if (keyop == KEY_SUB)
			ival0 = ival0 - ival1;
		    else if (keyop == KEY_MUL)
			ival0 = ival0 * ival1;
		    else if (keyop == KEY_DIV && ival1 != 0)
			ival0 = ival0 / ival1;
		    hputi4 (header, kwd[ikwd], ival0);
		    }

		/* Make floating point output */
		else {
		    dval0 = atof (string);
		    ndec = numdec (string);
		    dval1 = atof (kwv);
		    if (keyop == KEY_ADD)
			dval0 = dval0 + dval1;
		    else if (keyop == KEY_SUB)
			dval0 = dval0 - dval1;
		    else if (keyop == KEY_MUL)
			dval0 = dval0 * dval1;
		    else if (keyop == KEY_DIV && dval1 != 0)
			dval0 = dval0 / dval1;
		    hputnr8 (header, kwd[ikwd], ndec, dval0);
		    }
		continue;
		}

	    /* Write comment without quotes, filling spaces first */
	    if (!strcmp (kwd[ikwd], "COMMENT") || !strcmp (kwd[ikwd], "HISTORY")) {
		if (spchar)
		    stc2s (spchar, kwv);
		hputc (header, kwd[ikwd], kwv);
		}

	    /* Write numeric value to keyword */
	    else if (isnum (kwv)) {
		i = 21 - lkwv;
		for (v = kwv; v < kwv+lkwv; v++)
		    cval[i++] = *v;
		cval[21] = 0;
		if (hputc (header, kwd[ikwd], cval)) {
		    lhead = lhead + 28800;
		    if ((header =
			(char *)realloc(header,(unsigned int)lhead)) != NULL) {
			hlength (header, lhead);
			hputc (header,kwd[ikwd], cval);
			}
		    }
		}

	    /* Write boolean value to keyword */
	    else if (!strcmp (kwv,"T") || !strcmp (kwv,"t") ||
		     !strcmp (kwv,"YES") || !strcmp (kwv,"yes")) {
		if (hputl (header, kwd[ikwd], 1)) {
		    lhead = lhead + 28800;
		    if ((header =
			(char *)realloc(header,(unsigned int)lhead)) != NULL) {
			hlength (header, lhead);
			hputl (header,kwd[ikwd], 1);
			}
		    }
		}
	    else if (!strcmp (kwv,"F") || !strcmp (kwv,"f") ||
		     !strcmp (kwv,"NO") || !strcmp (kwv,"no")) {
		if (hputl (header, kwd[ikwd], 0)) {
		    lhead = lhead + 28800;
		    if ((header =
			(char *)realloc(header,(unsigned int)lhead)) != NULL) {
			hlength (header, lhead);
			hputl (header,kwd[ikwd], 0);
			}
		    }
		}

	    /* Write character string to keyword */
	    else {

		/* Remove double quotes from character string */
		if ((vq0 = strchr (kwv,dquote)) != NULL) {
		    vq0 = vq0 + 1;
		    vq1 = strchr (vq0,dquote);
		    if (vq0 != NULL && vq1 != NULL) {
			kwv = vq0;
			*vq1 = (char) 0;
			}
		    }

		/* Remove single quotes from character string */
		else if ((vq0 = strchr (kwv,squote)) != NULL) {
		    vq0 = vq0 + 1;
		    vq1 = strchr (vq0,squote);
		    if (vq0 != NULL && vq1 != NULL) {
			kwv = vq0;
			*vq1 = (char) 0;
			}
		    }

		/* Replace special characters with spaces if not quoted */
		else if (spchar)
		    stc2s (spchar, kwv);

		/* Remove trailing blanks */
		lval = strlen (kwv);
		while (kwv[lval-1] < (char) 33)
		    kwv[--lval] = (char) 0;

		if (lval < 69) {
		    if (hputs (header, kwd[ikwd], kwv)) {
			lhead = lhead + 14400;
		    	if ((header =
			(char *)realloc(header,(unsigned int)lhead)) != NULL) {
			    hlength (header, lhead);
			    hputs (header, kwd[ikwd], kwv);
			    }
			}
		    }

		/* If character string is longer than 68 characters, split it */
		else {
		    strcpy (keyroot, kwd[ikwd]);
		    lroot = strlen (keyroot);
		    if (lroot > 6) {
			*(keyroot+6) = (char) 0;
			lroot = 6;
			}
		    ii = '1';
		    lkwv = strlen (kwv);
		    while (lkwv > 0) {
			if (lkwv > 67)
			    lv = 67;
			else
			    lv = lkwv;
			strncpy (value, kwv, lv);
			ctemp = value[lv];
			value[lv] = (char) 0;
			strcpy (newkey, keyroot);
			strcat (newkey, "_");
			newkey[lroot+1] = ii;
			newkey[lroot+2] = (char) 0;
			ii++;
			if (hputs (header, newkey, value)) {
			    lhead = lhead + 28800;
			    if ((header =
			  (char *)realloc(header,(unsigned int)lhead))!=NULL) {
				hlength (header, lhead);
				hputs (header, newkey, value);
				}
			    }
			value[lv] = ctemp;
			kwv = kwv + lv;
			lkwv = lkwv - lv;
			}
		    kwv = kwv0 + 1;
		    }
		}
	    if (comment[ikwd]) {
		hputcom (header, kwd[ikwd], comment[ikwd]);
		if (verbose)
		    printf ("%s = %s / %s\n", kwd[ikwd], kwv, comment[ikwd]);
		}
	    else if (verbose)
		printf ("%s = %s\n", kwd[ikwd], kwv);
	    *kwv0 = '=';
	    }
	}

    /* Add history to header */
    if (keyset || histset) {
	if (hgets (header, "SETHEAD", 72, history))
	    hputc (header, "HISTORY", history);
	strcpy (history, RevMsg);
	endchar = strchr (history, ',');
	*endchar = (char) 0;
	strcat (history, " ");
	ltime = lt2fd ();
	strcat (history, ltime);
	endchar = strrchr (history,':');
	*endchar = (char) 0;
	strcat (history, " ");
	for (ikwd = 0; ikwd < nkwd; ikwd++) {
	    kwv0 = strchr (kwd[ikwd],'=');
	    if (kwv0 == NULL)
	    keyop = 0;
	    if (kwv0 == NULL) {
		kwv0 = strchr (kwd[ikwd], '+');
		if (kwv != NULL)
		    keyop = KEY_ADD;
		else {
		    kwv0 = strchr (kwd[ikwd], '/');
		    if (kwv != NULL)
			keyop = KEY_DIV;
		    else {
			kwv0 = strchr (kwd[ikwd], '*');
			if (kwv != NULL)
			    keyop = KEY_MUL;
			else {
			    kwv0 = strchr (kwd[ikwd], '-');
			    if (kwv != NULL)
				keyop = KEY_SUB;
			    }
			}
		    }
		}
	    if (kwv0)
		*kwv0 = (char) 0;
	    lhist = strlen (history);
	    lkwd = strlen (kwd[ikwd]);

	    /* If too may keywords, start a second history line */
	    if (lhist + lkwd + 10 > 71) {
		if (histset) {
		    strcat (history, " updated");
		    hputc (header, "HISTORY", history);
		    strcpy (history, RevMsg);
		    endchar = strchr (history, ',');
		    *endchar = (char) 0;
		    strcat (history, " ");
		    ltime = lt2fd ();
		    strcat (history, ltime);
		    endchar = strrchr (history,':');
		    *endchar = (char) 0;
		    strcat (history, " ");
		    }
		else
		    break;
		}
	    strcat (history, kwd[ikwd]);
	    if (kwv0)
		*kwv0 = ops[keyop];
	    if (nkwd == 2 && ikwd < nkwd-1)
		strcat (history, " and ");
	    else if (ikwd < nkwd-1)
		strcat (history, ", ");
	    }
	strcat (history, " updated");
	if (keyset)
	    hputs (header, "SETHEAD", history);
	if (histset)
	    hputc (header, "HISTORY", history);
	}

    /* Compare size of output header to size of input header */
    nbnew = fitsheadsize (header);
    if (nbnew > nbold  && naxis == 0 && bitpix != 0) {
	if (verbose)
	    fprintf (stderr, "Rewriting primary header, copying rest of file\n");
	newimage = 1;
	}

    /* Make up name for new FITS or IRAF output file */
    if (newimage) {

    /* Remove directory path and extension from file name */
	fname = strrchr (filename, '/');
	if (fname)
	    fname = fname + 1;
	else
	    fname = filename;
	ext = strrchr (fname, '.');
	if (ext != NULL) {
	    lext = (fname + strlen (fname)) - ext;
	    lroot = ext - fname;
	    if (lroot > MAXNEW)
		lroot = MAXNEW - 1;
	    strncpy (newname, fname, lroot);
	    newname[lroot] = (char) 0;
	    }
	else {
	    lext = 0;
	    lroot = strlen (fname);
	    if (lroot > MAXNEW)
		lroot = MAXNEW - 1;
	    strncpy (newname, fname, lroot);
	    newname[lroot] = (char) 0;
	    }
	imext = strchr (fname, ',');
	imext1 = NULL;
	if (imext == NULL) {
	    imext = strchr (fname, '[');
	    if (imext != NULL) {
		imext1 = strchr (fname, ']');
		*imext1 = (char) 0;
		}
	    }
	if (imext != NULL && *(imext+1) != '0') {
	    strcat (newname, "_");
	    strcat (newname, imext+1);
	    }
	if (fname)
	    fname = fname + 1;
	else
	    fname = filename;
	strcat (newname, "e");
	if (lext > 0) {
	    if (imext != NULL) {
		echar = *imext;
		*imext = (char) 0;
		strcat (newname, ext);
		*imext = echar;
		if (imext1 != NULL)
		    *imext1 = ']';
		}
	    else
		strcat (newname, ext);
	    }
	}
    else {
	strcpy (newname, filename);
	if (!imext && ksearch (header,"XTENSION")) {
	    strcat (newname, ",1");
	    imext = strchr (newname,',');
	    }
	}

    /* Write fixed header to output file */
    if (iraffile) {
	if (irafwhead (newname, lhead, irafheader, header) > 0 && verbose)
	    printf ("%s rewritten successfully.\n", newname);
	else if (verbose)
	    printf ("%s could not be written.\n", newname);
	free (irafheader);
	irafheader = NULL;
	}

    /* If there is no data, write header by itself */
    else if (bitpix == 0) {
	if ((fdw = fitswhead (newname, header)) > 0) {
	    if (verbose)
		printf ("%s: rewritten successfully.\n", newname);
	    close (fdw);
	    }
	}

    /* Rewrite only header if it fits into the space from which it was read */
    else if (nbnew <= nbold && !newimage) {
	if (!fitswexhead (newname, header)) {
	    if (verbose)
		printf ("%s: rewritten successfully.\n", newname);
	    }
	}

    /* Rewrite header and data to a new image file */
    else if (naxis > 0 && imageread) {
	if (fitswimage (newname, header, image) > 0 && verbose)
	    printf ("%s: rewritten successfully.\n", newname);
	else if (verbose)
	    printf ("%s could not be written.\n", newname);
	free (image);
	image = NULL;
	}

    else {
	if ((fdw = fitswhead (newname, header)) > 0) {
	    fdr = fitsropen (filename);
	    ipos = lseek (fdr, nbhead, SEEK_SET);
	    image = (char *) calloc (2880, 1);
	    while ((nbr = read (fdr, image, 2880)) > 0) {
		nbw = write (fdw, image, nbr);
		if (nbw < nbr)
		    fprintf (stderr,"SETHEAD: %d / %d bytes written\n",nbw,nbr);
		}
	    close (fdr);
	    close (fdw);
	    if (verbose)
		printf ("%s: rewritten successfully.\n", newname);
	    }
	}

    /* Log the processing of this file, if requested */
    if (logfile) {
	nproc++;
	fprintf (stderr, "%d: %s processed.\r", nproc, newname);
	}

    free (header);
    header = NULL;
    if (image != NULL) {
	free (image);
	image = NULL;
	}
    return;
}

/* Oct 11 1996	New program
 * Dec 12 1996	Move ISNUM subroutine to hget.c
 *
 * Feb 21 1997  Check pointers against NULL explicitly for Linux
 *
 * May 28 1998	Include fitsio.h instead of fitshead.h
 * Jun  2 1998	Fix bug in hput()
 * Jul 24 1998	Deal coorectly with logical T or F
 * Jul 24 1998	Make irafheader char instead of int
 * Jul 30 1998	Allow use of list of files and multiple files
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 * Aug 14 1998	Preserve extension when creating new file name
 * Aug 14 1998	If changing primary header, write out entire input file
 * Aug 31 1998	Add options to add HISTORY and/or SETHEAD keyword
 * Sep  1 1998	Add option to keep changed keywords with new names
 * Oct  5 1998	Allow header changes even if no data is present
 * Oct  5 1998	Determine assignment arguments by presence of equal sign
 * Oct  5 1998	Use isiraf() to determine if file is IRAF or FITS
 * Oct 29 1998	Fix history setting
 * Nov 30 1998	Add version and help commands for consistency
 * Dec 30 1998	Write header without image if no image is present
 *
 * Mar  4 1999	Reset = for each keyword after setting value in header
 * Apr  2 1999	Add warning if too many files or keywords on command line
 * Apr  2 1999	Add -f and -m to change maximum number of files or keywords
 * Jul 14 1999  Read lists of BOTH keywords and files simultaneously
 * Jul 14 1999  Reallocate keyword array if too many in file
 * Jul 15 1999	Add capability of writing multi-line keywords a la IRAF
 * Jul 15 1999	Reallocate keyword and file lists if default limits exceeded
 * Oct 14 1999	Reallocate header if length is exceeded
 * Oct 22 1999	Drop unused variables after lint
 * Nov 17 1999	Fix bug which wrote a second entry for character values
 * Dec 20 1999	Add -d option to change date to ISO format
 *
 * Feb 17 2000	Fix bug reading last of assignments from file
 * Mar 22 2000	Use lt2fd() instead of getltime()
 * Apr 21 2000	Drop trailing spaces from character strings
 * May  1 2000	Drop -d option; it was unneeded for date-valued keywords
 * Jun  8 2000	Print revision message on only first of multiple files
 * Jun  8 2000	If no files or keywords specified, say so
 *
 * Aug 24 2001	If argument contains an equal sign, assume it not a file
 * Oct 31 2001	Do not write image if BITPIX or NAXIS is zero
 *
 * Jan  4 2002	Add / to add comments to keywords
 * Jan  4 2002	Replace newlines with nulls in keyword file buffer
 * Jan  4 2002	Allow null COMMENT and HISTORY values for spacing
 * Jan  4 2002	Allow stdin to be keyword input file if image file is set
 * Jan  9 2002	Add -s command to replace char with space in input value strings
 * Feb  5 2002	Add -l command to log files as they are processed
 * Nov  7 2002	If writing to a multiextension FITS file, write a new file
 *
 * Feb  4 2003	Fix bug dealing with quotes in @ files
 * Aug 21 2003	Read N-dimensional FITS image using fitsrfull()
 * Oct 28 2003	Increase output file name length from 128 to 1024
 * Oct 29 2003	Clean up command line and list file keyword input
 *
 * May  3 2004	Add code to overwrite header if revised one fits
 * Jul  1 2004	Fix bug deciding when to write extension header in place
 * Jul  1 2004	Change first extension if no extension specified
 * Aug 30 2004	Add option to add, subtract, multiply, or divide constant
 * Sep  7 2004	Set buffer pointers to null after freeing buffers
 * Sep  7 2004	Change input filename length from 128 to 1024
 *
 * Jan 25 2005	Drop newimage1; it is not needed
 * Mar  1 2005	Print header if first_file, not first
 * Mar  7 2005	Update comment as well as keyword vector if max exceeded
 * Mar  7 2005	Update comment as well as keyword vector if max exceeded
 */
