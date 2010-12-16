/* File gethead.c
 * January 17, 2006
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
#include "libwcs/wcscat.h"

#define MAXKWD 100
#define MAXFILES 2000
static int maxnkwd = MAXKWD;
static int maxncond = MAXKWD;
static int maxnfile = MAXFILES;

#define FILE_FITS 1
#define FILE_IRAF 2
#define FILE_ASCII 3

static void usage();
static void strclean();
extern char *GetFITShead();
static char nextnsp();
static int PrintValues();

const char *RevMsg = "GETHEAD WCSTools 3.6.4, 3 May 2006, Doug Mink (dmink@cfa.harvard.edu)";

static int verbose = 0;		/* verbose/debugging flag */
static int nfile = 0;
static int ndec = -9;
static int maxlfn = 0;
static int listall = 0;
static int listpath = 0;
static int tabout = 0;
static int tabpad = 0;
static int logfile = 0;
static int printhead = 0;
static int forceascii = 0;
static int version = 0;		/* If 1, print only program name and version */
static int printfill=0;		/* If 1, print ___ for unfound keyword values */
static int printfile=1;		/* If 1, print filename first if >1 files */
static int fillblank=0;		/* If 1, replace blanks in strings with _ */
static int keyeqval=0;		/* If 1, print keyword=value */
static int keyeqvaln=0;		/* If 1, print keyword=value<nl> */
static char *rootdir=NULL;	/* Root directory for input files */
static int ncond=0;		/* Number of keyword conditions to check */
static int condand=1;		/* If 1, AND comparisons, else OR */
static int toeol = 0;		/* If 1, return values from ASCII file to EOL */
static char **cond;		/* Conditions to check */
static char **ccond;		/* Condition characters */
static int nproc = 0;
static char *extensions;	/* Extension number(s) or name to read */
static char *extension;		/* Extension number or name to read */

int
main (ac, av)
int ac;
char **av;
{
    char *str;
    char **kwd;		/* Keywords to read */
    char **kwdnew;		/* Keywords to read */
    int nkwd = 0;
    int nkwd1 = 0;
    char **fn, **newfn;
    int *ft, *newft;
    int ifile;
    int lfn;
    char filename[256];
    char *name;
    FILE *flist, *fdk;
    char *listfile;
    char *ilistfile;
    char *klistfile;
    int ikwd, lkwd, i, j;
    char *kw, *kwe;
    char string[80];
    int nbytes;
    int filetype;
    int icond;
    int nfext;
    int nrmax=10;
    struct Range *erange;

    ilistfile = NULL;
    klistfile = NULL;
    extension = NULL;
    extensions = NULL;
    ncond = 0;
    nfile = 0;
    fn = (char **)calloc (maxnfile, sizeof(char *));
    ft = (int *)calloc (maxnfile, sizeof(int));
    kwd = (char **)calloc (maxnkwd, sizeof(char *));
    cond = (char **)calloc (maxncond, sizeof(char *));
    ccond = (char **)calloc (maxncond, sizeof(char *));

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
	if (*(str = *av)=='-') {
	    char c;
	    while (c = *++str)
	    switch (c) {

		case 'a': /* list file even if the keywords are not found */
		    listall++;
		    break;
	
		case 'b': /* Replace blanks with underscores */
		    fillblank++;
		    break;
	
		case 'c': /* Force reading as an ASCII file */
		    forceascii++;
		    break;
	
		case 'd': /* Root directory for input */
		    if (ac < 2)
			usage();
		    rootdir = *++av;
		    ac--;
		    break;

		case 'e': /* list keyword=value for input to sethead */
		    keyeqval++;
		    break;
	
		case 'f': /* Do not print file names */
		    printfile = 0;
		    break;

		case 'g': /* list keyword=value<nl> */
		    keyeqvaln++;
		    break;
	
		case 'h': /* Output column headings */
		    printhead++;
		    break;

                case 'i':       /* Log files changed */
                    logfile++;
                    break;
	
		case 'l': /* Return values to end of line */
		    toeol++;
		    break;
	
		case 'n': /* Number of decimal places in output */
		    if (ac < 2)
			usage();
		    ndec = (int) (atof (*++av));
		    ac--;
		    break;
	
		case 'o': /* OR conditions insted of ANDing them */
		    condand = 0;
		    break;
	
		case 'p': /* List file pathnames, not just file names */
		    listall++;
		    listpath++;
		    break;
	
		case 's': /* Do not pad output tab table */
		    tabpad = 0;
		    break;
	
		case 't': /* Output tab table */
		    tabout++;
		    break;

		case 'u': /* Always print ___ if keyword not found */
		    printfill++;
		    break;
	
		case 'v': /* More verbosity */
		    verbose++;
		    break;
	
		case 'x': /* FITS extension to read */
		    if (ac < 2)
			usage();
		    if (isnum (*(av+1)))
			extensions = *++av;
		    else {
			extensions = calloc (16, 1);
			strcpy (extensions, "1-1000");
			}
		    listall++;
		    ac--;
		    break;

		default:
		    usage();
		    break;
		}
	    }

	/* File containing a list of keywords or files */
	else if (*av[0] == '@') {
	    listfile = *av + 1;
	    if (isimlistd (listfile, rootdir)) {
		ilistfile = listfile;
		nfile = getfilelines (ilistfile);
		}
	    else if (isfilelist (listfile, rootdir)) {
		ilistfile = listfile;
		nfile = getfilelines (ilistfile);
		}
	    else {
		klistfile = listfile;
		nkwd1 = getfilelines (klistfile);
		if (nkwd1 > 0) {
		    if (nkwd+nkwd1 > maxnkwd) {
			maxnkwd = nkwd + nkwd1 + 32;
			kwdnew = (char **) calloc (maxnkwd, sizeof (char **));
		 	for (ikwd = 0; ikwd < nkwd; ikwd++)
			    kwdnew[ikwd] = kwd[ikwd];
			free (kwd);
			kwd = kwdnew;
			}
		    if ((fdk = fopen (klistfile, "r")) == NULL) {
			fprintf (stderr,"GETHEAD: File %s cannot be read\n",
				 klistfile);
			nkwd = 0;
			}
		    else {
			for (ikwd = 0; ikwd < nkwd1; ikwd++) {
			    kwd[nkwd] = (char *) calloc (32, 1);
			    first_token (fdk, 31, kwd[nkwd++]);
			    }
			fclose (fdk);
			}
		    }
		}
	    }

	/* Image file */
	else if (isfits (*av) || isiraf (*av)) {
	    if (nfile >= maxnfile) {
		maxnfile = maxnfile * 2;
		nbytes = maxnfile * sizeof (char *);
		newfn = (char **) calloc (maxnfile, sizeof (char *));
		for (i = 0; i < nfile; i++)
		    newfn[i] = fn[i];
		free (fn);
		fn = newfn;
		newft = (int *) calloc (maxnfile*2, sizeof (int));
		for (i = 0; i < nfile; i++)
		    newft[i] = ft[i];
		free (ft);
		ft = newft;
		maxnfile = maxnfile * 2;
		}
	    fn[nfile] = *av;
	    if (forceascii)
		ft[nfile] = FILE_ASCII;
	    else if (isfits (*av))
		ft[nfile] = FILE_FITS;
	    else
		ft[nfile] = FILE_IRAF;

	    if (listpath || (name = strrchr (fn[nfile],'/')) == NULL)
		name = fn[nfile];
	    else
		name = name + 1;
	    lfn = strlen (name);
	    if (lfn > maxlfn)
		maxlfn = lfn;
	    nfile++;
	    }

	/* Text file */
	else if (isfile (*av)) {
	    if (nfile >= maxnfile) {
		maxnfile = maxnfile * 2;
		nbytes = maxnfile * sizeof (char *);
		fn = (char **) realloc ((void *)fn, nbytes);
		nbytes = maxnfile * sizeof (int);
		ft = (int *) realloc ((void *)ft, nbytes);
		}
	    fn[nfile] = *av;
	    ft[nfile] = FILE_ASCII;

	    if (listpath || (name = strrchr (fn[nfile],'/')) == NULL)
		name = fn[nfile];
	    else
		name = name + 1;
	    lfn = strlen (name);
	    if (lfn > maxlfn)
		maxlfn = lfn;
	    nfile++;
	    }

	/* Condition */
	else if (strchr (*av, '=') != NULL || strchr (*av, '#') != NULL ||
		 strchr (*av, '>') != NULL || strchr (*av, '<') != NULL ) {
	    if (ncond >= maxncond) {
		maxncond = maxncond * 2;
		cond = (char **)realloc((void *)cond, maxncond*sizeof(void *));
		ccond = (char **)realloc((void *)cond, maxncond*sizeof(void *));
		}
	    cond[ncond] = *av;
	    ccond[ncond] = strchr (*av, '=');
	    if (ccond[ncond] == NULL)
		ccond[ncond] = strchr (*av, '#');
	    if (ccond[ncond] == NULL)
		ccond[ncond] = strchr (*av, '>');
	    if (ccond[ncond] == NULL)
		ccond[ncond] = strchr (*av, '<');
	    kwe = ccond[ncond];
	    if (kwe != NULL) {
		for (kw = cond[ncond]; kw < kwe; kw++) {
		    if (*kw > 96 && *kw < 123)
			*kw = *kw - 32;
		    }
		}
	    ncond++;
	    }

	/* Keyword */
	else {
	    if (nkwd >= maxnkwd) {
		maxnkwd = maxnkwd * 2;
		kwdnew = (char **) realloc ((void *)kwd, maxnkwd);
	 	for (ikwd = 0; ikwd < nkwd; ikwd++)
		    kwdnew[ikwd] = kwd[ikwd];
		free (kwd);
		kwd = kwdnew;
		}
	    kwd[nkwd] = *av;
	    nkwd++;
	    }
	}

    if (nkwd <= 0 && nfile <= 0 )
	usage ();
    else if (nkwd <= 0) {
	fprintf (stderr, "GETHEAD: no keywords specified\n");
	exit (1);
	}
    else if (nfile <= 0 ) {
	fprintf (stderr, "GETHEAD: no files specified\n");
	exit (1);
	}

    if (nkwd > 1 && tabpad)
	printfill = 1;
    if (nfile < 2 && !listall)
	printfile = 0;

    /* Print column headings if tab table or headings requested */
    if (printhead) {

	/* Print conditions in header */
	for (icond = 0; icond < ncond; icond++) {
	    if (verbose) {
		if (condand || icond == 0)
		    printf ("%s\n",cond[icond]);
		else
		    printf (" or %s\n",cond[icond]);
		}
	    else if (tabout) {
		if (condand || icond == 0)
		    printf ("condition	%s\n", cond[icond]);
		else
		    printf ("condition	or %s\n", cond[icond]);
		}
	    }

	if (printfile) {
	    printf ("FILENAME");
	    if (maxlfn > 8) {
		for (i = 8; i < maxlfn; i++)
		    printf (" ");
		}
	    if (tabout)
	    	printf ("	");
	    else
		printf (" ");
	    }

	/* Make keyword names upper case and print keyword names in header */
	for (ikwd = 0; ikwd < nkwd; ikwd++) {
	    lkwd = strlen (kwd[ikwd]);
	    kwe = kwd[ikwd] + lkwd;
	    if (ft[0] != FILE_ASCII) {
		for (kw = kwd[ikwd]; kw < kwe; kw++) {
		    if (*kw > 96 && *kw < 123)
			*kw = *kw - 32;
		    }
		}
	    printf ("%s",kwd[ikwd]);
	    if (verbose || ikwd == nkwd - 1)
	    	printf ("\n");
	    else if (tabout)
	    	printf ("	");
	    else
		printf (" ");
	    }

	/* Print field-defining hyphens if tab table output requested */
	if (printhead && tabout) {
	    if (printfile) {
		strcpy (string, "-----------------------------------");
		if (maxlfn > 8)
		    string[maxlfn] = (char) 0;
		else
		    string[8] = (char) 0;
		printf ("%s",string);
		if (verbose || ikwd == nkwd - 1)
		    printf ("\n");
		else
		    printf ("	");
		}

	    for (ikwd = 0; ikwd < nkwd; ikwd++) {
		strcpy (string, "--------------");
		lkwd = strlen (kwd[ikwd]);
		string[lkwd] = (char) 0;
		printf ("%s",string);
		if (verbose || ikwd == nkwd - 1)
		    printf ("\n");
		else
		    printf ("	");
		string[lkwd] = '-';
		}
	    }
	}

    /* Open file containing a list of images, if there is one */
    if (ilistfile != NULL) {
	if ((flist = fopen (ilistfile, "r")) == NULL) {
	    fprintf (stderr,"GETHEAD: Image list file %s cannot be read\n",
		     ilistfile);
	    usage ();
	    }
	}

    /* Check extensions for range and set accordingly */
    if (extensions != NULL) {
	if (isrange (extensions)) {
	    erange = RangeInit (extensions, nrmax);
	    nfext = rgetn (erange);
	    extension = calloc (1, 8);
	    }
	else {
	    extension = extensions;
	    if (extension)
		nfext = 1;
	    }
	}
    else {
	extension = NULL;
	nfext = 0;
	}

    /* Read through headers of images */
    for (ifile = 0; ifile < nfile; ifile++) {
	if (ilistfile != NULL) {
	    first_token (flist, 254, filename);
	    if (forceascii)
		filetype = FILE_ASCII;
	    else if (isiraf (filename))
		filetype = FILE_IRAF;
	    else if (isfits (filename))
		filetype = FILE_FITS;
	    else
		filetype = FILE_ASCII;
	    if (nfext > 1) {
		rstart (erange);
		for (i = 0; i < nfext; i++) {
		    j = rgeti4 (erange);
		    sprintf (extension, "%d", j);
		    if (PrintValues (filename, filetype, nkwd, kwd))
			break;
		    }
		}
	    else
		PrintValues (filename, filetype, nkwd, kwd);
	    }
	else {
	    if (nfext > 1) {
		rstart (erange);
		for (i = 0; i < nfext; i++) {
		    j = rgeti4 (erange);
		    sprintf (extension, "%d", j);
		    if (PrintValues (fn[ifile], ft[ifile], nkwd, kwd))
			break;
		    }
		}
	    else
		(void) PrintValues (fn[ifile], ft[ifile], nkwd, kwd);
	    }

	if (verbose)
	    printf ("\n");

	/* Log the processing of this file, if requested */
	if (logfile) {
	    nproc++;
	    fprintf (stderr, "%d: %s processed.\r", nproc, filename);
	    }
	}
    if (ilistfile != NULL)
	fclose (flist);

    return (0);
}

static void
usage ()
{
    fprintf (stderr,"%s\n",RevMsg);
    if (version)
	exit (-1);
    fprintf (stderr,"Print FITS or IRAF header keyword values\n");
    fprintf(stderr,"Usage: gethead [-abhoptv][-d dir][-f num][-m num][-n num] file1.fit ... filen.fits kw1 kw2 ... kwn\n");
    fprintf(stderr,"  or : gethead [-abhoptv][-d dir][-f num][-m num][-n num] @filelist kw1 kw2 ... kwn\n");
    fprintf(stderr,"  or : gethead [-abhoptv][-d dir][-f num][-m num][-n num] file1.fit ... filen.fits @keywordlist\n");
    fprintf(stderr,"  or : gethead [-abhoptv][-d dir][-f num][-m num][-n num] @filelist @keywordlist\n");
    fprintf(stderr,"  -a: List file even if keywords are not found or null\n");
    fprintf(stderr,"  -b: Replace blanks in strings with underscores\n");
    fprintf(stderr,"  -c: Read all files as plain ASCII\n");
    fprintf(stderr,"  -d: Root directory for input files (default is cwd)\n");
    fprintf(stderr,"  -e: Output keyword=value's on one line per file\n");
    fprintf(stderr,"  -f: Never print filenames (default is print if >1)\n");
    fprintf(stderr,"  -g: Output keyword=value's on one line per keyword\n");
    fprintf(stderr,"  -h: Print column headings\n");
    fprintf(stderr,"  -n: Number of decimal places in numeric output\n");
    fprintf(stderr,"  -o: OR conditions instead of ANDing them\n");
    fprintf(stderr,"  -p: Print full pathnames of files\n");
    fprintf(stderr,"  -s: Do not pad tab-separated table with spaces\n");
    fprintf(stderr,"  -t: Output in tab-separated table format\n");
    fprintf(stderr,"  -u: Always print ___ if keyword not found or null value\n");
    fprintf(stderr,"  -v: Verbose\n");
    fprintf(stderr,"  -x [range]: Read header for these extensions (no arg=all)\n");
    exit (1);
}


static int
PrintValues (name, filetype, nkwd, kwd)

char	*name;		/* Name of FITS or IRAF image file */
int	filetype;	/* Type of file (FILE_FITS, FILE_IRAF, FILE_ASCII) */
int	nkwd;		/* Number of keywords for which to print values */
char	*kwd[];		/* Names of keywords for which to print values */

{
    char *header;	/* FITS image header or contents of ASCII file */
    char *cstr, *str;
    int iraffile;
    char fnform[8];
    char string[80];
    char temp[1028];
    char keyword[16];
    char *filename;
    char outline[1000];
    char padline[1000];
    char mstring[800];
    char ctab = (char) 9;
    char cspace = (char) 32;
    char pchar;
    char *kw, *kwe, *filepath;
    char *ext, *namext, cext;
    int ikwd, lkwd, nfound, notfound, nch;
    int jval, jcond, icond, i, j, lout, lstr;
    double dval, dcond, dnum;
    char tcond;
    char cvalue[64], *cval;
    char numstr[32], numstr1[32];
    int pass, iwcs, nwild;

    namext = NULL;
    ext = strchr (name, ',');
    if (extension && !ext) {
	nch = strlen (name) + 2 + strlen (extension);
	namext = (char *) calloc (1, nch);
	strcpy (namext, name);
	strcat (namext, ",");
	strcat (namext, extension);
	}
    else {
	nch = strlen (name) + 1;
	namext = (char *) calloc (1, nch);
	strcpy (namext, name);
	}
    ext = strchr (namext, ',');

    if (rootdir) {
	nch = strlen (rootdir) + strlen (namext) + 1;
	filepath = (char *) calloc (1, nch);
	strcat (filepath, rootdir);
	strcat (filepath, "/");
	strcat (filepath, namext);
	}
    else {
	nch = strlen (namext) + 1;
	filepath = (char *) calloc (1, nch);
	strcpy (filepath, namext);
	}

    if (!tabout && listall)
	printfill = 1;

    /* Read ASCII file into buffer */
    if (filetype == FILE_ASCII) {
	if ((header = getfilebuff (filepath)) == NULL)
	    return (-1);
	}

    /* Retrieve FITS header from FITS or IRAF .imh file */
    else if ((header = GetFITShead (filepath, verbose)) == NULL)
	return (-1);

    if (verbose) {
	fprintf (stderr,"%s\n",RevMsg);
	fprintf (stderr,"Print Header Parameter Values from ");
	hgeti4 (header, "IMHVER", &iraffile );
	if (filetype == FILE_ASCII)
	    fprintf (stderr,"ASCII file %s\n", namext);
	else if (filetype == FILE_IRAF)
	    fprintf (stderr,"IRAF image file %s\n", namext);
	else
	    fprintf (stderr,"FITS image file %s\n", namext);
	}

    /* Find file name */
    if (listpath || (filename = strrchr (namext,'/')) == NULL)
	filename = namext;
    else if (rootdir)
	filename = namext;
    else
	filename = filename + 1;

    if (printfile) {
	if (ext) {
	    cext = *ext;
	    *ext = (char) 0;
	    }
	sprintf (fnform, "%%-%ds", maxlfn);
	sprintf (outline, fnform, filename);
	if (ext) {
	    *ext = '[';
	    sprintf (string, "%s]", ext);
	    strcat (outline, string);
	    *ext = cext;
	    }
	if (tabout)
	    strcat (outline, "	");
	else
	    strcat (outline, " ");
	}
    else
	outline[0] = (char) 0;

    /* Check conditions */
    pass = 0;
    if (ncond > 0) {
	for (icond = 0; icond < ncond; icond++) {
	    if (condand)
		pass = 0;
	    tcond = *ccond[icond];

	    /* Extract test value from comparison string */
	    *ccond[icond] = (char) 0;
	    cstr = ccond[icond]+1;
	    if (strchr (cstr, ':')) {
		dnum = str2dec (cstr);
		num2str (numstr, dnum, 0, 7);
		cstr = numstr;
		}
	    strclean (cstr);

	    /* Read comparison value from header */
	    if (!hgets (header, cond[icond], 64, cvalue))
		continue;
	    cval = cvalue;
	    if (strchr (cval, ':')) {
		dnum = str2dec (cval);
		num2str (numstr1, dnum, 0, 7);
		cval = numstr1;
		}
	    strclean (cval);

	    /* Compare floating point numbers */
	    if (isnum (cstr) == 2 && isnum (cval)) {
		*ccond[icond] = tcond;
		dcond = atof (cstr);
		dval = atof (cval);
		if (tcond == '=' && dval == dcond)
		    pass = 1;
		if (tcond == '#' && dval != dcond)
		    pass = 1;
		if (tcond == '>' && dval > dcond)
		    pass = 1;
		if (tcond == '<' && dval < dcond)
		    pass = 1;
		}

	    /* Compare integers */
	    else if (isnum (cstr) == 1 && isnum (cval)) {
		*ccond[icond] = tcond;
		jcond = atoi (cstr);
		jval = atoi (cval);
		if (tcond == '=' && jval == jcond)
		    pass = 1;
		if (tcond == '#' && jval != jcond)
		    pass = 1;
		if (tcond == '>' && jval > jcond)
		    pass = 1;
		if (tcond == '<' && jval < jcond)
		    pass = 1;
		}

	    /* Compare strings (only equal or not equal */
	    else {
		*ccond[icond] = tcond;
		if (tcond == '=' && !strcmp (cstr, cval))
		    pass = 1;
		if (tcond == '#' && strcmp (cstr, cval))
		    pass = 1;
		}
	    if (condand && !pass)
		return (0);
	    }
	if (!pass)
	    return (0);
	}

    /* Read keywords from header */
    nfound = 0;
    notfound = 0;
    for (ikwd = 0; ikwd < nkwd; ikwd++) {

	/* IF FITS header, look for upper case keywords only */
	if (filetype != FILE_ASCII) {
	    lkwd = strlen (kwd[ikwd]);
	    kwe = kwd[ikwd] + lkwd;
	    for (kw = kwd[ikwd]; kw < kwe; kw++) {
		if (*kw > 96 && *kw < 123)
		    *kw = *kw - 32;
		}
	    }
	strcpy (keyword, kwd[ikwd]);

	/* Read keyword value from ASCII file */
	if (toeol)
	    lstr = -80;
	else
	    lstr = 80;
	if (filetype == FILE_ASCII &&
	    agets (header, keyword, lstr, string)) {
	    str = string;
	    strclean (str);
	    if (ndec > -9 && isnum (str) && strchr (str, '.'))
		num2str (str, atof(str), 0, ndec);
	    if (keyeqvaln)
		printf ("%s = %s\n", keyword, str);
	    else if (keyeqval) {
		sprintf (temp, " %s=%s", keyword, str);
		strcat (outline, temp);
		}
	    else if (verbose) {
		sprintf (temp, " %s = %s", keyword, str);
		strcat (outline, temp);
		}
	    else
		strcat (outline, str);
	    nfound++;
	    }

	/* Read all keywords from multiple WCS's if wildcarded with @ */
	else if (keyword[lkwd-1] == '@') {
	    keyword[lkwd-1] = (char) 0;
	    nwild = 0;
	    for (iwcs = -1; iwcs < 26; iwcs++) {
		if (iwcs < 0)
		    keyword[lkwd-1] = (char)0;
		else
		    keyword[lkwd-1] = (char)(iwcs+65);
		if (hgets (header, keyword, 80, string)) {
		    if (string[0] == '#' && isnum (string+1)) {
			lstr = strlen (string);
			for (i = 0; i < lstr; i++)
			    string[i] = string[i+1];
			}
		    str = string;
		    strclean (str);
		    if (ndec > -9 && isnum (str) && strchr (str, '.'))
			num2str (string, atof(str), 0, ndec);
		    if (verbose) {
			if (strchr (str,' '))
			    printf ("%s = \"%s\"\n", keyword, str);
			else
			    printf ("%s = %s\n", keyword, str);
			}
		    else if (keyeqvaln) {
			if (strchr (str,' '))
			    printf ("%s=\"%s\"\n", keyword, str);
			else
			    printf ("%s=%s\n", keyword, str);
			}
		    else if (keyeqval) {
			sprintf (temp, " %s=%s", keyword, str);
			strcat (outline, temp);
			}
		    else {
			if (nwild) {
			    if (tabout)
				strcat (outline, "	");
			    else
				strcat (outline, " ");
			    }
			strcat (outline, str);
			nwild++;
			}
		    nfound++;
		    }
		}
	    }

	/* Read one FITS keyword value */
	else if (hgets (header, keyword, 80, string)) {
	    if (string[0] == '#' && isnum (string+1)) {
		lstr = strlen (string);
		for (i = 0; i < lstr; i++)
		    string[i] = string[i+1];
		}
	    str = string;
	    strclean (str);
	    if (ndec > -9 && isnum (str) && strchr (str, '.'))
		num2str (string, atof(str), 0, ndec);
	    if (verbose) {
		if (strchr (str,' '))
		    printf ("%s = \"%s\"\n", keyword, str);
		else
		    printf ("%s = %s\n", keyword, str);
		}
	    else if (keyeqvaln) {
		if (strchr (str,' '))
		    printf ("%s=\"%s\"\n", keyword, str);
		else
		    printf ("%s=%s\n", keyword, str);
		}
	    else if (keyeqval) {
		sprintf (temp, " %s=%s", keyword, str);
		strcat (outline, temp);
		}
	    else if (strlen (str) < 1)
		strcat (outline, "___");
	    else
		strcat (outline, str);
	    nfound++;
	    }

	/* Read IRAF-style multiple-line keyword value */
	else if (hgetm (header, keyword, 600, mstring)) {
	    if (verbose) {
		if (strchr (mstring,' '))
		    printf ("%s = \"%s\"\n", keyword, mstring);
		else
		    printf ("%s = %s\n", keyword, mstring);
		}
	    else if (keyeqvaln) {
		if (strchr (mstring,' '))
		    printf ("%s=\"%s\"\n", keyword, mstring);
		else
		    printf ("%s=%s\n", keyword, mstring);
		}
	    else if (keyeqval) {
		sprintf (temp, " %s=%s", keyword, mstring);
		strcat (outline, temp);
		}
	    else
		strcat (outline, mstring);
	    nfound++;
	    }

	else if (keyeqvaln)
	    printf ("%s not found\n", keyword);
	else if (printfill)
	    strcat (outline, "___");
	else
	    notfound = 1;

	if (!keyeqvaln && ikwd < nkwd-1) {
	    if (tabout)
		strcat (outline, "	");
	    else
		strcat (outline, " ");
	    }
	}

    /* Print line of keywords */
    if (!verbose && !keyeqvaln && (nfound > 0 || printfill || listall) && 
	(nfile < 2 || nfound > 0 || listall)) {

	/* Remove spaces used to pad tab-separated tables for readability */
	if (tabout && !tabpad) {
	    lout = strlen (outline);
	    strcpy (padline, outline);
	    for (i = 0; i < 1000; i++)
		outline[i] = (char) 0;
	    j = 0;
	    pchar = ctab;
	    for (i = 0; i < lout; i++) {
		if (padline[i] == cspace && pchar == ctab)
		    continue;
		if (padline[i] == cspace && nextnsp (padline+i) == ctab)
		    continue;
		if (padline[i] == cspace && nextnsp (padline+i) == (char) 0)
		    continue;
		outline[j++] = padline[i];
		pchar = padline[i];
		}
	    }

	printf ("%s\n", outline);
	}

    free (header);
    return (0);
}


/* Return next character in string which is not a space */
static char
nextnsp (string)

char *string;

{
    int lstring;
    char *schar;

    schar = string;
    while (*schar == ' ')
	schar++;
    return (*schar);
}


/* Remove exponent, leading #, and/or trailing zeroes, if reasonable */
static void
strclean (string)

char *string;

{
    char *sdot, *s, *strend, *str, ctemp, *slast;
    int ndek, lstr, i;

    /* If number, ignore leading # and remove trailing non-numeric character */
    if (string[0] == '#') {
	strend = string + strlen (string);
	str = string + 1;
	strend = str + strlen (str) - 1;
	ctemp = *strend;
	if (!isnum (strend))
	    *strend = (char) 0;
	if (isnum (str)) {
	    strend = string + strlen (string);
	    for (str = string; str < strend; str++)
		*str = *(str + 1);
	    }
	else
	    *strend = ctemp;
	}

    /* Remove positive exponent if there are enough digits given */
    if (strsrch (string, "E+") != NULL) {
	lstr = strlen (string);
	ndek = (int) (string[lstr-1] - 48);
	ndek = ndek + (10 * ((int) (string[lstr-2] - 48)));
	if (ndek < lstr - 7) {
	    lstr = lstr - 4;
	    string[lstr] = (char) 0;
	    string[lstr+1] = (char) 0;
	    string[lstr+2] = (char) 0;
	    string[lstr+3] = (char) 0;
	    sdot = strchr (string, '.');
	    if (ndek > 0 && sdot != NULL) {
		for (i = 1; i <= ndek; i++) {
		    *sdot = *(sdot+1);
		    sdot++;
		    *sdot = '.';
		    }
		}
	    }
	}

    /* Remove trailing zeroes if they are not significant */
    if (strchr (string, '.') != NULL &&
	strsrch (string, "E-") == NULL &&
	strsrch (string, "E+") == NULL &&
	strsrch (string, "e-") == NULL &&
	strsrch (string, "e+") == NULL) {
	lstr = strlen (string);
	s = string + lstr - 1;
	while (*s == '0' && lstr > 1) {
	    if (*(s - 1) != '.') {
		*s = (char) 0;
		lstr --;
		}
	    s--;
	    }
	}

    /* Remove trailing decimal point */
    lstr = strlen (string);
    s = string + lstr - 1;
    if (*s == '.')
	*s = (char) 0;

    /* Replace embedded blanks with underscores, if requested to */
    if (fillblank) {
	lstr = strlen (string);
	slast = string + lstr;
	for (s = string; s < slast; s++) {
	    if (*s == ' ') *s = '_';
	    }
	}

    return;

}

/* Sep  4 1996	New program
 * Oct  8 1996	Add newline after file name in verbose mode
 *
 * Feb 21 1997  Check pointers against NULL explicitly for Linux
 * Nov  2 1997	Allow search of multiple files from command line or list file
 * Nov 12 1997	Add option to print tab tables or column headings
 * Dec 12 1997	Read IRAF version 2 .imh files
 *
 * Jan  5 1998	Do not print line unless keyword is found
 * Jan 27 1998	Clean up scientific notation and trailing zeroes
 * Mar 12 1998	Read IRAF multi-line keyword values
 * May 27 1998	Include fitsio.h instead of fitshead.h
 * Jun  2 1998  Fix bug in hput()
 * Jun  3 1998	Add -p option
 * Jun 18 1998	Print tab table heading only if -h option used
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 * Sep  1 1998	Set number of decimal places in floating output with -n
 * Oct  5 1998	Check more carefully for FITS and IRAF files on command line
 * Nov 30 1998	Add version and help commands for consistency
 *
 * Feb 12 1999	Print null string if single keyword is not found
 * Feb 17 1999	Add -d option to set root input directory
 * Apr  1 1999	Add warning if too many files or keywords on command line
 * Apr  2 1999	Add -f and -m to change maximum number of files or keywords
 * Jun  9 1999	Initialize outline so Linux doesn't print garbage
 * Jun 21 1999	Fix bug so that -a option works, always printing filename
 * Jul 13 1999	Use only first token from line of list as filename
 * Jul 14 1999	Read lists of BOTH keywords and files simultaneously
 * Jul 15 1999	Reallocate keyword and file lists if default limits exceeded
 * Oct 22 1999	Drop unused variables after lint
 * Nov 16 1999	Add -u to always print underscores if keyword not found
 * Nov 19 1999	Add -f to never print filenames
 * Nov 30 1999	Fix so no file name is printed if only one file unless -a
 * Nov 30 1999	Cast realloc's
 *
 * Feb 24 2000	Add option to output keyword=value
 * Mar  1 2000	Add option to read ASCII files with keyword=value in them
 * Mar 17 2000	Add conditions
 * Mar 20 2000	Drop leading # from numbers
 * Mar 21 2000	Add -b option to replace blanks with underscores
 * Jun  8 2000	If no keywords or files specified, say so
 * Jun 12 2000	If -p is set, print all file names (-a)
 *
 * Feb 21 2001	Add @ wildcard option for multiple WCS keywords
 * Feb 27 2001	Add space or tab between wildcard multiple WCS keyword values
 * Apr 23 2001	Add -g for keyword=val<lf> and print same way if -v
 * Sep 25 2001	Allow file and command line lists of ASCII files
 *
 * Jan 31 2002	Fix bug to always add underlines if printing tab table headers
 * Feb 26 2002	Return values to end of line if -l option set
 * Apr  2 2002	Add option to log processing to stderr
 * Apr 24 2002	Add -c option to force files to be read as plain ASCII
 * May 23 2002	Add -x option to specify extension to read
 * May 31 2002	Add -s option to drop space-padding for tab table readability
 * Jun  6 2002	Allow -x to specify range of numbered extensions
 * Jun 18 2002	List filename(s) if -x used
 * Jun 19 2002	Add verbose argument to GetFITShead()
 * Jun 20 2002	If -x and no argument, read all extensions
 *
 * Feb  5 2003	Set nfext tp zero if no extensions
 * Mar 25 2003	If null keyword value and padding on, print ___
 * Jul 17 2003	Add root directory argumeht to isfilelist()
 * Oct 29 2003	Allow keyword specification from both list file and command line
 * Nov 18 2003	Fix strclean() to keep all of exponents (found by Anthony Miceli)
 *
 * Apr 15 2004	Allow e as well as E for exponents
 * Dec 20 2004	If printing filename for tab-separated file, head with PathName
 *
 * Jan 14 2005	Change PathName back to filename
 * Jun  9 2005	Fix bugs dealing with large numbers of files or keywords
 * Jun 15 2005	Write one-per-line output in SETHEAD input file format
 * Jul 18 2005	Do not write one-per-line verbose output in SETHEAD input file format
 *
 * Jan 17 2006	ALWAYS print line if -a
 */
