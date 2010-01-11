/* File delhead.c
 * April 24, 2006
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

#define MAXKWD 500
#define MAXFILES 1000
static int maxnkwd = MAXKWD;
static int maxnfile = MAXFILES;
#define MAXNEW 1024

static void usage();
static void DelKeywords();
extern char *fitserrmsg;

static int verbose = 0;		/* verbose/debugging flag */
static int newimage = 0;
static char *RevMsg = "DELHEAD WCSTools 3.6.4, 3 May 2006, Doug Mink (dmink@cfa.harvard.edu)";
static int version = 0;		/* If 1, print only program name and version */
static int logfile = 0;
static int nproc = 0;
static int first_file = 1;

int
main (ac, av)
int ac;
char **av;
{
    char *str;
    char **kwd;
    int nkwd = 0;
    int nkwd1 = 0;
    int ikwd;
    char **fn;
    int nfile = 0;
    int ifile;
    char filename[128];
    FILE *flist, *fdk;
    char *listfile;
    char *ilistfile;
    char *klistfile;
    char **kwdnew;

    ilistfile = NULL;
    klistfile = NULL;
    nfile = 0;
    fn = (char **)calloc (maxnfile, sizeof(char *));
    kwd = (char **)calloc (maxnkwd, sizeof(char *));

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage ();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage ();
	}

    /* crack arguments */
    for (av++; --ac > 0; av++) {
	if (*(str = *av) == '-') {
	    char c;
	    while (c = *++str)
	    switch (c) {
	
		case 'l':	/* Log files changed */
		    logfile++;
		    break;

		case 'n':	/* write new file */
		    newimage++;
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
	    listfile = *av + 1;
	    if (isimlist (listfile)) {
		ilistfile = listfile;
		nfile = getfilelines (ilistfile);
		}
	    else {
		klistfile = listfile;
		nkwd1 = getfilelines (klistfile);
		if (nkwd1 > 0) {
		    if (nkwd1 + nkwd > maxnkwd) {
			maxnkwd = maxnkwd + nkwd1 + 32;
			kwdnew = (char **)calloc (maxnkwd, sizeof(char *));
			for (ikwd = 0; ikwd < nkwd; ikwd++)
			    kwdnew[ikwd] = kwd[ikwd];
			free (kwd);
			kwd = kwdnew;
			}
		    if ((fdk = fopen (klistfile, "r")) == NULL) {
			fprintf (stderr,"DELHEAD: File %s cannot be read\n",
				 klistfile);
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
		fn = (char **) realloc ((void *)fn, maxnfile);
		}
	    fn[nfile] = *av;
	    nfile++;
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
	fprintf (stderr, "DELHEAD: no keywords specified\n");
	exit (1);
	}
    else if (nfile <= 0 ) {
	fprintf (stderr, "DELHEAD: no files specified\n");
	exit (1);
	}

    /* Delete keyword values one file at a time */

    /* Open file containing a list of images, if there is one */
    if (ilistfile != NULL) {
	if ((flist = fopen (ilistfile, "r")) == NULL) {
	    fprintf (stderr,"DELHEAD: Image list file %s cannot be read\n",
		     ilistfile);
	    usage ();
	    }
	}

    /* Read through headers of images */
    for (ifile = 0; ifile < nfile; ifile++) {
	if (ilistfile != NULL) {
	    first_token (flist, 254, filename);
	    DelKeywords (filename, nkwd, kwd);
	    }
	else
	    DelKeywords (fn[ifile], nkwd, kwd);
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
    fprintf (stderr,"Delete FITS or IRAF header keyword entries\n");
    fprintf(stderr,"Usage: [-nv] file1.fits [ ... filen.fits] kw1 [... kwn]\n");
    fprintf(stderr,"  or : [-nv] @listfile kw1 [... kwn]\n");
    fprintf(stderr,"  or : [-nv] file1.fits [ ... filen.fits] @keylistfile\n");
    fprintf(stderr,"  or : [-nv] @listfile @keylistfile\n");
    fprintf(stderr,"  -n: write new file\n");
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}


static void
DelKeywords (filename, nkwd, kwd)

char	*filename;	/* Name of FITS or IRAF image file */
int	nkwd;		/* Number of keywords to delete */
char	*kwd[];		/* Names of those keywords */

{
    char *header;	/* FITS image header */
    char *image;	/* FITS image */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    int nblold, nblnew;	/* Number of FITS blocks (=2880 bytes) in header */
    char *irafheader;	/* IRAF image header */
    int iraffile;	/* 1 if IRAF image, 0 if FITS image */
    int lext, lroot, naxis;
    char newname[MAXNEW];
    char *ext, *fname, *imext, *imext1;
    char *kw, *kwl;
    char echar;
    int ikwd;
    int ibhead;		/* Number of bytes to skip to header */
    int fdr, fdw, ipos, nbr, nbw, bitpix, imageread;
    int nbold, nbnew;

    image = NULL;
    header = NULL;

    /* Open IRAF image if .imh extension is present */
    if (isiraf (filename)) {
	iraffile = 1;
	if ((irafheader = irafrhead (filename, &lhead)) != NULL) {
	    if ((header = iraf2fits (filename, irafheader, lhead, &nbhead)) == NULL) {
		fprintf (stderr, "Cannot translate IRAF header %s/n",filename);
		free (irafheader);
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
	fprintf (stderr,"Delete Header Parameter Entries from ");
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s\n", filename);
	else
	    fprintf (stderr,"FITS image file %s\n", filename);
	first_file = 0;
	}

    if (nkwd < 1)
	return;

    nbold = fitsheadsize (header);

    /* Remove directory path and extension from file name */
    fname = strrchr (filename, '/');
    if (fname)
	fname = fname + 1;
    else
	fname = filename;

    if (strchr (fname, ',') || strchr (fname,'['))
	setheadshrink (0);

    /* Delete keywords one at a time */
    for (ikwd = 0; ikwd < nkwd; ikwd++) {

	/* Make keyword all upper case */
	kwl = kwd[ikwd] + strlen (kwd[ikwd]);
	for (kw = kwd[ikwd]; kw < kwl; kw++) {
	    if (*kw > 96 && *kw < 123)
		*kw = *kw - 32;
	    }

	/* Delete keyword */
	if (hdel (header, kwd[ikwd]) && verbose)
	    printf ("%s: %s deleted\n", filename, kwd[ikwd]);
	}

    /* Compare size of output header to size of input header */
    nbnew = fitsheadsize (header);
    nblnew = (int) (0.98 + (double) nbnew / 2880.0);
    nblold = (int) (0.98 + (double) nbold / 2880.0);
    if (nbnew > nbold && naxis == 0 && bitpix != 0) {
	if (verbose)
	    fprintf (stderr, "Rewriting primary header, copying rest of file\n");
	newimage = 1;
	}

    /* Make up name for new FITS or IRAF output file */
    if (newimage) {
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

	/* Keep name */
	strcpy (newname, filename);

	/* Set image extension if there is one */
	imext = strchr (filename, ',');
	if (imext == NULL)
	    imext = strchr (filename, '[');

	/* Add extension if modifying extension header */
	if (imext == NULL && ksearch (header,"XTENSION")) {
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
    else if (nblnew == nblold && !newimage) {
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

    if (header != NULL)
	free (header);
    if (image != NULL)
	free (image);
    return;
}

/* Jul 27 1998	New program
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 * Aug 14 1998	If changing primary header, write out entire input file
 * Oct  5 1998	Allow header changes even if no data is present
 * Oct  5 1998	Use isiraf() and isfits() to check for data file
 * Nov 30 1998	Add version and help commands for consistency
 *
 * Mar  2 1999	Add option to delete list of keyword names from file
 * Apr  2 1999	Add warning if too many files or keywords on command line
 * Apr  2 1999	Add -f and -m to change maximum number of files or keywords
 * Jul 14 1999	Read lists of BOTH keywords and files simultaneously
 * Jul 14 1999	Reallocate keyword array if too many in file
 * Jul 15 1999	Reallocate keyword and file lists if default limits exceeded
 * Sep 29 1999	Change maximum number of keywords from 100 to 500
 * Oct 21 1999	Drop unused variables after lint
 * Nov 29 1999	Fix usage command list
 * Nov 30 1999	Cast realloc's
 *
 * Jun  8 2000	If no files or keywords specified, say so
 *
 * Dec 16 2002	Fix bug so arbitrary number of keywords can be deleted
 *
 * Aug 21 2003	Use fitsrfull() to deal with n dimensional FITS images
 * Oct 29 2003	Keep count of keywords correctly when reading them from file
 *
 * May  6 2004	Allow keywords to be deleted from extension headers
 * Jul  1 2004	Do not drop lines from multi-extension headers
 * Jul  1 2004	Change first extension if no extension specified
 *
 * Jan 12 2005	Write over unread image only if number of header blocks same
 * Mar  1 2005	Print program version only on first file if looping
 * Jun 10 2005	Fix bug dealing with large numbers of keywords
 *
 * Apr 26 2006	Avoid freeing alread-freed image buffers
 */
