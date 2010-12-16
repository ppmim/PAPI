/* mosaic.c -- use correlate.c to find translation offsets between images 
 *             If the offsets from the first image are above a certain limit,
 *             it tryes to correlate with the previous image in the list */
	

/*
 * Input filelist is a list of SExtractor OBJECTS output frames (FITS files), 
 * with the first file called the reference frame.  Input frames are expected
 * to be similar except for a translation.  readwcs() assumes that North is
 * up and East is left.
 *
 * Procedure:  
 *   Read the object pixels from the reference frame into x,y,pixel lists
 *   For each additional frame in filelist:
 *     Guess the pixel offsets verse the reference frame using header RA,DEC
 *     Correlate with reference frame pixel list to measure precise offset
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "correlate.h"
#include "eprintf.h"
#include "fitsIO.h"
#include "listIO.h"

#define HWID 10.0     /* default half-width of cross-corr search box, arcsec */
#define MAXNFILES 399
#define MAXNLIST 99999L
#define MINFRAC  0.1      /* 0.15; correlation failed 
							 if overlap*comarea < MINFRAC */
#define MINAREA0 0.25	/* min fractional common area to use the first image */
#define MINAREA  0.02	/* min fractional common area to make a computation */

static char *fn[MAXNFILES];              /* FITS OBJECTS frames filenames */
static int n, x[MAXNLIST], y[MAXNLIST];       /* list of object pixels */
static float p[MAXNLIST];
static int xprev[MAXNLIST], yprev[MAXNLIST];  /* list of object pixels */
static float pprev[MAXNLIST];

static float readwcs(char *fn, int *ixoff, int *iyoff);
static char *chomp(char *fn);
static void usage(void);

int main(int argc, char *argv[])
{
    int i, nx, ny, ixoff, iyoff, nfiles, prev, ix, iy, j, m;
	int hwid, hwid2, hwid4;
    float xoff, yoff, frac, scale, *img0, *img, *pimg, corrx, corry;
    float offx[MAXNFILES], offy[MAXNFILES];	/* list of all offsets */
	float comarea0, comarea[MAXNFILES], maxcomarea;	/* expected common area */

    if (argc != 2 && argc != 3)
        usage();

    if ((nfiles = readlist(argv[1], fn, NULL, NULL, NULL, MAXNFILES)) < 1)
        eprintf("%s: no valid files\n", argv[0]);

    img0 = readfits(fn[0], &nx, &ny, NULL, NULL);

    scale = readwcs(fn[0], NULL, NULL);              /* arcsec per pixel */

    hwid = (argc == 3) ? (atof(argv[2]) / scale + 0.5) : (HWID / scale + 0.5);

    for (i = 0, n = 0; i < nx * ny; i++) {    /* make pixel list from refimg */
        if (img0[i] > 0 && n < MAXNLIST) {
            x[n] = i % nx;
            y[n] = i / nx;
            p[n++] = img0[i];
        }
    }

    free(img0);

    if (n >= MAXNLIST)
        eprintf("%s: n >= MAXNLIST, increase detection threshold\n", argv[0]);

    if (n < 1)
        eprintf("%s: found no object pixels in %s\n", argv[0], fn[0]);

    /* printf("%s %f %f %f\n", chomp(fn[0]), 0.0, 0.0, 1.0); reference frame */
    printf("%s %f %f %f\n", fn[0], 0.0, 0.0, 1.0); /* reference frame */
    fprintf(stderr, "%s %f %f %f\n", fn[0], 0.0, 0.0, 1.0); /* reference frame */

    corrx = 0.;
    corry = 0.;

    offx[0] = 0.0 ;
    offy[0] = 0.0 ;

    for (i = 1; i < nfiles; i++) {                        /* other frames... */
		img = readfits(fn[i], &nx, &ny, NULL, NULL);

		(void) readwcs(fn[i], &ixoff, &iyoff);           /* get offset guess */

		ixoff = ixoff+corrx;	/* use the previous corrections */
		iyoff = iyoff+corry;
		
		ix = abs(ixoff);
		iy = abs(iyoff);
		if (ix<nx && iy<ny)
			comarea0  = (float) ((nx-ix)*(ny-iy)) / (nx*ny) ;
		else 
			comarea0  = 0.0 ;
		
		
		/* printf ("IX= %d  IY= %d  COM= %f \n", ixoff, iyoff, comarea0); */

		if (comarea0 > MINAREA0 || i==1) {	/* enough overlap: using first image */

			/* printf ("Using First image; Comarea = %f \n ", comarea); */

        	frac = correlate(x, y, p, n, img, nx, ny, ixoff, iyoff, 
                         &xoff, &yoff, hwid);

        	if (frac < MINFRAC) {             /* cross-corr failed? */
				hwid2 = 2*hwid;
            	fprintf(stderr, "-> increase search radius to %d pix\n", \
								hwid2);
            	frac = correlate(x, y, p, n, img, nx, ny, ixoff, iyoff, 
                             &xoff, &yoff, hwid2);

        		if (frac < MINFRAC) {        /* cross-corr failed again? */
					hwid4 = 2*hwid2;
            		fprintf(stderr, "-> increase search radius to %d pix\n", \
								hwid4);
            		frac = correlate(x, y, p, n, img, nx, ny, ixoff, iyoff, 
                             &xoff, &yoff, hwid4);
				}
        	}
		
			if (frac>MINFRAC) {
				corrx = corrx + xoff - ixoff ;	/* cumulate the corrections */
				corry = corry + yoff - iyoff ;
			}

			offx[i] = xoff ;
			offy[i] = yoff ;

        	/* printf("%s %f %f %f (IM=1  %d %d) \n", \
						chomp(fn[i]), xoff, yoff, frac, ixoff, iyoff); */
        	printf("%s %f %f %f (IM=0  %d %d) \n", \
						fn[i], xoff, yoff, frac, ixoff, iyoff);
        	fprintf(stderr, "%s %f %f %f (IM=0  %d %d) \n", \
						fn[i], xoff, yoff, frac, ixoff, iyoff);


		} else {	/* not enough overlap: trying with previous image */

			maxcomarea = MINAREA;
			prev = 0.;
			for (m=1; m<i; m++) {	/* loop over previous images */

				ix = abs(ixoff-offx[m]) ;	/* offsets to the n-th image */
				iy = abs(iyoff-offy[m]) ;
				if (ix<nx && iy<ny)
					comarea[m]  = (float) ((nx-ix)*(ny-iy)) / (nx*ny) ;
				else 
					comarea[m]  = 0.0 ;
				if (comarea[m]>maxcomarea) {
						prev = m;				/* image with max overlap */
						maxcomarea = comarea[m];
				}
				/* printf ("M= %d COM= %f MAX= %f PREV= %d \n", \
								m, comarea[m], maxcomarea, prev); */
			}

			if (maxcomarea<=MINAREA) {
					eprintf ("STOP: NOT ENOUGH OVERLAP for image %s \n", fn[i]);
			}

			ixoff = ixoff - offx[prev] ;	/* offsets to the previous image */
			iyoff = iyoff - offy[prev] ;

			/* printf ("IX= %d  IY= %d  COM= %f \n", ixoff, iyoff, comarea); */

			/* printf ("Using image %d ; Comarea = %f \n ", prev, comarea[prev]); */

			pimg = readfits(fn[prev], &nx, &ny, NULL, NULL);
    		for (j = 0, n = 0; j < nx * ny; j++) {    /* make pixel list */
        		if (pimg[j] > 0 && n < MAXNLIST) {
            		xprev[n] = j % nx;
            		yprev[n] = j / nx;
            		pprev[n++] = pimg[j];
        		}
    		}

        	frac = correlate(xprev, yprev, pprev, n, img, nx, ny, ixoff, iyoff, 
                         &xoff, &yoff, hwid);

        	if (frac < MINFRAC) {            /* cross-corr failed? */
			  hwid2=2*hwid;
              fprintf(stderr, "-> increase search radius to %d pix\n", \
								hwid2);
              frac = correlate(xprev, yprev, pprev, n, img, nx, ny, ixoff, iyoff, 
                             &xoff, &yoff, hwid2);

        	  if (frac < MINFRAC) {            /* cross-corr failed? */
			    hwid4=2*hwid2;
                fprintf(stderr, "-> increase search radius to %d pix\n", \
			  				hwid4);
                frac = correlate(xprev, yprev, pprev, n, img, nx, ny, ixoff, iyoff, 
                             &xoff, &yoff, hwid4);
        	    }
        	}

			if (frac>MINFRAC) {
				corrx = corrx + xoff - ixoff ;	/* cumulate the corrections */
				corry = corry + yoff - iyoff ;
			}

			offx[i] = xoff + offx[prev];
			offy[i] = yoff + offy[prev];

        	/* printf("%s %f %f %f (IM=%d: %d %d - %f %f) \n", \
				chomp(fn[i]), offx[i], offy[i], frac, prev, ixoff, iyoff, xoff, yoff); */
        	printf("%s %f %f %f (IM=%d: %d %d - %f %f) \n", \
				fn[i], offx[i], offy[i], frac, prev, ixoff, iyoff, xoff, yoff);
        	fprintf(stderr, "%s %f %f %f (IM=%d: %d %d - %f %f) \n", \
				fn[i], offx[i], offy[i], frac, prev, ixoff, iyoff, xoff, yoff);
		}

        free(img);  img = NULL;
    }

    return 0;
}

/* use WCS information from FITS header to guess image offset */
static float readwcs(char *fn, int *ixoff, int *iyoff)
{
    static int init = 0;
    static double ra0, dec0, scale0, posang0, arcsecperdeg = 3600.0;
    double ra1, dec1, scale1, posang1, ax, ay;

    if (!init) {
        init = 1;
        if (get_wcs(fn, &ra0, &dec0, &scale0, &posang0) < 0)
            eprintf("offsets: RA/DEC/SCALE missing from %s\n", fn);
        posang0 = posang0/57.2958;
        /* fprintf(stderr, "RA=%f  DEC=%f  SCL=%f  ANG=%f\n", ra0, dec0, scale0, posang0); */
        return (float)scale0;
    }


    if (get_wcs(fn, &ra1, &dec1, &scale1, &posang1) < 0)
        eprintf("offsets: RA/DEC/SCALE missing from %s\n", fn);
      	/* fprintf(stderr, "RA=%f  DEC=%f  SCL=%f  ANG=%f\n", ra1, dec1, scale1, posang1);  */

    if (scale1 / scale0 > 1.01 || scale1 / scale0 < 0.99)
        eprintf("offsets: pixel scale changed: %s %f %f\n", fn, scale0, scale1);

    if (scale1 <= 0.1 || scale1 > 10)
        eprintf("offsets: image scale [arcsec/pix] seems wrong %f\n", scale1);

    ax = arcsecperdeg * (ra1 - ra0) * cos(dec0/57.2958)/ scale0 + 0.5;
    ay = arcsecperdeg * (dec0 - dec1) / scale0 + 0.5;

    *ixoff = (int)(ax*cos(posang0) + ay*sin(posang0));
    *iyoff = (int)(ay*cos(posang0) - ax*sin(posang0));

    /* fprintf (stderr, "%s  AX=%f  AY=%f  IX=%d  IY=%d\n", fn, ax, ay, *ixoff, *iyoff);  */

    return (float)scale1;
}

/* add filename suffix ".skysub" */
static char *chomp(char *fn)
{
    char *p;
    static char buf[256];
    static char *s = ".skysub";

    if (strlen(fn) + strlen(s) + 2 > 256)
        eprintf("chomp: increase buffer size\n");

    strcpy(buf, fn);

    if ((p = strstr(buf, ".fits")) == NULL)
        return fn;

    p[5] = '\0';

    strcat(buf, ".skysub");

    return buf;
}

/* print out usage and exit */
static void usage(void)
{
    static char *usage = "\n"
    "offsets - find translation offsets between images by cross-correlation\n\n"
    "usage: offsets listfn [search_hwidth_arcsec]\n\n"
    "where listfn - list of SExtractor OBJECTS images for a dither set\n"
    "      search - half-width of search box in arcsec (def. 10)\n\n"
    "example: offsets filelist 5.0\n\n"
    "note: need search hwidth > relative error in the header RA,DEC keywords\n"
    "      scale (arcsec/pix) read from SCALE, SECPIX, CDELT, or CDELT1 key\n"
    "      offsets.c assumes that North is up and East is left\n\n";

    printf("%s", usage);
    exit(0);
}
