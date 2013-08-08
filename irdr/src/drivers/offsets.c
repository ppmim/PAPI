/* offsets.c -- use correlate.c to find translation offsets between images */

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
 *
 * IMPORTANT: The output of this program is dumped to an output file used as input for the PAPI,
 *            thus any 'printf' could cause malfunction on the next PAPI step.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "correlate.h"
#include "eprintf.h"
#include "fitsIO.h"
#include "listIO.h"

#define HWID 10.0     /* default half-width of cross-corr search box, arcsec */
#define MAXNFILES 999
#define MAXNLIST 999999L
#define MINFRAC  0.1      /* 0.15; correlation failed if overlap < MINFRAC */

static char *fn[MAXNFILES];              /* FITS OBJECTS frames filenames */
static int n, x[MAXNLIST], y[MAXNLIST];       /* list of object pixels */
static float p[MAXNLIST];

static float readwcs(char *fn, int *ixoff, int *iyoff);
static char *chomp(char *fn);
static void usage(void);

int main(int argc, char *argv[])
{
    int i, nx, ny, hwid, ixoff, iyoff, nfiles;
    float xoff, yoff, frac, scale, *img0, *img, corrx, corry;

    if (argc != 2 && argc != 3)
        usage();

    if ((nfiles = readlist(argv[1], fn, NULL, NULL, NULL, MAXNFILES)) < 1)
        eprintf("%s: no valid files\n", argv[0]);

    img0 = readfits(fn[0], &nx, &ny, NULL, NULL);

    scale = readwcs(fn[0], NULL, NULL);              /* arcsec per pixel */
  
    /*fprintf(stderr, "\nDEBUGGIN ....\n->offsets_SCALE=%f", scale);*/

    /* Initial half-width of cross-corr search box, in pixels, but the argv value is given in arcsec*/
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
        eprintf("%s: n >= MAXNLIST, too much object pixel found. Please, increase detection threshold\n", argv[0]);

    if (n < 1)
        eprintf("%s: found no object pixels in %s\n", argv[0], fn[0]);

    printf("%s %f %f %f\n", fn[0], 0.0, 0.0, 1.0); /* reference frame */

    corrx = 0.;
    corry = 0.;

    for (i = 1; i < nfiles; i++) 
    {   
        /* other frames... */
        img = readfits(fn[i], &nx, &ny, NULL, NULL);

        (void) readwcs(fn[i], &ixoff, &iyoff);           /* get offset guess (in pixels)*/
    
        /* ixoff,iyoff --> estimated offset */
	    ixoff = ixoff+corrx;	/* use the previous corrections */
	    iyoff = iyoff+corry;
    
        frac = correlate(x, y, p, n, img, nx, ny, ixoff, iyoff, 
                            &xoff, &yoff, hwid);
    
	    fprintf(stderr, "\n[iter_0] corr_frac=%f, n=%d, nx=%d, ny=%d, ixoff=%d, iyoff=%d, xoff=%f, yoff=%f, hwid=%d\n",
            	   frac, n, nx, ny, ixoff, iyoff, xoff, yoff, hwid);
        
        /*fprintf(stderr, "-->First correlation overlap computed is : %f", frac);*/

        if (frac < MINFRAC) {                          /* cross-corr failed? */
            int maxhwid = MAXNCC / 2 - 1;
    
            fprintf(stderr, "-> Bad correlation (%f) overlap, increasing search radius to %d pix\n", frac, maxhwid);
    
            frac = correlate(x, y, p, n, img, nx, ny, 0, 0, /* no estimation are given this time */
                                &xoff, &yoff, maxhwid);
            
            fprintf(stderr, "\n[iter_1] corr_frac=%f, n=%d, nx=%d, ny=%d, ixoff=%d, iyoff=%d, xoff=%f, yoff=%f, hwid=%d\n",
                	   frac, n, nx, ny, ixoff, iyoff, xoff, yoff, maxhwid);
        
            /*fprintf(stderr, "-->Second correlation overlap computed is : %f", frac);*/

            if (frac < MINFRAC) {
                fprintf(stderr, "-> Still Bad correlation (%f) overlap; cannot find translation offsets overlap for this image\n", frac);
                
            }    
        }
    
        corrx = corrx + xoff - ixoff ;	/* cumulate the corrections , to take into account with the next frame */
        corry = corry + yoff - iyoff ;
        fprintf(stderr, "\nCORRX = %f  CORRY=%f \n",corrx, corry);

        /*printf("%s %f %f %f (%d %d) \n", chomp(fn[i]), xoff, yoff, frac, ixoff, iyoff);*/
        /* print out the results (in pixels) to the std output, which will be dumped with ">" to a text file used in the next step in the pipeline */
        printf("%s %f %f %f (%d %d) \n", fn[i], xoff, yoff, frac, ixoff, iyoff);
    
        free(img);  img = NULL;
    }

    return 0;
}

/* use WCS information from FITS header to guess image offset (in pixels) */
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
        /*printf("RA=%f  DEC=%f  SCL=%f  ANG=%f\n", ra0, dec0, scale0, posang0);*/ 
        return (float)scale0;
    }


    if (get_wcs(fn, &ra1, &dec1, &scale1, &posang1) < 0)
        eprintf("offsets: RA/DEC/SCALE missing from %s\n", fn);
      /* printf("RA=%f  DEC=%f  SCL=%f  ANG=%f\n", ra1, dec1, scale1, posang1);  */

    if (scale1 / scale0 > 1.01 || scale1 / scale0 < 0.99)
        eprintf("offsets: pixel scale changed: %s %f %f\n", fn, scale0, scale1);

    if (scale1 < 0.1 || scale1 > 10)
        eprintf("offsets: image scale [arcsec/pix] seems wrong %f\n", scale1);

    ax = arcsecperdeg * (ra1 - ra0) * cos(dec0/57.296)/ scale0 + 0.5;
    ay = arcsecperdeg * (dec0 - dec1) / scale0 + 0.5;

    *ixoff = (int)(ax*cos(posang0) + ay*sin(posang0));
    *iyoff = (int)(ay*cos(posang0) - ax*sin(posang0));

    /* printf ("%s  AX=%f  AY=%f  IX=%d  IY=%d\n", fn, ax, ay, *ixoff, *iyoff); */

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
