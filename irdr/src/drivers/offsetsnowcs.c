/* offsetsnowcs.c -- use correlate.c to find translation offsets */

/* like offsets.c but don't call readwcs(), assume no knowledge of RA DEC */

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
#define MAXNFILES 99
#define MAXNLIST 99999L
#define MINFRAC  0.35            /* correlation failed if overlap < MINFRAC */

static char *fn[MAXNFILES];              /* FITS OBJECTS frames filenames */
static int n, x[MAXNLIST], y[MAXNLIST];       /* list of object pixels */
static float p[MAXNLIST];

static char *chomp(char *fn);
static void usage(void);

int main(int argc, char *argv[])
{
    int i, nx, ny, hwid, ixoff, iyoff, nfiles;
    float xoff, yoff, frac, scale, *img0, *img;

    if (argc != 3 && argc != 4)
        usage();

    if ((nfiles = readlist(argv[1], fn, NULL, NULL, NULL, MAXNFILES)) < 1)
        eprintf("%s: no valid files\n", argv[0]);

    img0 = readfits(fn[0], &nx, &ny, NULL, NULL);

    scale = atof(argv[2]);

    hwid = (argc == 4) ? (atof(argv[3]) / scale + 0.5) : (HWID / scale + 0.5);

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

    printf("%s %f %f %f\n", chomp(fn[0]), 0.0, 0.0, 1.0); /* reference frame */

    for (i = 1; i < nfiles; i++) {                        /* other frames... */
        img = readfits(fn[i], &nx, &ny, NULL, NULL);

ixoff = 0;
iyoff = 0;

        frac = correlate(x, y, p, n, img, nx, ny, ixoff, iyoff, 
                         &xoff, &yoff, hwid);

        if (frac < MINFRAC) {                          /* cross-corr failed? */
            int maxhwid = MAXNCC / 2 - 1;

            fprintf(stderr, "-> increase search radius to %d pix\n", maxhwid);

            frac = correlate(x, y, p, n, img, nx, ny, 0, 0, 
                             &xoff, &yoff, maxhwid);
        }

        printf("%s %f %f %f\n", chomp(fn[i]), xoff, yoff, frac);

        free(img);  img = NULL;
    }

    return 0;
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
    "usage: offsetsnowcs listfn imagescale [search_hwidth_arcsec]\n\n"
    "where listfn - list of SExtractor OBJECTS images for a dither set\n"
    "      search - half-width of search box in arcsec (def. 10)\n\n"
    "example: offsets filelist 5.0\n\n"
    "note: need search hwidth > relative error in the header RA,DEC keywords\n"
    "      scale (arcsec/pix) read from SCALE, SECPIX, CDELT, or CDELT1 key\n"
    "      offsets.c assumes that North is up and East is left\n\n";

    printf("%s", usage);
    exit(0);
}
