/* gainmap.c - produce a gain map and set bad pixels set to 0 */

/*
 * The input flatfield image is divided by its mode (normalized) to produce
 * the gain map.  A pixel is considered bad (and set to 0.0 in gain map) if 
 * the pixel gain < MINGAIN or gain > MAXGAIN or gain > NSIG sigma deviant 
 * from the local background in the gain map.
 */

#include <stdio.h>
#include <stdlib.h>
#include "irdr.h"

#define MINGAIN 0.67   /* pixels with sensitivity < MINGAIN are assumed bad */
#define MAXGAIN 1.5    /* pixels with sensitivity > MAXGAIN are assumed bad */
#define NXBLOCK 16     /* image size should be multiple of block size */
#define NYBLOCK 16
#define NSIG 5.0       /* badpix if sensitivity > NSIG sigma from local bkg */

static void usage(void);

int main(int argc, char *argv[])
{
    int i, j, k, l, nx, ny, nbad = 0, nxb = NXBLOCK, nyb = NYBLOCK;
    float scale, mode, med, sig, lo, hi, nsig = NSIG;
    float mingain = MINGAIN, maxgain = MAXGAIN;
    float *dev, *gain, *flat, *buf;

    if (argc != 3 && argc != 8)
        usage();

    if (argc == 8) {
        nsig = atof(argv[3]);
        nxb = atoi(argv[4]);
        nyb = atoi(argv[5]);
        mingain = atof(argv[6]);
        maxgain = atof(argv[7]);
    }

    flat = readfits(argv[1], &nx, &ny, &mode, NULL);

    if (mode <= 0.0)
        eprintf("%s: mode <= 0\n", argv[0]);

    scale = 1.0 / mode;
    fprintf(stderr, "  flat mode: %f\n", mode);

    gain = (float *) emalloc(nx * ny * sizeof(float));

    for (i = 0; i < nx * ny; i++) {
        gain[i] = (float)flat[i] * scale;         /* scale flatfield by mode */

        if (gain[i] < mingain || gain[i] > maxgain) {      /* bad pixel? */
            gain[i] = 0.0;
            nbad++;
        }
    }
    
    dev = (float *) emalloc(nx * ny * sizeof(float));      /* local dev map */
    buf = (float *) emalloc(nxb * nyb * sizeof(float));

    for (i = 0; i < ny; i += nyb) {
        for (j = 0; j < nx; j += nxb) {               /* foreach image block */
            int n = 0;

            for (k = i; k < i + nyb; k++)            /* foreach pix in block */
                for (l = j; l < j + nxb; l++)
                    if (gain[k*nx+l] > 0.0)                /* if good pixel */
                        buf[n++] = gain[k*nx+l];

            med = (n > 0) ? median(buf, n) : 0.0;          /* block median */
              
            for (k = i; k < i + nyb; k++)            /* foreach pix in block */
                for (l = j; l < j + nxb; l++)
                    if (gain[k*nx+l] > 0.0)               /* subtract median */
                        dev[k*nx+l] = gain[k*nx+l] - med;
                    else
                        dev[k*nx+l] = 0.0;           /* already known badpix */
        }
    }

    med = median(dev, nx * ny);
    sig = median_absdev(dev, med, nx * ny) / 0.6745;
    lo  = med - nsig * sig;
    hi  = med + nsig * sig;

    /*printf("\nMED=%f SIG=%f LO=%f HI=%f", med, sig, lo, hi);*/

    for (i = 0; i < nx * ny; i++)           /* find more badpix by local dev */
        if (dev[i] < lo || dev[i] > hi) {
            gain[i] = 0.0;
            nbad++;
        }

    writefits(argv[2], argv[1], (char*)gain, -32, nx, ny);

    put_key_float(argv[2], "DATAMODE", 1.0);
    put_key_float(argv[2], "DATASIG", sig);

    fprintf(stderr, "\n%s processed: %s\n", argv[0], argv[2]);
    fprintf(stderr, "  bad pixel count: %d\n", nbad);
    fprintf(stderr, "  sigma threshold: %f\n", nsig);
    fprintf(stderr, "  min gain threshold: %f\n", mingain);
    fprintf(stderr, "  max gain threshold: %f\n", maxgain);
    fprintf(stderr, "  block size in x: %d\n", nxb);
    fprintf(stderr, "  block size in y: %d\n", nyb);

    return 0;
}

/* print out usage and exit */
static void usage(void)
{
    static char *usage = "\n"
    "gainmap - normalize a flatfield to make a gainmap; set bad pixels to 0\n\n"
    "usage: gainmap flatfn gainfn [nsig nxblock nyblock mingain maxgain]\n\n"
    "where flatfn  - filename of input flatfield FITS image\n"
    "      gainfn  - filename of output gainmap FITS image\n"
    "      nsig    - # of stddev from local bkg to be bad pixel (default 5.0)\n"
    "      nxblock - local bkg blocks of nxblock pixels across (def. 16) *img size should be multiple\n"
    "      nyblock - local bkg blocks of nyblock pixels down (def. 16)   *img size should be multiple\n"
    "      mingain - pixels with gain < mingain are assumed bad (def. 0.7)\n"
    "      maxgain - pixels with gain > maxgain are assumed bad (def. 1.3)\n\n"
    "example: gainmap flat.fits gain.fits\n\n";

    printf("%s", usage);
    exit(0);
}
