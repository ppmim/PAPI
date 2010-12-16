/* cubemean.c -- given a list of FITS files, find mean image plane */

/* 
 * Procedure:
 * Read FITS data planes specified in input filelist
 * Calculate scale factors to normalize image bkg levels
 * Generate a weight map image per plane, scalar weight, or no weights
 * Calculate the robust mean (or median or sigma) plane
 *
 * Warning:
 * Need enough memory for all image planes (and optional weight planes)
 * Expects images and gain map of same dimensions
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "irdr.h"

static char *fn [MAXNPLANES];          /* FITS image file per image plane */
static float *data [MAXNPLANES];       /* pointers to image planes */
static float *wdata [MAXNPLANES];      /* pointers to weight map planes */
static float scale [MAXNPLANES];       /* scale factor per image plane */
static float weights [MAXNPLANES];     /* scalar weight per image plane */

static void usage(void);

int main(int argc, char *argv[])
{
    int i, nx, ny, nplanes;
    float bkg, sig = 0.0, avgscale = 0.0;
    float *sumwplanes, *meanplane = NULL, *gainmap = NULL;
    int noweight = 0, scalarweight = 0, mapweight = 0, doshort = 0;
    int domean = 0, domedian = 0, dosigma = 0, dooffset = 0;

    if (argc != 8 && argc != 9)
        usage();

    if (!strcmp(argv[4], "offset"))
        dooffset = 1;
    else if (!strcmp(argv[4], "scale"))
        dooffset = 0;
    else if (!strcmp(argv[4], "none"))
        dooffset = -1;
    else
        usage();

    if (!strcmp(argv[5], "mean"))
        domean = 1;
    else if (!strcmp(argv[5], "median"))
        domedian = 1;
    else if (!strcmp(argv[5], "sigma"))
        dosigma = 1;
    else
        usage();

    if (!strcmp(argv[6], "noweight"))
        noweight = 1;
    else if (!strcmp(argv[6], "scalarweight") && domean)
        scalarweight = 1;
    else if (!strcmp(argv[6], "mapweight") && domean && argc == 9) {
        mapweight = 1;
        gainmap = readfits(argv[8], &nx, &ny, NULL, NULL);
    } else
        usage();

    if (!strcmp(argv[7], "short"))
        doshort = 1;
    else if (!strcmp(argv[7], "float"))
        doshort = 0;
    else
        usage();
  
    if ((nplanes = readlist(argv[1], fn, NULL, NULL, NULL, MAXNPLANES)) < 1)
        eprintf("%s: no valid image planes", argv[0]);

    for (i = 0; i < nplanes; i++) {                  /* for each image plane */
        data[i] = readfits(fn[i], &nx, &ny, &bkg, &sig);

        if (bkg <= 0 || sig <= 0)
            eprintf("bkg %f, sig %f for %s\n", bkg, sig, fn[i]);

        avgscale += (scale[i] = bkg);                /* store scale factors */

        if (scalarweight)
            weights[i] = 1.0 / (sig * sig);
        else if (mapweight)
            wdata[i] = getwmap(fn[i], nx, ny, gainmap, sig);

        fprintf(stderr, "  %s %f %f\n", fn[i], bkg, sig);
    }

    avgscale /= nplanes;

    for (i = 0; i < nplanes; i++)
        if (dooffset==1)
            scale[i] = avgscale - scale[i];
        else if (dooffset==0)
            scale[i] = avgscale / scale[i];
        else 
            scale[i] = 1.0;

    if (dooffset==-1) 
        dooffset=0;

    if (mapweight)
        meanplane = cube_mean(data, wdata, nplanes, nx, ny, &sumwplanes, 
                              scale, dooffset);
    else if (scalarweight)
        meanplane = cube_mean_sw(data, weights, nplanes, nx, ny, 
                                 scale, dooffset);
    else if (domean && noweight)
        meanplane = cube_mean_nw(data, nplanes, nx, ny, scale, dooffset);
    else if (domedian)
        meanplane = cube_median(data, nplanes, nx, ny, scale, dooffset);
    else if (dosigma)
        meanplane = cube_sigma(data, nplanes, nx, ny, scale, dooffset);

    if (doshort)
        writefits(argv[2], fn[0], (char*)shortint(meanplane,nx,ny), 16, nx, ny);
    else
        writefits(argv[2], fn[0], (char*)meanplane, -32, nx, ny);

    put_key_int(argv[2], "NCOMBINE", nplanes);
    put_key_float(argv[2], "DATAMODE", avgscale);

    if (mapweight)
        writefits(argv[3], argv[2], (char*)sumwplanes, -32, nx, ny);

    return 0;
}

/* print out program usage and exit */
static void usage(void) 
{
    static char *usage = "\n"
    "cubemean - calculate the robust mean plane of an image stack\n\n"
    "usage: cubemean listfn outfn outwfn scale|offset mean|median|sigma\n"
    "                noweight|scalarweight|mapweight short|float [gainfn]\n\n"
    "where listfn - contains a list of FITS images to coadd\n"
    "      outfn  - filename for output coadded FITS image\n"
    "      outwfn - filename for output coadded FITS weight image\n"
    "      scale  - image normalize by multiplicative scale or zero offset\n"
    "      mean   - calculate mean, median, or stddev image plane\n"
    "      noweig - use no weights, a weight per data plane, or map weights\n"
    "      short  - output image type is either short int or floating point\n"
    "      gainfn - filename of FITS image gainmap (normalized flatfield)\n\n"
    "example: cubemean filelist coadd.fits weight.fits offset mean mapweight\n"
    "                  short gain.fits\n\n"
    "note: gainfn is required for mapweight option, noweight is required for\n"
    "      median or sigma plane calculation\n\n";

    printf("%s", usage);
    exit(0);
}
