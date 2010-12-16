/* dilate.c -- grow an object mask by dilation */

/* 
 * input image is overwritten with dilated image.  mask region to be expanded
 * is indicated by pixel values > 0 (eg, a SExtractor OBJECTS image)
 */

#include <stdio.h>
#include <stdlib.h>
#include "irdr.h"

#define SCALE 0.5    /* default mult. scale factor to expand mask regions */

static void usage(void);

int main(int argc, char *argv[])
{
    int nx, ny;
    float *mask, scale;

    if (argc != 2 && argc != 3)
        usage();

    scale = (argc == 3) ? atof(argv[2]) : SCALE;

    mask = readfits(argv[1], &nx, &ny, NULL, NULL);

    dilate(mask, nx, ny, scale);

    writefits(argv[1], argv[1], (char*)shortint(mask, nx, ny), 16, nx, ny);

    return 0;
}

/* print out usage and exit */
static void usage(void)
{
    static char *usage = "\n"
    "dilate - dilate an object mask by a multiplicative expansion factor\n\n"
    "usage: dilate objmaskfn [scale]\n\n"
    "where objmaskfn - filename of FITS object mask (SExtractor OBJECTS img)\n"
    "      scale     - mult. scale factor to expand object regions; default\n"
    "                  is 0.5 (ie, make 50%% larger)\n\n"
    "example: dilate mask.fits\n\n";

    printf("%s", usage);
    exit(0);
}
