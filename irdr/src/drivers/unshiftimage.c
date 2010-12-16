/* unshiftimage.c -- test shift.c */

#include <stdio.h>
#include <stdlib.h>
#include "irdr.h"

int main(int argc, char *argv[])
{
    int nx, ny, border;
    float *img, *wimg, *imgout, *wimgout, xshift, yshift;

    if (argc != 8)
        eprintf("Usage: %s shift.fits wshift.fits unshift.fits wunshift.fits "
                "border xshift yshift\n", argv[0]);

    img = readfits(argv[1], &nx, &ny, NULL, NULL);
    wimg = readfits(argv[2], &nx, &ny, NULL, NULL);
    border = atoi(argv[5]);

    xshift = (float)atof(argv[6]);
    yshift = (float)atof(argv[7]);

    imgout = unshift_image(img, wimg, nx, ny, border, xshift, yshift, &wimgout);

    nx = nx - 2 * border;
    ny = ny - 2 * border;

    writefits(argv[3], argv[1], (char*)imgout, -32, nx, ny);
    writefits(argv[4], argv[2], (char*)wimgout, -32, nx, ny);

    return 0;
}
