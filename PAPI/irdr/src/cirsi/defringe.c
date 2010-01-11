/* defringe.c -- subtract background image to defringe */

#include <stdio.h>
#include <stdlib.h>
#include "irdr.h"

int main(int argc, char *argv[])
{
    int i, nx, ny;
    float *img, *bkg, mode, sigma;

    if (argc != 3)
        eprintf("Usage: %s img.fits bkg.fits\n", argv[0]);

    img = readfits(argv[1], &nx, &ny, NULL, NULL);
    bkg = readfits(argv[2], &nx, &ny, &mode, &sigma);

    for (i = 0; i < nx * ny; i++)
        img[i] += (mode - bkg[i]);          /* subtract out bkg structure */

    writefits(argv[1], argv[1], (char*)img, -32, nx, ny);

    return 0;
}
