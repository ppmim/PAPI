/* imagegrid.c - use hist.c to make background map */

/* updates IRDR_BKG and IRDR_SIG header keywords */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "irdr.h"

int main(int argc, char *argv[])
{
    int nx, ny;
    float mode, sigma, *img;

    if (argc != 4)
        eprintf("Usage: %s NXblock NYblock file.fits\n", argv[0]);

    img = readfits(argv[3], &nx, &ny, &mode, &sigma);
    printf("readfits: mode %f, sigma %f\n", mode, sigma);

    mode = histcalcf(img, nx, ny, atoi(argv[1]), atoi(argv[2]), &sigma);
    printf("histcalcf: mode %f, sigma %f\n", mode, sigma);

/*
    histcalca(img, nx * ny);
*/

    return 0;
}
