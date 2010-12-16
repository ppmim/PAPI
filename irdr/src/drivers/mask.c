/* mask.c - produce a weight map with object pixels masked */

#include <stdio.h>
#include <stdlib.h>
#include "irdr.h"

static char *usage = 
    "Usage: %s img.fits gain.fits mask.fits obj.fits xshift yshift";

int main(int argc, char *argv[])
{
    int nx, ny;
    float *gain, *wmap, *img, *mask;
    float xoff, yoff, mode, sigma;

    if (argc != 7)
        eprintf(usage, argv[0]);

    img = readfits(argv[1], &nx, &ny, &mode, &sigma);

    gain = readfits(argv[2], &nx, &ny, NULL, NULL);

    wmap = getwmap(argv[1], nx, ny, gain, sigma);

    writefits("wmap.fits", argv[1], (char*)wmap, -32, nx, ny);

    xoff = atof(argv[5]);
    yoff = atof(argv[6]);

    mask = getmask(wmap, nx, ny, argv[4], xoff, yoff);

    writefits(argv[3], argv[1], (char*)mask, -32, nx, ny);

    return 0;
}
