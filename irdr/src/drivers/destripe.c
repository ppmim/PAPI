/* destripe.c -- call stripe.c to do row/column destriping */

#include <stdio.h>
#include <stdlib.h>
#include "irdr.h"

int main(int argc, char *argv[])
{
    int nx, ny;
    float bkg, *bpm, *img;

    fprintf(stderr,"\nDestripe driver......");
    if (argc != 5)
        eprintf("Usage: %s in.fits bpm.fits out.fits"
                " row|col|rowcol|colrow|none\n", argv[0]);

    img = readfits(argv[1], &nx, &ny, &bkg, NULL);
    bpm = readfits(argv[2], &nx, &ny, NULL, NULL);

    destripe(img, bpm, nx, ny, bkg, argv[4]);

    writefits(argv[3], argv[1], (char*)shortint(img, nx, ny), 16, nx, ny);

    return 0;
}
