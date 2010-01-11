/* sky.c -- driver for skysub() */

/* mask.fits is the merged objmask and gainmap, bpm.fits is the gainmap */

#include <stdio.h>
#include <stdlib.h>
#include "irdr.h"

int main(int argc, char *argv[])
{
    int nx, ny;
    float bkg, *img, *sky, *skyw, *bpm, *mask, *imgout;

    if (argc < 7)
        eprintf("Usage: %s in.fits sky.fits skyw.fits bpm.fits mask.fits"
                " out.fits [row|col]\n", argv[0]);

    img  = readfits(argv[1], &nx, &ny, &bkg, NULL);
    sky  = readfits(argv[2], &nx, &ny, NULL, NULL);
    skyw = readfits(argv[3], &nx, &ny, NULL, NULL);
    bpm  = readfits(argv[4], &nx, &ny, NULL, NULL);
    mask = readfits(argv[5], &nx, &ny, NULL, NULL);

    imgout = skysub(img, nx, ny, bkg, bpm, sky, skyw, mask, 
                    ((argc == 8) ? argv[7] : "none"));

    writefits(argv[6], argv[1], (char*)shortint(imgout, nx, ny), 16, nx, ny);

    return 0;
}
