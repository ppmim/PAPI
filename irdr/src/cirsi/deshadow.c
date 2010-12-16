/* deshadow.c -- call stripe.c to subtract column/row modes */

#include <stdio.h>
#include <stdlib.h>
#include "irdr.h"

#define NMAX 99

static char *fn[NMAX];              /* FITS image file per image */
static char *mfn[NMAX];             /* FITS objmask file per image */
static float xshift[NMAX];          /* dither offset in x direction */
static float yshift[NMAX];          /* dither offset in y direction */

int main(int argc, char *argv[])
{
    int i, nx, ny, nimg;
    float bkg, *bpm, *img, *mask;
    unsigned short *shortimg;

    if (argc != 4)
        eprintf("Usage: %s filelist bpm.fits rowcol|colrow|none\n", argv[0]);

    if ((nimg = readlist(argv[1], fn, mfn, xshift, yshift, NMAX)) < 1)
        eprintf("%s: no valid images", argv[0]);

    bpm = readfits(argv[2], &nx, &ny, NULL, NULL);

    for (i = 0; i < nimg; i++) {
        img = readfits(fn[i], &nx, &ny, &bkg, NULL);

        mask = getmask(bpm, nx, ny, mfn[i], xshift[i], yshift[i]);

        destripe(img, mask, nx, ny, bkg, argv[3]);

        shortimg = shortint(img, nx, ny);

        writefits(fn[i], fn[i], (char*)shortimg, 16, nx, ny);

        free(shortimg);  free(img);  free(mask);

        printf("%s: %s %s %s\n", argv[0], fn[i], argv[2], argv[3]);
    }

    return 0;
}
