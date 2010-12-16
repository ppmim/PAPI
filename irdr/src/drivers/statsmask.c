/* stats_mask.c -- calculate DATAMODE and DATASIG and write to hdr */

/* uses input weight map (with bad pixels set to 0.0) and object mask */

#include <stdio.h>
#include <stdlib.h>
#include "irdr.h"

#define MIN_WEIGHT 0.0

int main(int argc, char *argv[])
{
    int i, n, nx, ny;
    float *img, *buf, *wmap, *mask, mode, sigma;

    if (argc != 4)
        eprintf("Usage: %s file.fits wmap.fits mask.fits\n", argv[0]);

    img  = readfits(argv[1], &nx, &ny, NULL, NULL);
    wmap = readfits(argv[2], &nx, &ny, NULL, NULL);
    mask = readfits(argv[3], &nx, &ny, NULL, NULL);

    buf = (float *) emalloc(nx * ny * sizeof(float));

    for (i = 0, n = 0; i < nx * ny; i++)
        if (wmap[i] > MIN_WEIGHT && mask[i] == 0.0)
            buf[n++] = img[i];
        
    mode = histcalcf(buf, n, 1, n, 1, &sigma);

    printf("%s, BKG: %f, SIG: %f\n", argv[1], mode, sigma);

    put_key_float(argv[1], "DATAMODE", mode);
    put_key_float(argv[1], "DATASIG", sigma);

    return 0;
}
