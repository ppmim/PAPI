/* normalize.c -- offset or scale image to specified input level */

/*
 * method argument is either offset or scale.  norm is the signal level to
 * scale to.  overwrites input image with scaled image.  reads DATAMODE
 * from header.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "irdr.h"

int main(int argc, char *argv[])
{
    int i, nx, ny;
    float *img, *out, norm, mode, sigma, scale = 0;

    if (argc != 4)
        eprintf("Usage: %s file.fits norm method\n", argv[0]);

    img = readfits(argv[1], &nx, &ny, NULL, NULL);

    norm = (float) atof(argv[2]);

    out = (float *) emalloc(nx * ny * sizeof(float));

    if (get_key_float(argv[1], "DATAMODE", &mode) < 0)
        eprintf("%s: failed reading DATAMODE\n", argv[0]);

    if (mode <= 0.0)
        eprintf("%s: mode <= 0\n", argv[0]);

    if (!strcmp(argv[3], "scale")) {
        scale = norm / mode;

        for (i = 0; i < nx * ny; i++)
            out[i] = scale * img[i];

    } else if (!strcmp(argv[3], "offset")) {
        scale = norm - mode;

        for (i = 0; i < nx * ny; i++)
            out[i] = scale + img[i];

    } else {
        eprintf("method should be offset or scale\n");
    }

    writefits(argv[1], argv[1], (char*)out, -32, nx, ny);

    put_key_float(argv[1], "DATAMODE", norm);

    if (!strcmp(argv[3], "scale"))                           /* fix DATASIG */
        if (!get_key_float(argv[1], "DATASIG", &sigma))
            put_key_float(argv[1], "DATASIG", scale * sigma);

    return 0;
}
