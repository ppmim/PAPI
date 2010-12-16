/* floatimage.c - read a u_short image and write float image */

#include <stdio.h>
#include <stdlib.h>
#include "irdr.h"

int main(int argc, char *argv[])
{
    int nx, ny;
    float *img;

    if (argc != 3)
        eprintf("Usage: %s input.fits output.fits\n", argv[0]);

    img = readfits(argv[1], &nx, &ny, NULL, NULL);

    writefits(argv[2], argv[1], (char *)img, -32, nx, ny);

    return 0;
}
