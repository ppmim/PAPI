/* printwcs.c -- print WCS info to stdout */

#include <stdio.h>
#include <stdlib.h>
#include "irdr.h"

int main(int argc, char *argv[])
{
    int nx, ny;
    double ra, dec, scale, secpix1, secpix2, posang;

    if (argc != 2)
        eprintf("Usage: %s file.fits\n", argv[0]);

    if (get_wcs(argv[1], &ra, &dec, &scale, &posang) < 0)
        eprintf("%s: get_wcs failed\n", argv[0]);

    if (get_key_int(argv[1], "NAXIS1", &nx) < 0)
        eprintf("%s: failed reading NAXIS1\n", argv[0]);

    if (get_key_int(argv[1], "NAXIS2", &ny) < 0)
        eprintf("%s: failed reading NAXIS2\n", argv[0]);

    /*printf("RA=%f  DEC=%f  SCALE=%f  NX=%d  NY=%d  ANG=%f\n", ra, dec, scale, nx, ny, posang);*/
    printf("%f  %f  %f  %d  %d  %f\n", ra, dec, scale, nx, ny, posang);

    return 0;
}
