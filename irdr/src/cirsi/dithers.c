/* dithers.c -- group a list of FITS files into dither sets */

/*
 * a dither set is defined as consecutive frames with RA/DEC position offsets
 * smaller than 1/4 the chip size.  reads RA/DEC FITS header keywords.  should
 * be run on loop 1 raw data (ie, RA/DEC keywords refer to telescope position,
 * not chip center).  writes ASCII list of dither sets called "dithers".
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "irdr.h"

#define THRESH 256   /* offset tolerance to be in same dither set [pixels] */

int main(int argc, char *argv[])
{
    int i, run_beg = 0, run_curr = 0, run_prev = 0;
    double ra_prev = 0.0, dec_prev = 0.0;
    double ra, dec, dra, ddec, scale, posang;
    char *outfn = "dithers";
    FILE *fp;

    if (argc < 2)
        eprintf("Usage: %s irx*001.fits\n", argv[0]);

    if ((fp = fopen(outfn, "w")) == NULL)
        eprintf("%s: open %s for write failed\n", argv[0], outfn); 

    for (i = 1; i < argc; i++) {
        if (get_wcs(argv[i], &ra, &dec, &scale, &posang) < 0)
            eprintf("%s: get_wcs failed\n", argv[0]);

        if (get_key_int(argv[i], "NRUN", &run_curr) < 0)
            if (get_key_int(argv[i], "FILE", &run_curr) < 0)
                eprintf("%s: failed reading NRUN/FILE keyword\n", argv[0]);

        if (i == 1) {
            run_beg = run_curr;
            ra_prev = ra;
            dec_prev = dec;
        }

        dra = 3600.0 * fabs(ra - ra_prev) / scale;           /* in pixels */
        ddec = 3600.0 * fabs(dec - dec_prev) / scale;

        if (dra > THRESH || ddec > THRESH) {
            fprintf(fp, "%d %d\n", run_beg, run_prev);
            run_beg = run_curr;
        }

        run_prev = run_curr;
        dec_prev = dec;
        ra_prev = ra;
    }

    fprintf(fp, "%d %d\n", run_beg, run_prev);

    fclose(fp);

    return 0;
}
