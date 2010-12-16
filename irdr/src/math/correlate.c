/* correlate.c -- cross-correlate images to find translation to align */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "eprintf.h"
#include "parabola.h"
#include "correlate.h"

static float ccimg[MAXNCC][MAXNCC];           /* cross-correlation image */
static float xprof[MAXNCC], yprof[MAXNCC];          /* 1-D profiles */

/*
 * cross-correlate the reference frame list of object pixels (input x, y, p)
 * with an image (input img) and fit a parabola to the CC peak to measure the 
 * image offset.  (ixoff, iyoff) is the input offset guess, and (xoff, yoff) 
 * are the output measured image offset.
 */

extern float correlate(int *x, int *y, float *p, int nlist, float *img, int nx,
              int ny, int ixoff, int iyoff, float *xoff, float *yoff, int hwid)
{
    int i, j, k, l, ncc, ix0 = 0, iy0 = 0, count = 0;
    float maxv;

    if ((ncc = 2 * hwid + 1) > MAXNCC)
        eprintf("correlate: increase MAXNCC: %d\n", ncc);

    if (nlist < 1)
        return 0.0;

    for (i = 0; i < ncc; i++)         /* initialize cross-correlation image */
        for (j = 0; j < ncc; j++)
            ccimg[i][j] = 0.0;

    for (i = 0; i < nlist; i++) {
        int ipix = (y[i] + iyoff) * nx + x[i] + ixoff;
        
        for (k = -hwid; k <= hwid; k++) {           /* cross-correlate */
            int m = ipix + k * nx;

            for (l = -hwid; l <= hwid; l++) {
                int n = m + l;

                if (n >= 0 && n < nx * ny)
                    ccimg[k+hwid][l+hwid] += p[i] * img[n];
            }
        }
    }

    for (i = 0, maxv = ccimg[0][0]; i < ncc; i++)        /* locate CC peak */
        for (j = 0; j < ncc; j++)
            if (ccimg[i][j] > maxv)
                maxv = ccimg[iy0 = i][ix0 = j];

    if (ix0 < 1 || ix0 > ncc-2 || iy0 < 1 || iy0 > ncc-2) {
        fprintf(stderr, "correlate: CC peak near edge\n");
        return 0.0;
    }

    for (i = 0; i < ncc; i++) {                     /* extract 1-D profiles */
        yprof[i] = ccimg[i][ix0];
        xprof[i] = ccimg[iy0][i];
    }

    *xoff = fit_parabola(xprof, ncc) - hwid + ixoff;    /* fit parabolas */
    *yoff = fit_parabola(yprof, ncc) - hwid + iyoff;    /* to get center */

    ixoff = ixoff + (ix0 - hwid);
    iyoff = iyoff + (iy0 - hwid);               /* calculate match fraction */

    for (i = 0; i < nlist; i++) {
        int ipix = (y[i] + iyoff) * nx + x[i] + ixoff;

        if (ipix > 0 && ipix < nx * ny && img[ipix] > 0)
            count++;
    }

    return (float)count / nlist;                 /* return match fraction */
}
