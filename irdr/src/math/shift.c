/* shift.c -- shift an image using bilinear interpolation */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "minmax.h"
#include "shift.h"
#include "eprintf.h"

/* shift an image and weight map using bilinear interpolation */
extern float *shift_image(float *img, float *wimg, int nx, int ny, int border, 
                          float xshift, float yshift, float **wimgout)
{
    int i, j, nxout, nyout, ixoffset, iyoffset;
    float xfrac, yfrac, area[4], *imgout;

    if (border < 1 + fabs(xshift) || border < 1 + fabs(yshift))
        eprintf("shift_image: border must be larger than shift\n");

    nxout = nx + 2 * border;
    nyout = ny + 2 * border;

    ixoffset = (int)(border + xshift + 1);         /* round up shifts */
    iyoffset = (int)(border + yshift + 1);

    xfrac = (float)ixoffset - (border + xshift);   /* shift down by fraction */
    yfrac = (float)iyoffset - (border + yshift);

    area[0] = (1.0 - xfrac) * (1.0 - yfrac);
    area[1] = xfrac         * (1.0 - yfrac);
    area[2] = (1.0 - xfrac) * yfrac;               /* pixel area overlap */
    area[3] = xfrac         * yfrac;

    imgout = (float *) ecalloc(nxout * nyout, sizeof(float));
    *wimgout = (float *) ecalloc(nxout * nyout, sizeof(float));

    for (i = 0; i < ny - 1; i++) {
        float *crow = img + i * nx;                   /* current row */
        float *cwrow = wimg + i * nx;  
        float *nrow = img + (i + 1) * nx;              /* next row */
        float *nwrow = wimg + (i + 1) * nx;
        float *outrow = imgout + (i + iyoffset) * nxout + ixoffset;
        float *outwrow = *wimgout + (i + iyoffset) * nxout + ixoffset;
 
        for (j = 0; j < nx - 1; j++) {
            outwrow[j] = area[0] * cwrow[j] + area[1] * cwrow[j+1] +
                         area[2] * nwrow[j] + area[3] * nwrow[j+1] ;

            outrow[j] = (outwrow[j] <= 0) ? 0.0 : 
                        (area[0] * cwrow[j] * crow[j] + 
                         area[1] * cwrow[j+1] * crow[j+1] +
                         area[2] * nwrow[j] * nrow[j] + 
                         area[3] * nwrow[j+1] * nrow[j+1]) / outwrow[j];
        }
    }

    return imgout;
}
/* shift an image and weight map using bilinear interpolation */
extern float *new_shift_image(float *img, float *wimg, int nx, int ny, int xbelow, int xabove, int ybelow, int yabove,  float xshift, float yshift, float **wimgout)
{
    int i, j, nxout, nyout, ixoffset, iyoffset;
    float xfrac, yfrac, area[4], *imgout;

    if (xbelow>-xshift || xabove<-xshift || ybelow>-yshift || yabove<-yshift)
        eprintf("shift_image: border must be larger than shift\n");

    nxout = nx + xabove - xbelow;
    nyout = ny + yabove - ybelow;

    ixoffset = (int)(xshift + xabove + 1);         /* round up shifts */
    iyoffset = (int)(yshift + yabove + 1);

    xfrac = (float)ixoffset - (xshift + xabove);   /* shift down by fraction */
    yfrac = (float)iyoffset - (yshift + yabove);

    printf ("X: %f %d %f    Y: %f %d %f \n", 
		    xshift, ixoffset, xfrac, yshift, iyoffset, yfrac);

    area[0] = (1.0 - xfrac) * (1.0 - yfrac);
    area[1] = xfrac         * (1.0 - yfrac);
    area[2] = (1.0 - xfrac) * yfrac;               /* pixel area overlap */
    area[3] = xfrac         * yfrac;

    imgout = (float *) ecalloc(nxout * nyout, sizeof(float));
    *wimgout = (float *) ecalloc(nxout * nyout, sizeof(float));

    for (i = 0; i < ny - 1; i++) {
        float *crow = img + i * nx;                   /* current row */
        float *cwrow = wimg + i * nx;  
        float *nrow = img + (i + 1) * nx;              /* next row */
        float *nwrow = wimg + (i + 1) * nx;
        float *outrow = imgout + (i + iyoffset) * nxout + ixoffset;
        float *outwrow = *wimgout + (i + iyoffset) * nxout + ixoffset;
 
        for (j = 0; j < nx - 1; j++) {
            outwrow[j] = area[0] * cwrow[j] + area[1] * cwrow[j+1] +
                         area[2] * nwrow[j] + area[3] * nwrow[j+1] ;

            outrow[j] = (outwrow[j] <= 0) ? 0.0 : 
                        (area[0] * cwrow[j] * crow[j] + 
                         area[1] * cwrow[j+1] * crow[j+1] +
                         area[2] * nwrow[j] * nrow[j] + 
                         area[3] * nwrow[j+1] * nrow[j+1]) / outwrow[j];
        }
    }

    return imgout;
}

/* reverse the effects of shift_image (unshift, remove border) */
extern float *unshift_image(float *img, float *wimg, int nx, int ny, int border,
                            float xshift, float yshift, float **wimgout)
{
    int i, j, nout, nxunshift, nyunshift, n = 0;
    float *unshift, *wunshift, *imgout;

    unshift = 
        shift_image(img, wimg, nx, ny, border, -xshift, -yshift, &wunshift);

    nxunshift = nx + 2 * border;
    nyunshift = ny + 2 * border;

    nout = (nxunshift - 4 * border) * (nyunshift - 4 * border);

    imgout = (float *) emalloc(nout * sizeof(float));
    *wimgout = (float *) emalloc(nout * sizeof(float));

    for (i = 2 * border; i < nyunshift - 2 * border; i++) {
        float *row = unshift + i * nxunshift;
        float *wrow = wunshift + i * nxunshift;

        for (j = 2 * border; j < nxunshift - 2 * border; j++) {
            imgout[n] = row[j];
            (*wimgout)[n++] = wrow[j];
        }
    }

    free(unshift);  free(wunshift);

    return imgout;
}

/* get_border: calculate border size from max image shift in list of shifts */
extern int get_border(float *xshift, float *yshift, int n)
{
    int i, imax, pad = 4;
    float dx, dy, maxv, maxshift = 0.0;

    for (i = 0; i < n; i++) {
        dx = (float) fabs((double)xshift[i]);
        dy = (float) fabs((double)yshift[i]);

        if ((maxv = max(dx, dy)) > maxshift)
            maxshift = maxv;
    }

    maxshift += 2;                /* padding to avoid edge in shift_image() */

    imax = maxshift + 1.0;                        /* round up */

    imax += (pad - (imax % pad));   /* enlarge border to padding multiple */

    return imax;
}

/* new_get_border: calculate border sizes (xabove, xbelow, yabove, ybelow)
 * from max image shift in list of shifts */
extern int new_get_border(float *xshift, float *yshift, int n, int *xbelow, int *xabove, int *ybelow, int *yabove)
{
    int i;
    float x1=0.0, x2=0.0, y1=0.0, y2=0.0;

    for (i = 0; i < n; i++) {
	x1 = min(x1, xshift[i]);
	x2 = max(x2, xshift[i]);
	y1 = min(y1, yshift[i]);
	y2 = max(y2, yshift[i]);
    }

    *xbelow = x1 - 2.0;
    *xabove = x2 + 2.0;
    *ybelow = y1 - 2.0;
    *yabove = y2 + 2.0;

     printf (" PASSED: %d %d %d %d \n", *xbelow, *xabove, *ybelow, *yabove);

    return 1;
}
