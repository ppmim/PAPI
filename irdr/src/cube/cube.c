/* cube.c -- find median, mean, or stdev image plane of an image cube */

#include <stdio.h>
#include <stdlib.h>
#include "kselect.h"
#include "eprintf.h"
#include "median.h"
#include "mean.h"
#include "cube.h"

#define NSIG 2.0   /* clipping threshold */

/* influye bastante en el modo skyfilter_on_off/off_on 
#define NSIG 2.0
*/

/* 
 * cube_median_cl: find clipped median image plane of image cube 
 */

extern float *
cube_median_cl(float *planes[MAXNPLANES], int np, int nx, int ny, float *scale, 
            int offset)
{
    int i, j, is_even = !(np & 1);
    static float buf[MAXNPLANES];      /* values of a pixel in all planes */
    float *medplane;
    float med, sig, lcut, hcut,a=0;
    int k,p;
    int nsig = (int) NSIG;
    static float buf2[MAXNPLANES];

    medplane = (float *) emalloc(nx * ny * sizeof(float));

    if (offset) {                          /* zero offset to normalize */
        for (i = 0; i < nx * ny; i++) {    /* median combine cube to plane */
            for (j = 0; j < np; j++)
                buf[j] = *(planes[j] + i) + scale[j];
                
            med = median(buf, np);                              /* clipping */
            sig = median_absdev(buf, med, np) / 0.6745;

            lcut = med - nsig * sig;
            hcut = med + nsig * sig;
            p=0;
            for (k=0;k<np; k++)
                if ((a = buf[k]) >= lcut && a <= hcut) buf2[p++]= a;
                    
            medplane[i] = kselect(buf2, p, p/2 - (!(p & 1)));
        }
    } else {
        for (i = 0; i < nx * ny; i++) {   /* mult. scale to normalize */
            for (j = 0; j < np; j++)
                buf[j] = *(planes[j] + i) * scale[j];

            medplane[i] = kselect(buf, np, np/2 - is_even);
        }
    }

    return medplane;
}

/* 
 * cube_median: find median image plane of image cube 
 */

extern float *
cube_median(float *planes[MAXNPLANES], int np, int nx, int ny, float *scale, 
            int offset)
{
    int i, j, is_even = !(np & 1);
    static float buf[MAXNPLANES];      /* values of a pixel in all planes */
    float *medplane;

    medplane = (float *) emalloc(nx * ny * sizeof(float));

    if (offset) {                          /* zero offset to normalize */
        for (i = 0; i < nx * ny; i++) {    /* median combine cube to plane */
            for (j = 0; j < np; j++)
                buf[j] = *(planes[j] + i) + scale[j];

            medplane[i] = kselect(buf, np, np/2 - is_even);
        }
    } else {
        for (i = 0; i < nx * ny; i++) {   /* mult. scale to normalize */
            for (j = 0; j < np; j++)
                buf[j] = *(planes[j] + i) * scale[j];

            medplane[i] = kselect(buf, np, np/2 - is_even);
        }
    }

    return medplane;
}


/* 
 * cube_mean_nw: find clipped mean image plane of image cube (no weights) 
 */

extern float *
cube_mean_nw(float *planes[MAXNPLANES], int np, int nx, int ny, float *scale, 
             int offset)
{
    int i, j;
    static float buf[MAXNPLANES];       /* values of a pixel in all planes */
    float *meanplane;

    meanplane = (float *) emalloc(nx * ny * sizeof(float));

    if (offset) {                             /* zero offset to normalize */
        for (i = 0; i < nx * ny; i++) {       /* mean combine cube to plane */
            for (j = 0; j < np; j++)
                buf[j] = *(planes[j] + i) + scale[j];

            meanplane[i] = mean_nw(buf, np, NSIG);
        }
    } else {
        for (i = 0; i < nx * ny; i++) {      /* mult. scale to normalize */
            for (j = 0; j < np; j++)
                buf[j] = *(planes[j] + i) * scale[j];
                                     
            meanplane[i] = mean_nw(buf, np, NSIG);
        }
    }

    return meanplane;
}


/* 
 * cube_mean_sw: find clipped, mean image plane of cube (scalar weights) 
 */

extern float *
cube_mean_sw(float *planes[MAXNPLANES], float *weights, int np, int nx, int ny,
             float *scale, int offset)
{
    int i, j;
    static float buf[MAXNPLANES];        /* values of a pixel in all planes */
    float wsum, *meanplane;

    meanplane = (float *) emalloc(nx * ny * sizeof(float));

    if (offset) {                            /* zero offset to normalize */
        for (i = 0; i < nx * ny; i++) {      /* mean combine cube to plane */
            for (j = 0; j < np; j++)
                buf[j] = *(planes[j] + i) + scale[j];

            meanplane[i] = mean(buf, weights, np, NSIG, &wsum);
        }
    } else {
        for (i = 0; i < nx * ny; i++) {      /* mult. scale to normalize */
            for (j = 0; j < np; j++)
                buf[j] = *(planes[j] + i) * scale[j];

            meanplane[i] = mean(buf, weights, np, NSIG, &wsum);
        }
    }

    return meanplane;
}


/* 
 * cube_mean: find clipped, mean image plane of cube (full map weights) 
 */

extern float *
cube_mean(float *planes[MAXNPLANES], float *wplanes[MAXNPLANES], int np, 
          int nx, int ny, float **wplanesum, float *scale, int offset)
{
    int i, j, nval;
    static float buf[MAXNPLANES], wbuf[MAXNPLANES];
    float *sumwplanes, *meanplane, wval;

    meanplane = (float *) emalloc(nx * ny * sizeof(float));
    sumwplanes = (float *) emalloc(nx * ny * sizeof(float));

    if (offset) {                            /* zero offset to normalize */
        for (i = 0; i < nx * ny; i++) {      /* mean combine cube to plane */
            nval = 0;

            for (j = 0; j < np; j++)
                if ((wval = *(wplanes[j] + i)) > 0.0) {     /* if not masked */
                    buf[nval] = *(planes[j] + i) + scale[j];
                    wbuf[nval++] = wval;
                }

            meanplane[i] = mean(buf, wbuf, nval, NSIG, &sumwplanes[i]);
        }
    } else {                                /* mult. scale to normalize */
        for (i = 0; i < nx * ny; i++) {
            nval = 0;

            for (j = 0; j < np; j++)
                if ((wval = *(wplanes[j] + i)) > 0.0) {     /* if not masked */
                    buf[nval] = *(planes[j] + i) * scale[j];
                    wbuf[nval++] = wval;
                }

            meanplane[i] = mean(buf, wbuf, nval, NSIG, &sumwplanes[i]);
        }
    }

    *wplanesum = sumwplanes;

    return meanplane;
}


/* 
 * cube_sigma: find robust standard deviation plane of cube 
 */

extern float *
cube_sigma(float *planes[MAXNPLANES], int np, int nx, int ny, float *scale, 
           int offset)
{
    int i, j;
    static float buf[MAXNPLANES];       /* values of a pixel in all planes */
    float *sigplane;

    sigplane = (float *) emalloc(nx * ny * sizeof(float));

    if (offset) {                                /* zero offset to normalize */
        for (i = 0; i < nx * ny; i++) {            /* calculate sigma plane */
            for (j = 0; j < np; j++)
                buf[j] = *(planes[j] + i) + scale[j];

            sigplane[i] = median_absdev(buf, median(buf, np), np) / 0.6745;
        }
    } else {
        for (i = 0; i < nx * ny; i++) {          /* mult. scale to normalize */
            for (j = 0; j < np; j++)
                buf[j] = *(planes[j] + i) * scale[j];

            sigplane[i] = median_absdev(buf, median(buf, np), np) / 0.6745;
        }
    }

    return sigplane;
}
