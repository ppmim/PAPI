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

    if (offset) {                          /* zero offset to normalize, it is much more "safe" than a multiplicative normalization  */
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

    if (offset) {                          /* zero offset to normalize; it is much more "safe" than a multiplicative normalization */
        for (i = 0; i < nx * ny; i++) {    /* median combine cube to plane */
            for (j = 0; j < np; j++)
            {
                buf[j] = *(planes[j] + i) + scale[j];
            }
            
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
 * cube_sum: find out the sum image plane of image cube 
 * Added by jmiguel 16-Feb-2015;
 */

extern float *
cube_sum(float *planes[MAXNPLANES], int np, int nx, int ny)
{
    int i, j, is_even = !(np & 1);
    static float buf[MAXNPLANES];      /* values of a pixel in all planes */
    float *sumplane;

    sumplane = (float *) emalloc(nx * ny * sizeof(float));

    for (i = 0; i < nx * ny; i++) {    /* sum combine cube to plane */
        sumplane[i] = 0;    
        for (j = 0; j < np; j++)
            sumplane[i] += (*(planes[j] + i));
    }

    return sumplane;
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

    if (offset) {                             /* zero offset to normalize , it is much more "safe" than a multiplicative normalization */
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

    if (offset) {                            /* zero offset to normalize, it is much more "safe" than a multiplicative normalization  */
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

    if (offset) {                            /* zero offset to normalize, it is much more "safe" than a multiplicative normalization  */
        for (i = 0; i < nx * ny; i++) {      /* mean combine cube to plane */
            nval = 0;

            for (j = 0; j < np; j++)
                if ((wval = *(wplanes[j] + i)) > 0.0) {     /* if not masked ==> neither object nor bad pixel */
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

    if (offset) {                                /* zero offset to normalize , it is much more "safe" than a multiplicative normalization */
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

/* 
 * cube_mean_min: find mean value of the N-smallest values of the sorted plane 
 * of cube (no weights). It is thought for sky background of crowded fields.
 * Implemented by jmiguel@iaa.es on 2011-Oct-20
 */

extern float *
cube_mean_min(float *planes[MAXNPLANES], int np, int nx, int ny, float *scale, 
            int offset, int N)

{
    int i, j;
    static float buf[MAXNPLANES];      /* values of a pixel in all planes */
    float *meanplane, mean=0;
    int n=0;
    
    meanplane = (float *) emalloc(nx * ny * sizeof(float));

    if (offset) {                          /* zero offset to normalize, it is much more "safe" than a multiplicative normalization  */
        for (i = 0; i < nx * ny; i++) {    /* median combine cube to plane */
            mean = 0;
            for (j = 0; j < np; j++)
                buf[j] = *(planes[j] + i) + scale[j];

            for (n=0; n< N; n++)            /* take the n-smallest values */
                mean+= kselect(buf, np, n );   
                            
            meanplane[i] = mean/N;
        }
    } else {
        for (i = 0; i < nx * ny; i++) {   /* mult. scale to normalize */
            mean = 0;
            for (j = 0; j < np; j++)      /* take the n-smallest values */
                buf[j] = *(planes[j] + i) * scale[j];

            for (n=0; n< N; n++)
                mean+= kselect(buf, np, n );   
                            
            meanplane[i] = mean/N;

        }
    }

    return meanplane;
}

/* 
 * cube_mean_min: find mean value of the N-biggest values of the sorted plane 
 * of cube (no weights). It is thought to detect stars in to well aligned
 * dithered images due to field distortion.
 * Implemented by jmiguel@iaa.es on 2015-March-12.
 */

extern float *
cube_mean_max(float *planes[MAXNPLANES], int np, int nx, int ny, float *scale, 
            int offset, int N)

{
    int i, j;
    static float buf[MAXNPLANES];      /* values of a pixel in all planes */
    float *meanplane, mean = 0;
    int n = 0, st = 0;
    int is_even = !(np & 1);
    
    meanplane = (float *) emalloc(nx * ny * sizeof(float));

    if (offset) {                          /* zero offset to normalize, it is much more "safe" than a multiplicative normalization  */
        for (i = 0; i < nx * ny; i++) {    /* median combine cube to plane */
            mean = 0;
            for (j = 0; j < np; j++)
                buf[j] = *(planes[j] + i) + scale[j];
            if (!is_even) st = N + 1;
            else st = N;
            for (n=N; n < np; n++)            /* take the n-biggest values */
                mean+= kselect(buf, np, n );
                            
            meanplane[i] = mean/N;
        }
    } else {
        for (i = 0; i < nx * ny; i++) {   /* mult. scale to normalize */
            mean = 0;
            for (j = 0; j < np; j++)      /* take the n-biggest values */
                buf[j] = *(planes[j] + i) * scale[j];

            if (!is_even) st = N + 1;
            else st = N;
            for (n=st; n < np; n++)            /* take the n-biggest values */
                mean+= kselect(buf, np, n );
                            
            meanplane[i] = mean/N;

        }
    }

    return meanplane;
}



/* 
 * cube_mean_min_w: find mean value of the N-smallest values of the sorted plane 
 * of cube (full map weights) 
 * Implemented by jmiguel@iaa.es on 2011-Oct-25
 */

extern float *
cube_mean_min_w(float *planes[MAXNPLANES], float *wplanes[MAXNPLANES], int np, 
          int nx, int ny, float **wplanesum, float *scale, int offset, int N)
{

    int i, j, nval = 0, n = 0;
    static float buf[MAXNPLANES], wbuf[MAXNPLANES]; /* values of a pixel in all planes */
    static float buf2[MAXNPLANES]; /* wbuf2[MAXNPLANES];  values after masking */
    float mean = 0;
    float *sumwplanes, *meanplane, wval;

    meanplane = (float *) emalloc(nx * ny * sizeof(float));
    sumwplanes = (float *) emalloc(nx * ny * sizeof(float));
    

    if (offset) {                          /* zero offset to normalize, it is much more "safe" than a multiplicative normalization  */
        for (i = 0; i < nx * ny; i++) {    /* median combine cube to plane */
            mean = 0;
            nval = 0;
            for (j = 0; j < np; j++)
                if ((wval = *(wplanes[j] + i)) > 0.0) {     /* if not masked */
                    buf[nval] = *(planes[j] + i) + scale[j];
                    wbuf[nval++] = wval;
                }
            
            for (n=0; n< N; n++){
                /*fprintf(stderr, "NVAL=%d -----", nval);*/
                buf2[n] = kselect(buf, nval, n );   /* take the n-smallest values */
                mean += buf2[n];
                /* Here, we should get the concerning weights of the n-smallest 
                values in order to compute the weighted mean
                */
                /*fprintf(stderr, "MEAN=%f", mean);*/
            }                
            meanplane[i] = mean/N; /* not weighted mean ! */
            sumwplanes[i] = 1.0;
            /*fprintf(stderr, "I=%d !!", i);*/
        }
    } else {
        for (i = 0; i < nx * ny; i++) {   /* mult. scale to normalize */
            mean = 0;
            for (j = 0; j < np; j++)
                buf[j] = *(planes[j] + i) * scale[j];

            for (n=0; n< N; n++){
                buf2[n] = kselect(buf, nval, n );   /* take the n-smallest values */
                mean += buf2[n];
                /* Here, we should get the concerning weights of the n-smallest 
                values in order to compute the weighted mean
                */    
            }                
            meanplane[i] = mean/N; /* not weighted mean ! */
            sumwplanes[i] = 1.0;

        }
    }
    
    *wplanesum = sumwplanes;
    return meanplane;


}

/* 
 * cube_median_min: find median value of the N-smallest values of the sorted plane 
 * of cube (no weights). It is thought to compute sky background of crowded fields.
 * Implemented by jmiguel@iaa.es on 2011-Oct-25 
 */

extern float *
cube_median_min(float *planes[MAXNPLANES], int np, int nx, int ny, float *scale, 
            int offset, int N)
{
    int i, j;
    static float buf[MAXNPLANES];      /* values of a pixel in all planes */
    static float buf2[MAXNPLANES];     /* values of a pixel in minimal planes */
    float *medianplane;
    int n=0, nvalues=0;
    
    medianplane = (float *) emalloc(nx * ny * sizeof(float));

    if (offset) {                          /* zero offset to normalize, it is much more "safe" than a multiplicative normalization  */
        for (i = 0; i < nx * ny; i++) {    /* median combine cube to plane */
            for (j = 0; j < np; j++)
                buf[j] = *(planes[j] + i) + scale[j];

            nvalues = 0;
            for (n=0; n< N; n++)          /* take the n-smallest values */
                buf2[nvalues++] = kselect(buf, np, n );
                            
            medianplane[i] = median(buf2, nvalues);
        }
    } else {
        for (i = 0; i < nx * ny; i++) {   /* mult. scale to normalize */
            for (j = 0; j < np; j++)
                buf[j] = *(planes[j] + i) * scale[j];

            nvalues = 0;
            for (n=0; n< N; n++)          /* take the n-smallest values */
                 buf2[nvalues++] = kselect(buf, np, n );   
                            
            medianplane[i] = median(buf2, nvalues);

        }
    }

    return medianplane;
}


