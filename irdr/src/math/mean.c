/* mean.c -- compute weighted clipped means and other statistics */

/*
 * If there are 5 or more values to be averaged then the stdev is estimated
 * using the median absolute deviation and bad points are rejected
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "median.h"
#include "mean.h"

#define BLANK 0.0
#define MINCLIP 5               /* need MINCLIP or more points for clipping */

/* mean_nw: calculate robust mean of array (no weights) */
extern float mean_nw(float *arr, int n, float nsig)
{
    int i, count = 0;
    float a, med, sig, lcut, hcut, sum = 0.0; 

    if (n < MINCLIP) {                                 /* no clipping */
        for (i = 0; i < n; i++)
            sum += arr[i];

        return (n > 0) ? (sum / n) : BLANK;
    }
        
    med = median(arr, n);                              /* clipping */
    sig = median_absdev(arr, med, n) / 0.6745;

    lcut = med - nsig * sig;
    hcut = med + nsig * sig;

    for (i = 0; i < n; i++) {
        if ((a = arr[i]) >= lcut && a <= hcut) {
            sum += a;
            count++;
        }
    }
    /*printf ("\nDEBUG- SUM=%f, COUNT=%d LCUT=%f HCUT=%f", sum, count, lcut, hcut);*/
    return (count > 0) ? (sum / count) : BLANK;
}

/* mean: calculate robust mean of array (using weights w) */
extern float mean(float *arr, float *w, int n, float nsig, float *sumw)
{
    int i;
    float a, med, sig, lcut, hcut, sum = 0.0, wsum = 0.0;

    if (n < MINCLIP) {                                   /* no clipping */
        for (i = 0; i < n; i++) {
            sum += w[i] * arr[i];
            wsum += w[i];
        }

        *sumw = wsum;

        return (wsum > 0.0) ? (sum / wsum) : BLANK;
    }

    med = median(arr, n);                                /* clipping */
    sig = median_absdev(arr, med, n) / 0.6745;

    lcut = med - nsig * sig;
    hcut = med + nsig * sig;

    for (i = 0; i < n; i++) {
        if ((a = arr[i]) > lcut && a < hcut) {     /* reject deviant points */
            sum += w[i] * a;
            wsum += w[i];
        }
    }
 
    *sumw = wsum;

    return (wsum > 0.0) ? (sum / wsum) : BLANK;
}

/* stdev: calculate the standard deviation of array arr */
extern float stdev(float *arr, int n)
{
    int i;
    float mean = 0.0, sig = 0.0;

    if (n < 2)
        return 0.0;

    for (i = 0; i < n; i++) {
        mean += arr[i];
        sig += arr[i] * arr[i];
    }

    mean /= n;

    sig = sig / n - mean * mean;
    
    return (sig > 0.0) ? sqrt(sig) : 0.0;
}
