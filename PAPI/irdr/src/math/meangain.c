/* meangain.c -- compute weighted clipped means and other statistics */

/*
 * If there are MINMAD or more values to be averaged then the stdev is 
 * estimated using the median absolute deviation, otherwise Poisson noise
 * assumption and gain value are used to estimate the noise level
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "eprintf.h"
#include "median.h"
#include "mean.h"

#define BLANK 0.0
#define MINMAD 5          /* need MINMAD or more points for medabsdev calc. */

static float getsigma(float *arr, float med, int n);

/* calculate robust standard deviation */
static float getsigma(float *arr, float med, int n)
{
    char *gainstr;
    static float gain = 0.0;
    float sig;

    if (!gain) {
        if ((gainstr = getenv("GAIN")) == NULL)
            eprintf("getsigma: ERR: GAIN environment variable not set\n");
        else {
            gain = (float)atof(gainstr);
            printf("GAIN = %f\n", gain);
        }
    
        if (gain <= 0)
            eprintf("getsigma: ERR: gain = %f\n", gain);
    }

    if (n < MINMAD) {
        sig = sqrt(fabs(med) / gain);       /* incorrect if previous coadds */
    } else {
        sig = 1.48 * median_absdev(arr, med, n);
        sig = (sig <= 0.0) ? sqrt(fabs(med) / gain) : sig;
    }

    return sig;
}

/* calculate robust mean of array (no weights) */
extern float mean_nw(float *arr, int n, float nsig)
{
    int i, count = 0;
    float a, med, sig, lcut, hcut, sum = 0.0;

    if (n < 1)
        return BLANK;

    med = median(arr, n);
    sig = getsigma(arr, med, n);

    lcut = med - nsig * sig;
    hcut = med + nsig * sig;

    for (i = 0; i < n; i++) {
        if ((a = arr[i]) > lcut && a < hcut) {
            sum += a;
            count++;
        }
    }

    return (count > 0) ? (sum / count) : BLANK;
}

/* calculate robust mean of array (using weights w) */
extern float mean(float *arr, float *w, int n, float nsig, float *sumw)
{
    int i;
    float a, med, sig, lcut, hcut, sum = 0.0, wsum = 0.0;

    if (n < 1) {
        *sumw = 0.0;
        return BLANK;
    }

    med = median(arr, n);                                /* clipping */
    sig = getsigma(arr, med, n);

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

/* calculate the standard deviation of array arr */
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
