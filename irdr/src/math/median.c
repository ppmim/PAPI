/* median.c -- median statistics */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "eprintf.h"
#include "kselect.h"
#include "median.h"

#define IS_EVEN(n) (!(n&1))               /* 1 if n is even, else 0 */

/* median: return median value of float array of n elements */
extern float median(const float *arr, int n)
{
    float *buf, median;

    buf = (float *) emalloc(n * sizeof(float));

    memcpy(buf, arr, n * sizeof(float)); /* make copy, kselect reorders data */

    median = kselect(buf, n, n/2 - IS_EVEN(n));

    free(buf);

    return median;
}

/* median_double: return median value of double precision array of n elements */
extern float median_double(const double *arr, int n)
{
    double *buf, median;

    buf = (double *) emalloc(n * sizeof(double));

    memcpy(buf, arr, n * sizeof(double)); /* make copy, kselect reorders data */

    median = kselect_double(buf, n, n/2 - IS_EVEN(n));

    free(buf);

    return median;
}

/* sample_median: sample the image to quickly estimate the median and stddev */
extern float sample_median(const float *img, float *stddev, int n)
{
    int i, count = 0, step = 5;
    float *buf, med;

    buf = (float*) emalloc((n/step + 1) * sizeof(float));
    
    for (i = 0; i < n; i += step)
        buf[count++] = (float)img[i];

    med = kselect(buf, count, count/2 - IS_EVEN(count));

    if (stddev != NULL)
        *stddev = median_absdev(buf, med, count) / 0.6745;   /* for gaussian */

    free(buf);

    return med;
}

/* median_absdev: calculate med. abs. dev. of array arr from signal level m */
extern float median_absdev(const float *arr, float m, int n)
{
    int i;
    float *absdevs, mabsdev;

    if (n <= 1)
        return 0.0;

    absdevs = (float *) emalloc(n * sizeof(float));

    for (i = 0; i < n; i++)
        absdevs[i] = fabs((double)arr[i] - m);

    mabsdev = kselect(absdevs, n, n/2);         /* want larger of middle two */
                                                   /* don't need IS_EVEN() */
    free(absdevs);

    return mabsdev;
}

/* median_absdev_double: median_absdev for double precision array */
extern double median_absdev_double(const double *arr, float m, int n)
{
    int i;
    double *absdevs, mabsdev;

    if (n <= 1)
        return 0.0;

    absdevs = (double *) emalloc(n * sizeof(double));

    for (i = 0; i < n; i++)
        absdevs[i] = fabs((double)arr[i] - m);

    mabsdev = kselect_double(absdevs, n, n/2);
                                                   
    free(absdevs);

    return mabsdev;
}
